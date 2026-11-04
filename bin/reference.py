#!/usr/bin/env python3
"""Download the frozen reference if requested, then verify its archive and extracted files."""
import argparse
import hashlib
from pathlib import Path
import tarfile
import urllib.request


def digest(path):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def expected(manifest, suffix):
    matches = [line.split()[0] for line in Path(manifest).read_text().splitlines()
               if line.split()[-1].endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f'Expected one checksum for {suffix}')
    return matches[0]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index', default='index.idx')
    p.add_argument('--t2g', default='t2g.txt')
    p.add_argument('--inputs-manifest', required=True)
    p.add_argument('--downloads-manifest')
    p.add_argument('--url')
    p.add_argument('--retries', type=int, default=5)
    a = p.parse_args()
    if a.url:
        archive = Path('human_index_standard.tar.xz')
        for attempt in range(a.retries):
            try:
                urllib.request.urlretrieve(a.url, archive)
                if digest(archive) != expected(a.downloads_manifest, archive.name):
                    raise ValueError('Reference archive checksum mismatch')
                break
            except Exception:
                if attempt == a.retries-1:
                    raise
        with tarfile.open(archive) as handle:
            for target in [a.index, a.t2g]:
                member = handle.getmember(Path(target).name)
                with handle.extractfile(member) as source, open(target, 'wb') as dest:
                    import shutil
                    shutil.copyfileobj(source, dest)
        archive.unlink()
    for path, suffix in [(a.index, '/index.idx'), (a.t2g, '/t2g.txt')]:
        if digest(path) != expected(a.inputs_manifest, suffix):
            raise ValueError(f'Reference checksum mismatch: {path}')
    print('Verified frozen Ensembl 108 index and t2g checksums')


if __name__ == '__main__':
    main()
