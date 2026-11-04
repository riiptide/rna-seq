#!/usr/bin/env bash
# Run from the repository root with Nextflow, Python and seqtk 1.4 on PATH.
set -euo pipefail
mkdir -p work
nextflow config -o json > work/test-config.json
python3 - <<'PY'
import csv, gzip, json, os, subprocess
from pathlib import Path
p = json.loads(Path('work/test-config.json').read_text())['params']
output = Path(p['test_dir'])
output.mkdir(parents=True, exist_ok=True)
rows = []
with open(p['samplesheet']) as handle:
    for number, row in enumerate(csv.DictReader(handle), start=1):
        lane = row.get('lane') or f'L{number:03d}'
        pair = []
        for mate in [1, 2]:
            dest = output / f"{row['sample']}_{lane}_R{mate}.fastq.gz"
            temporary = dest.with_suffix('.tmp.gz')
            # The same seed selects the same record positions in synchronized mates.
            with open(temporary, 'wb') as handle_out:
                sampled = subprocess.Popen(['seqtk', 'sample', f"-s{p['seed']}", row[f'fastq_{mate}'],
                                            str(p['test_pairs_per_lane'])], stdout=subprocess.PIPE)
                zipped = subprocess.run(['gzip', '-n'], stdin=sampled.stdout, stdout=handle_out)
                sampled.stdout.close()
                if sampled.wait() or zipped.returncode:
                    raise RuntimeError('seqtk/gzip failed')
            os.replace(temporary, dest)
            pair.append(str(dest))
        subprocess.run(['python3', 'bin/validate_fastq_pairs.py', *pair, '--expected-pairs', str(p['test_pairs_per_lane'])], check=True)
        rows.append([row['sample'], *pair])
with open('tests/samplesheet_test.csv', 'w') as handle:
    writer = csv.writer(handle)
    writer.writerow(['sample', 'fastq_1', 'fastq_2'])
    writer.writerows(rows)
PY
