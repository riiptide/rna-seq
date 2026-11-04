#!/usr/bin/env python3
"""Validate synchronized, four-line FASTQ files without loading them into memory."""

import argparse
import gzip
import json
from itertools import zip_longest
from pathlib import Path


def records(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        number = 0
        while header := handle.readline():
            number += 1
            sequence = handle.readline().rstrip('\r\n')
            separator = handle.readline()
            quality = handle.readline().rstrip('\r\n')
            if (not header.startswith('@') or not separator.startswith('+')
                    or not sequence or len(sequence) != len(quality)):
                raise ValueError(f'{path}: invalid four-line FASTQ record {number}')
            name = header.split()[0][1:].removesuffix('/1').removesuffix('/2')
            if not name:
                raise ValueError(f'{path}: empty read identifier at record {number}')
            yield name


def validate(read1, read2, expected_pairs=None):
    count = 0
    for count, (left, right) in enumerate(
            zip_longest(records(read1), records(read2)), start=1):
        if left is None or right is None or left != right:
            raise ValueError(f'Pair {count} is not synchronized: {left!r} / {right!r}')
    if count == 0:
        raise ValueError('FASTQ inputs are empty')
    if expected_pairs is not None and count != expected_pairs:
        raise ValueError(f'Expected {expected_pairs} pairs; observed {count}')
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('read1', type=Path)
    parser.add_argument('read2', type=Path)
    parser.add_argument('--expected-pairs', type=int)
    args = parser.parse_args()
    count = validate(args.read1, args.read2, args.expected_pairs)
    print(json.dumps({'read1': str(args.read1), 'read2': str(args.read2),
                      'read_pairs': count, 'synchronized': True}))


if __name__ == '__main__':
    main()
