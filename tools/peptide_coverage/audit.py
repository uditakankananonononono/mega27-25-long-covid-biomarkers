#!/usr/bin/env python3
"""Exact reference coverage audit. No assay values, detection or safety claims."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

AA = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$')
REF = re.compile(r'^[A-Z*]+$')
LIMITS = [
    'Reference coverage is not peptide detection, presentation, biological absence or safety.',
    'A hash-verified file is the expected file, not proof that every experimental search used it.',
    'Exact matching does not collapse I/L, remove modifications or translate nucleotide sequences.',
    'HLA, donor, cell model, condition, independence and comparator eligibility are not evaluated.',
    'Unknown residues and sequence breaks cannot support a full absence conclusion.',
]

def records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        header, sequence = None, []
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    if not sequence:
                        raise ValueError(f'Empty sequence at {header}')
                    yield header, ''.join(sequence)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f'Empty header at line {line_number}')
                sequence = []
            else:
                if header is None:
                    raise ValueError(f'Sequence before first header at line {line_number}')
                if not REF.fullmatch(line):
                    raise ValueError(f'Invalid protein sequence at line {line_number}')
                sequence.append(line)
        if header is None or not sequence:
            raise ValueError('Empty FASTA or final record without sequence')
        yield header, ''.join(sequence)


def audit(path, peptides, expected_sha256=None, expected_bytes=None, max_hits=5):
    path = Path(path)
    peptides = list(dict.fromkeys(peptides))
    if not peptides or any(not AA.fullmatch(p) for p in peptides):
        raise ValueError('Supply unmodified uppercase canonical amino-acid peptides')
    if max_hits < 0:
        raise ValueError('max_hits must be nonnegative')
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    sha, size = h.hexdigest(), path.stat().st_size
    if expected_sha256 and sha != expected_sha256.lower():
        raise ValueError('SHA-256 mismatch; no coverage conclusions produced')
    if expected_bytes is not None and size != expected_bytes:
        raise ValueError('Byte-size mismatch; no coverage conclusions produced')
    hits = {p: [] for p in peptides}
    counts = dict.fromkeys(peptides, 0)
    n_records, ambiguous_records = 0, 0
    for header, seq in records(path):
        n_records += 1
        ambiguous_records += int(any(a not in 'ACDEFGHIKLMNPQRSTVWY*' for a in seq))
        # Do not join across stop symbols, ambiguous residues or record boundaries.
        for p in peptides:
            start = 0
            while True:
                position = seq.find(p, start)
                if position < 0:
                    break
                counts[p] += 1
                if len(hits[p]) < max_hits:
                    hits[p].append({'record_index': n_records, 'header': header,
                                    'start_1based': position + 1,
                                    'end_1based_inclusive': position + len(p)})
                start = position + 1
    verified = expected_sha256 is not None
    report = {'schema_version': 1, 'reference': {'filename': path.name,
              'sha256': sha, 'bytes': size, 'records': n_records,
              'ambiguous_records': ambiguous_records,
              'expected_hash_verified': verified}, 'peptides': [], 'limitations': LIMITS}
    for p in peptides:
        status = ('exact_reference_match' if counts[p] else
                  'not_in_verified_reference' if verified and not ambiguous_records else
                  'no_exact_match_in_unverified_or_ambiguous_reference')
        report['peptides'].append({'sequence': p, 'status': status,
             'occurrences': counts[p], 'matches': hits[p],
             'matches_truncated': counts[p] > len(hits[p]),
             'interpretation': 'Reference inclusion only; no assay or clinical inference.'})
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', required=True, type=Path)
    parser.add_argument('--peptide', required=True, action='append')
    parser.add_argument('--expected-sha256')
    parser.add_argument('--expected-bytes', type=int)
    parser.add_argument('--max-hits', type=int, default=5)
    args = parser.parse_args(argv)
    try:
        result = audit(args.fasta, args.peptide, args.expected_sha256,
                       args.expected_bytes, args.max_hits)
    except (ValueError, OSError, UnicodeError, EOFError) as exc:
        print(f'Audit stopped: {exc}', file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 0

if __name__ == '__main__':
    sys.exit(main())
