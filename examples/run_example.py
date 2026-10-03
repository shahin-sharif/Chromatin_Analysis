#!/usr/bin/env python3
"""Run the synthetic ChromHMM comparison and verify known results."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / 'data'

def matrix(path):
    with path.open() as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    return [[int(value) for value in row[1:]] for row in rows[1:]]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True, type=Path, help='New output directory')
    args = parser.parse_args()
    out = args.outdir.resolve()
    if out.exists():
        parser.error('Output directory already exists; choose a new directory.')
    command = [sys.executable, str(ROOT / 'ChromHMMTools.py'), 'run',
               '--wt', str(DATA / 'wt_states.bed'), '--mt', str(DATA / 'mt_states.bed'),
               '--state-mode', 'shared', '--anno-bed', str(DATA / 'annotations.bed'),
               '--write-change-bed', '--outdir', str(out)]
    env = dict(os.environ, MPLBACKEND='Agg')
    result = subprocess.run(command, env=env)
    if result.returncode:
        return result.returncode
    try:
        if matrix(out / 'overlap.bp.tsv') != [[200, 200], [0, 400]]:
            raise ValueError('Unexpected overlap counts')
        report = (out / 'report.html').read_text()
        if '75.00%' not in report or 'data:image/png;base64,' not in report:
            raise ValueError('Retention or embedded plot missing')
        meta = json.loads((out / 'RUNINFO.json').read_text())
        totals = {label: sum(sum(row) for row in matrix(out / filename))
                  for label, filename in meta['annotation_files'].items()}
        if totals != {'genic': 300, 'promoter': 200, 'intergenic': 300}:
            raise ValueError('Unexpected annotation allocation')
        if (out / 'state_changes.bed').read_text() != 'chr1\t200\t400\tE1->E2\t0\t.\n':
            raise ValueError('Unexpected changed interval')
    except (OSError, ValueError, KeyError) as error:
        print('FAIL:', error, file=sys.stderr)
        return 1
    (out / 'example-validation.json').write_text(json.dumps({
        'status':'PASS', 'shared_bp':800, 'matching_label_percent':75,
        'annotation_bp':totals, 'synthetic':True}, indent=2) + '\n')
    print('PASS: exact overlaps, annotation allocation, changed interval and embedded report.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
