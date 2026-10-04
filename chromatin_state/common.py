import csv
import gzip
import hashlib
import json
import math
from pathlib import Path


def require(value, message):
    if not value:
        raise ValueError(message)


def text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_tsv(path, required=()):
    with text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames)
                and set(required) <= set(reader.fieldnames), f'{path}: missing/duplicate columns; require {required}')
        rows = list(reader)
    require(all(None not in r and None not in r.values() for r in rows), f'{path}: malformed TSV row')
    return rows


def write_tsv(path, rows, columns):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: 'NA' if row.get(k) is None else row[k] for k in columns})


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def number(value, missing=True):
    if missing and str(value).lower() in ('', 'na', 'nan', 'none'):
        return None
    value = float(value)
    require(math.isfinite(value), 'Nonfinite numeric value')
    return value


def keys(obj, allowed, required=()):
    require(isinstance(obj, dict) and not (obj.keys()-set(allowed)) and set(required)<=obj.keys(),
            f'Invalid configuration keys: {obj}; allowed {allowed}, required {required}')


def safe_name(name):
    import re
    require(isinstance(name, str) and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name), 'Unsafe/empty name: '+str(name))
    return name
