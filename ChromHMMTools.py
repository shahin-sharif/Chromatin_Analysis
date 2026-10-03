#!/usr/bin/env python3
"""Compare ChromHMM segmentations and create a self-contained HTML report."""

import argparse
import base64
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import html
import io
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import csv
import re
import numpy as np

VERSION = "1.0.0"


def open_text(path):
    return (gzip.open if str(path).endswith(".gz") else open)(
        path, "rt", encoding="utf-8"
    )


def natural(s):
    return [int(v) if v.isdigit() else v for v in re.split(r"(\d+)", s)]


def read_bed(path, labeled=True):
    out = defaultdict(list)
    with open_text(path) as f:
        for no, line in enumerate(f, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            p = line.rstrip("\r\n").split("\t")
            if len(p) < (4 if labeled else 3):
                raise ValueError(f"{path} line {no}: missing BED columns")
            s, e = int(p[1]), int(p[2])
            if s < 0 or e <= s or not p[0] or (labeled and not p[3]):
                raise ValueError(f"{path} line {no}: invalid interval or label")
            out[p[0]].append((s, e, p[3]) if labeled else (s, e))
    return {c: sorted(rows) for c, rows in out.items()}


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(e, out[-1][1]))
        else:
            out.append((s, e))
    return out


def segmentation(path):
    result = read_bed(path)
    out = {}
    for chrom, rows in result.items():
        clean = []
        for s, e, state in rows:
            if clean and s < clean[-1][1]:
                raise ValueError(f"{path}: overlapping segmentation on {chrom}:{s}-{e}")
            if clean and s == clean[-1][1] and state == clean[-1][2]:
                clean[-1] = (clean[-1][0], e, state)
            else:
                clean.append((s, e, state))
        out[chrom] = clean
    if not out:
        raise ValueError(f"empty segmentation: {path}")
    return out


def clip_regions(data, include=None, exclude=None):
    from bisect import bisect_right

    result = {}
    for chrom, rows in data.items():
        inc = None if include is None else merge(include.get(chrom, []))
        exc = merge((exclude or {}).get(chrom, []))
        incends = [] if inc is None else [e for _, e in inc]
        excends = [e for _, e in exc]
        kept = []
        for s, e, label in rows:
            pieces = [(s, e)] if inc is None else []
            if inc is not None:
                j = bisect_right(incends, s)
                while j < len(inc) and inc[j][0] < e:
                    pieces.append((max(s, inc[j][0]), min(e, inc[j][1])))
                    j += 1
            for left, right in pieces:
                cur = left
                j = bisect_right(excends, left)
                while j < len(exc) and exc[j][0] < right:
                    a, b = exc[j]
                    if a > cur:
                        kept.append((cur, min(a, right), label))
                    cur = max(cur, b)
                    j += 1
                if cur < right:
                    kept.append((cur, right, label))
        if kept:
            result[chrom] = kept
    return result


def transitions(data, states, bin_size):
    index = {s: i for i, s in enumerate(states)}
    c = np.zeros((len(states), len(states)), dtype=np.int64)
    used = 0
    for rows in data.values():
        previous = None
        for s, e, label in rows:
            first = (s + bin_size - 1) // bin_size
            stop = e // bin_size
            n = stop - first
            if n <= 0:
                continue
            used += n
            i = index[label]
            c[i, i] += n - 1
            if previous and previous[0] == first:
                c[previous[1], i] += 1
            previous = (stop, i)
    return c, used


def overlap_tiles(wt, mt):
    for chrom in sorted(set(wt) & set(mt)):
        a = wt[chrom]
        b = mt[chrom]
        i = j = 0
        while i < len(a) and j < len(b):
            s = max(a[i][0], b[j][0])
            e = min(a[i][1], b[j][1])
            if e > s:
                yield chrom, s, e, a[i][2], b[j][2]
            if a[i][1] < b[j][1]:
                i += 1
            elif b[j][1] < a[i][1]:
                j += 1
            else:
                i += 1
                j += 1


def normalize(c, axis=1):
    sums = c.sum(axis=axis, keepdims=True)
    return np.divide(c, sums, out=np.full(c.shape, np.nan, dtype=float), where=sums > 0)


def write_matrix(path, m, rows, cols):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["from", *cols])
        for label, row in zip(rows, m):
            w.writerow([label, *["NA" if np.isnan(x) else str(x) for x in row]])


def read_matrix(path):
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        cols = header[1:]
        rows = []
        values = []
        for r in reader:
            if len(r) != len(header):
                raise ValueError(f"malformed matrix row: {path}")
            rows.append(r[0])
            values.append([float(x) if x != "NA" else np.nan for x in r[1:]])
    if not cols or len(set(cols)) != len(cols) or len(set(rows)) != len(rows):
        raise ValueError("empty or duplicate state labels")
    m = np.asarray(values, dtype=float).reshape(len(rows), len(cols))
    if not np.isfinite(m).all() or (m < 0).any():
        raise ValueError(f"invalid count matrix: {path}")
    return m, rows, cols


def annotations(a):
    if a.anno_bed:
        return read_bed(a.anno_bed)
    if not a.gtf:
        return None
    genes = defaultdict(list)
    with open_text(a.gtf) as f:
        for no, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                raise ValueError(f"GTF line {no}: expected 9 fields")
            if p[2] != "gene":
                continue
            s, e = int(p[3]) - 1, int(p[4])
            if s < 0 or e <= s or p[6] not in ("+", "-"):
                raise ValueError(f"invalid GTF gene line {no}")
            genes[p[0]].append((s, e, "genic"))
            tss = s if p[6] == "+" else e - 1
            left, right = (
                (tss - a.promoter_up, tss + a.promoter_down + 1)
                if p[6] == "+"
                else (tss - a.promoter_down, tss + a.promoter_up + 1)
            )
            genes[p[0]].append((max(0, left), right, "promoter"))
    if not genes:
        raise ValueError(
            "GTF contains no gene features; use annotation BED for other feature definitions"
        )
    return genes


def annotation_partition(data, priority):
    result = {}
    for chrom, rows in data.items():
        events = defaultdict(list)
        for s, e, label in rows:
            events[s].append((label, 1))
            events[e].append((label, -1))
        points = sorted(events)
        active = Counter()
        out = []
        for i, s in enumerate(points[:-1]):
            for label, delta in events[s]:
                active[label] += delta
            present = [label for label, n in active.items() if n > 0]
            if present:
                label = min(
                    present,
                    key=lambda x: (
                        priority.index(x) if x in priority else len(priority),
                        x,
                    ),
                )
                e = points[i + 1]
                if out and out[-1][1] == s and out[-1][2] == label:
                    out[-1] = (out[-1][0], e, label)
                else:
                    out.append((s, e, label))
        result[chrom] = out
    return result


def split_annotation(chrom, s, e, partition):
    from bisect import bisect_right

    rows, ends = partition.get(chrom, ([], []))
    j = bisect_right(ends, s)
    cur = s
    while j < len(rows) and rows[j][0] < e:
        a, b, label = rows[j]
        if a > cur:
            yield cur, min(e, a), "intergenic"
        left, right = max(cur, a), min(e, b)
        if right > left:
            yield left, right, label
        cur = max(cur, right)
        j += 1
    if cur < e:
        yield cur, e, "intergenic"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def compare(a, out):
    if a.bin_size <= 0 or a.promoter_up < 0 or a.promoter_down < 0:
        raise ValueError("invalid bin/promoter sizes")
    wt = segmentation(a.wt)
    mt = segmentation(a.mt)
    if a.state_mode == "mapped":
        if not a.state_map:
            raise ValueError(
                "mapped mode requires --state-map (MT state TAB matching WT state)"
            )
        mapping = {}
        with open_text(a.state_map) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) != 2 or parts[0] in mapping:
                    raise ValueError("invalid/duplicate state-map entry")
                mapping[parts[0]] = parts[1]
        if len(set(mapping.values())) != len(mapping):
            raise ValueError("state mapping must be one-to-one")
        source = {r[2] for v in mt.values() for r in v}
        target = {r[2] for v in wt.values() for r in v}
        if set(mapping) != source or not set(mapping.values()) <= target:
            raise ValueError(
                "state map must cover exactly all MT states and map to observed WT states"
            )
        mt = {
            c: [(s, e, mapping[label]) for s, e, label in rows]
            for c, rows in mt.items()
        }
    elif a.state_map:
        raise ValueError("--state-map requires mapped mode")
    inc = read_bed(a.include_bed, False) if a.include_bed else None
    exc = read_bed(a.exclude_bed, False) if a.exclude_bed else None
    if a.autosomes_only:
        keep = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}
        wt = {c: r for c, r in wt.items() if c in keep}
        mt = {c: r for c, r in mt.items() if c in keep}
    wt = clip_regions(wt, inc, exc)
    mt = clip_regions(mt, inc, exc)
    ws = sorted({r[2] for v in wt.values() for r in v}, key=natural)
    ms = sorted({r[2] for v in mt.values() for r in v}, key=natural)
    if not ws or not ms:
        raise ValueError("region selection leaves an empty condition")
    if a.state_mode != "distinct":
        ws = ms = sorted(set(ws) | set(ms), key=natural)
    wi = {s: i for i, s in enumerate(ws)}
    mi = {s: i for i, s in enumerate(ms)}
    matrix = np.zeros((len(ws), len(ms)), dtype=np.int64)
    anno = annotations(a)
    part = (
        annotation_partition(anno, a.anno_priority.split(","))
        if anno is not None
        else None
    )
    partition = {c: (r, [x[1] for x in r]) for c, r in (part or {}).items()}
    by_annotation = {}
    changes = []
    for chrom, s, e, w, m in overlap_tiles(wt, mt):
        matrix[wi[w], mi[m]] += e - s
        if a.write_change_bed and (a.state_mode == "distinct" or w != m):
            changes.append((chrom, s, e, w, m))
        if part is not None:
            for left, right, label in split_annotation(chrom, s, e, partition):
                if label not in by_annotation:
                    by_annotation[label] = np.zeros_like(matrix)
                by_annotation[label][wi[w], mi[m]] += right - left
    if matrix.sum() == 0:
        raise ValueError(
            "conditions have no shared covered territory; check chromosome naming and masks"
        )
    write_matrix(out / "overlap.bp.tsv", matrix, ws, ms)
    write_matrix(out / "overlap.row.tsv", normalize(matrix), ws, ms)
    write_matrix(out / "overlap.col.tsv", normalize(matrix, 0), ws, ms)
    write_matrix(out / "overlap.global.tsv", matrix / matrix.sum(), ws, ms)
    annotation_files = {}
    for n, (label, m) in enumerate(sorted(by_annotation.items()), 1):
        filename = f"annotation_{n}.overlap.bp.tsv"
        write_matrix(out / filename, m, ws, ms)
        annotation_files[label] = filename
    transition_meta = {}
    for name, data, states in [("WT", wt, ws), ("MT", mt, ms)]:
        counts, bins = transitions(data, states, a.bin_size)
        write_matrix(out / f"{name}.trans.counts.tsv", counts, states, states)
        write_matrix(out / f"{name}.trans.probs.tsv", normalize(counts), states, states)
        transition_meta[name] = {
            "full_bins": bins,
            "transitions": int(counts.sum()),
            "covered_bp": sum(e - s for rows in data.values() for s, e, _ in rows),
            "bin_bp": a.bin_size,
        }
    if a.write_change_bed:
        # Merge adjacent tiles with identical label pairs, then apply minimum length.
        merged = []
        for c, s, e, w, m in changes:
            if (
                merged
                and merged[-1][0] == c
                and merged[-1][2] == s
                and merged[-1][3:] == (w, m)
            ):
                merged[-1] = (c, merged[-1][1], e, w, m)
            else:
                merged.append((c, s, e, w, m))
        with open(out / "state_changes.bed", "w") as f:
            for c, s, e, w, m in merged:
                if e - s >= a.min_change_len:
                    f.write(f"{c}\t{s}\t{e}\t{w}->{m}\t0\t.\n")
    inputs = {
        name: {"path": str(path), "sha256": sha256(path)}
        for name, path in [
            ("wt", a.wt),
            ("mt", a.mt),
            ("include", a.include_bed),
            ("exclude", a.exclude_bed),
            ("gtf", a.gtf),
            ("annotation", a.anno_bed),
            ("mapping", a.state_map),
        ]
        if path
    }
    info = {
        "version": VERSION,
        "time_utc": datetime.now(timezone.utc).isoformat(),
        "state_mode": a.state_mode,
        "settings": vars(a),
        "inputs": inputs,
        "shared_bp": int(matrix.sum()),
        "conditions": transition_meta,
        "annotation_files": annotation_files,
    }
    for name in ("WT", "MT"):
        info["conditions"][name]["unshared_bp"] = transition_meta[name][
            "covered_bp"
        ] - int(matrix.sum())
    (out / "RUNINFO.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


def heatmap_image(m, rows, cols, title, signed=False):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(max(5, 0.45 * len(cols)), max(4, 0.4 * len(rows))))
    try:
        options = {"cmap": "coolwarm" if signed else "viridis"}
        if signed:
            finite = np.abs(m[np.isfinite(m)])
            lim = max(float(finite.max()) if finite.size else 1.0, 1e-9)
            options.update(vmin=-lim, vmax=lim)
        im = ax.imshow(np.ma.masked_invalid(m), aspect="auto", **options)
        ax.set_xticks(range(len(cols)), cols, rotation=90)
        ax.set_yticks(range(len(rows)), rows)
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=120)
        return (
            '<img style="max-width:100%" alt="'
            + html.escape(title, quote=True)
            + '" src="data:image/png;base64,'
            + base64.b64encode(buffer.getvalue()).decode()
            + '">'
        )
    finally:
        plt.close(fig)


def table_html(headers, rows):
    return (
        "<table><tr>"
        + "".join("<th>" + html.escape(str(x)) + "</th>" for x in headers)
        + "</tr>"
        + "".join(
            "<tr>"
            + "".join("<td>" + html.escape(str(x)) + "</td>" for x in r)
            + "</tr>"
            for r in rows
        )
        + "</table>"
    )


def report(indir, out, a):
    info = json.loads((indir / "RUNINFO.json").read_text())
    c, rows, cols = read_matrix(indir / "overlap.bp.tsv")
    if info["state_mode"] not in ("shared", "mapped", "distinct"):
        raise ValueError("unknown state mode")
    total = c.sum()
    if total <= 0:
        raise ValueError("no overlap counts for report")
    p = normalize(c)
    reverse = normalize(c, 0)
    expected = c.sum(1, keepdims=True) * c.sum(0, keepdims=True) / total
    enrich = np.log2((c + 1) / (expected + 1))
    shared = info["state_mode"] in ("shared", "mapped")
    conservation = (
        sum(c[i, cols.index(s)] for i, s in enumerate(rows) if s in cols) / total
        if shared
        else None
    )
    summary = []
    for i, state in enumerate(rows):
        if c[i].sum() > 0:
            nonzero = p[i][p[i] > 0]
            entropy = float(-np.sum(nonzero * np.log2(nonzero)))
            targets = 2**entropy
        else:
            entropy = targets = None
        retain = (
            float(p[i, cols.index(state)])
            if shared and state in cols and c[i].sum() > 0
            else None
        )
        summary.append([state, int(c[i].sum()), retain, entropy, targets])
    with open(out / "state_summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["WT_state", "shared_bp", "retention", "entropy_bits", "effective_targets"]
        )
        w.writerows([["NA" if x is None else x for x in r] for r in summary])
    pairs = []
    for i, w in enumerate(rows):
        for j, m in enumerate(cols):
            if shared and w == m:
                continue
            if (
                c[i, j] > 0
                and c[i, j] >= a.min_bp
                and np.isfinite(p[i, j])
                and p[i, j] >= a.min_fraction
            ):
                pairs.append(
                    [
                        w,
                        m,
                        int(c[i, j]),
                        float(p[i, j]),
                        float(reverse[i, j]),
                        float(enrich[i, j]),
                    ]
                )
    pairs.sort(key=lambda r: (-r[2], -r[3], r[0], r[1]))
    pairs = pairs[: a.top_k]
    headers = [
        "WT_state",
        "MT_state",
        "bp",
        "P_MT_given_WT",
        "P_WT_given_MT",
        "log2_enrichment",
    ]
    with open(out / "top_pairs.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(headers)
        w.writerows(pairs)
    content = [
        '<!doctype html><html><head><meta charset="utf-8"><title>ChromHMM comparison</title><style>body{font:16px sans-serif;max-width:1100px;margin:32px auto;padding:16px}td,th{border:1px solid #ddd;padding:6px}table{border-collapse:collapse}figure{margin:20px 0}</style></head><body>',
        "<h1>ChromHMM comparison</h1>",
        f'<p>Shared covered territory: {int(total):,} bp. State mode: {html.escape(info["state_mode"])}.</p>',
    ]
    if shared:
        content.append(
            f"<p>Matching-label retention: {conservation:.2%}. This interpretation assumes the declared state correspondence is valid.</p>"
        )
    else:
        content.append(
            "<p>States were declared distinct. Matching numbers are not treated as conservation; all state pairs are eligible for ranking.</p>"
        )
    content.append(
        "<p>Fractions refer to shared covered territory after masks, not the whole genome. Missing-support rows are NA. Rankings are descriptive (descending bp, then conditional probability), not significance tests. Enrichment uses a one-bp pseudocount.</p>"
    )
    if not a.no_plots:
        for m, title, sign in [
            (np.log10(c + 1), "Overlap bp (log10 + 1)", False),
            (p, "P(MT | WT)", False),
            (reverse, "P(WT | MT)", False),
            (enrich, "log2 enrichment over independence", True),
        ]:
            content.append(heatmap_image(m, rows, cols, title, sign))
    content += [
        "<h2>State summaries</h2>",
        table_html(
            [
                "WT state",
                "Shared bp",
                "Retention",
                "Entropy (bits)",
                "Effective targets",
            ],
            [["NA" if x is None else x for x in r] for r in summary],
        ),
        "<h2>Top state pairs</h2>",
        table_html(headers, pairs),
    ]
    transition_summary = []
    for name in ("WT", "MT"):
        t, rs, cs = read_matrix(indir / f"{name}.trans.counts.tsv")
        if set(rs) != set(cs):
            raise ValueError("transition row/column state sets differ")
        t = t[:, [cs.index(x) for x in rs]]
        tp = normalize(t)
        for i, state in enumerate(rs):
            support = int(t[i].sum())
            stay = float(tp[i, i]) if support else None
            transition_summary.append(
                [name, state, support, "NA" if stay is None else stay]
            )
        content.append(
            f"<h2>{name} spatial adjacency</h2><p>Only complete bins aligned to coordinate zero count; transitions never cross gaps. These are empirical segmentation transitions, not learned HMM parameters.</p>"
        )
        if not a.no_plots:
            content.append(heatmap_image(tp, rs, rs, name + " adjacency probabilities"))
    with open(out / "transition_summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            ["condition", "state", "outgoing_bin_pairs", "self_transition_probability"]
        )
        w.writerows(transition_summary)
    content += [
        table_html(
            ["Condition", "State", "Outgoing bin pairs", "Self-transition probability"],
            transition_summary,
        )
    ]
    for label, filename in info.get("annotation_files", {}).items():
        if Path(filename).name != filename:
            raise ValueError("invalid annotation matrix filename")
        m, rs, cs = read_matrix(indir / filename)
        content.append(
            "<h2>Annotation: "
            + html.escape(label)
            + "</h2><p>Assigned overlap: "
            + str(int(m.sum()))
            + " bp. Annotation boundaries split tiles; declared priority resolves overlapping classes.</p>"
        )
        if not a.no_plots:
            content.append(
                heatmap_image(
                    np.log10(m + 1), rs, cs, "Annotation " + label + " (log10 bp + 1)"
                )
            )
    content += [
        "<h2>Run metadata</h2><pre>"
        + html.escape(json.dumps(info, indent=2))
        + "</pre></body></html>"
    ]
    (out / "report.html").write_text("\n".join(content), encoding="utf-8")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--version", action="version", version=VERSION)
    subs = p.add_subparsers(dest="command", required=True)
    for name in ("compare", "report", "run"):
        s = subs.add_parser(name)
        s.add_argument("--outdir", required=True)
        s.add_argument("--force", action="store_true")
        if name in ("compare", "run"):
            s.add_argument("--wt", required=True)
            s.add_argument("--mt", required=True)
            s.add_argument(
                "--state-mode",
                choices=["shared", "mapped", "distinct"],
                required=True,
                help="Declare state semantics; use distinct when correspondence is unknown",
            )
            s.add_argument("--state-map")
            s.add_argument("--bin-size", type=int, default=200)
            s.add_argument("--include-bed")
            s.add_argument("--exclude-bed")
            s.add_argument("--autosomes-only", action="store_true")
            g = s.add_mutually_exclusive_group()
            g.add_argument("--gtf")
            g.add_argument("--anno-bed")
            s.add_argument("--promoter-up", type=int, default=1000)
            s.add_argument("--promoter-down", type=int, default=200)
            s.add_argument("--anno-priority", default="promoter,genic,intergenic")
            s.add_argument("--write-change-bed", action="store_true")
            s.add_argument("--min-change-len", type=int, default=1)
        if name in ("report", "run"):
            if name == "report":
                s.add_argument("--indir", required=True)
            s.add_argument("--no-plots", action="store_true")
            s.add_argument("--min-bp", type=int, default=0)
            s.add_argument("--min-fraction", type=float, default=0.0)
            s.add_argument("--top-k", type=int, default=20)
    return p


def main(argv=None):
    a = parser().parse_args(argv)
    try:
        if a.command in ("report", "run") and (
            a.min_bp < 0 or not 0 <= a.min_fraction <= 1 or a.top_k < 1
        ):
            raise ValueError("invalid report thresholds")
        if a.command != "report" and a.min_change_len < 1:
            raise ValueError("minimum change length must be positive")
        target = Path(a.outdir)
        target.mkdir(parents=True, exist_ok=True)
        protected = []
        if a.command == "report":
            source_dir = Path(a.indir)
            metadata = json.loads((source_dir / "RUNINFO.json").read_text())
            protected = [
                source_dir / x
                for x in [
                    "RUNINFO.json",
                    "overlap.bp.tsv",
                    "WT.trans.counts.tsv",
                    "MT.trans.counts.tsv",
                    *metadata.get("annotation_files", {}).values(),
                ]
            ]
        else:
            protected = [
                Path(x)
                for x in (
                    a.wt,
                    a.mt,
                    a.include_bed,
                    a.exclude_bed,
                    a.gtf,
                    a.anno_bed,
                    a.state_map,
                )
                if x
            ]
        with tempfile.TemporaryDirectory(prefix=".chromhmm-", dir=target) as tmp:
            stage = Path(tmp)
            if a.command in ("compare", "run"):
                compare(a, stage)
            if a.command in ("report", "run"):
                report(stage if a.command == "run" else Path(a.indir), stage, a)
            for source in stage.iterdir():
                dest = target / source.name
                if any(
                    dest.resolve() == x.resolve()
                    or (dest.exists() and x.exists() and os.path.samefile(dest, x))
                    for x in protected
                ):
                    raise ValueError(f"output aliases an input: {dest}")
                if dest.exists() and (not a.force or not dest.is_file()):
                    raise ValueError(f"output exists: {dest}; use --force")
            for source in stage.iterdir():
                os.replace(source, target / source.name)
        print(f"Wrote {a.command} outputs to {target}", file=sys.stderr)
        return 0
    except (ValueError, OSError, ImportError, KeyError, StopIteration) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
