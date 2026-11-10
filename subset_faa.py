#!/usr/bin/env python3
"""
Pull proteins out of a FASTA using the hits in an hmmsearch --tblout file.

The target names in the tblout are protein IDs (e.g. ctg1711_20), and those are
the same IDs as the FASTA headers, so this is a direct ID lookup.

Examples
--------
# every protein with a hit
python subset_faa.py hits.tbl proteins.faa -o hits.faa

# only strong hits
python subset_faa.py hits.tbl proteins.faa -e 1e-20 -o hits.faa

# only hits to one HMM (the query name column of the tblout)
python subset_faa.py hits.tbl proteins.faa -q MCP_144 -o mcp144.faa

# just the ID list, no FASTA needed
python subset_faa.py hits.tbl --ids-only > ids.txt

# every protein on any contig that had a hit (strips the trailing _N)
python subset_faa.py hits.tbl proteins.faa --contigs -o hit_contigs.faa
"""

import argparse
import re
import sys


def read_tblout(path, max_evalue=None, min_score=None, queries=None):
    """Yield (target, query, evalue, score) for each hit line that passes the filters."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 6:
                continue
            target, query = f[0], f[2]
            try:
                evalue, score = float(f[4]), float(f[5])
            except ValueError:
                continue                       # header or malformed line
            if max_evalue is not None and evalue > max_evalue:
                continue
            if min_score is not None and score < min_score:
                continue
            if queries and query not in queries:
                continue
            yield target, query, evalue, score


def to_contig(protein_id):
    """ctg1711_20 -> ctg1711   (only used with --contigs)"""
    return re.sub(r"_\d+$", "", protein_id)


def iter_fasta(path):
    """Stream (header, [seq_lines]) so a big .faa stays out of memory."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, seq
                header, seq = line.rstrip("\n"), []
            else:
                seq.append(line.rstrip("\n"))
    if header is not None:
        yield header, seq


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tblout", help="hmmsearch --tblout file")
    p.add_argument("fasta", nargs="?", help="protein FASTA to subset")
    p.add_argument("-o", "--out", help="output FASTA (default: stdout)")
    p.add_argument("-e", "--evalue", type=float, help="keep hits with E-value <= this")
    p.add_argument("-s", "--score", type=float, help="keep hits with bit score >= this")
    p.add_argument("-q", "--query", action="append",
                   help="only this HMM; repeat the flag for several")
    p.add_argument("--contigs", action="store_true",
                   help="keep every protein from the hit contigs, not just the hit proteins")
    p.add_argument("--ids-only", action="store_true", help="print the IDs and stop")
    p.add_argument("-v", "--invert", action="store_true",
                   help="keep the proteins that did NOT hit")
    args = p.parse_args()

    if not args.ids_only and not args.fasta:
        p.error("a FASTA file is required unless you pass --ids-only")

    hits = list(read_tblout(args.tblout, args.evalue, args.score,
                            set(args.query) if args.query else None))
    proteins = {t for t, _, _, _ in hits}
    if not proteins:
        sys.exit("No hits passed the filters. Check the thresholds or the -q name.")

    wanted = {to_contig(x) for x in proteins} if args.contigs else set(proteins)
    msg = f"{len(hits)} hit lines -> {len(proteins)} distinct proteins"
    if args.contigs:
        msg += f" on {len(wanted)} contigs"
    print(msg, file=sys.stderr)

    if args.ids_only:
        for x in sorted(wanted):
            print(x)
        return

    out = open(args.out, "w") if args.out else sys.stdout
    seen, emitted, written, dups = set(), set(), 0, 0
    try:
        for header, seq in iter_fasta(args.fasta):
            seq_id = header[1:].split()[0]
            key = to_contig(seq_id) if args.contigs else seq_id
            hit = key in wanted
            if hit:
                seen.add(key)
            if hit != args.invert:
                if seq_id in emitted:      # duplicate header in the input FASTA
                    dups += 1
                    continue
                emitted.add(seq_id)
                out.write(header + "\n")
                out.write("\n".join(seq) + "\n")
                written += 1
    finally:
        if out is not sys.stdout:
            out.close()

    print(f"wrote {written} sequences ({len(emitted)} unique IDs)", file=sys.stderr)
    if dups:
        print(f"NOTE: skipped {dups} repeated header(s) in {args.fasta}; "
              "each ID was written once", file=sys.stderr)
    missing = wanted - seen
    if missing:
        print(f"WARNING: {len(missing)} ID(s) from the table were not found in the FASTA, "
              f"e.g. {sorted(missing)[0]}", file=sys.stderr)


if __name__ == "__main__":
    main()
