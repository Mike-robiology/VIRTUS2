#!/usr/bin/env python3
"""
Filter a master FASTA (e.g., viruses.fasta) based on a list of accessions/names
from a text file (e.g., sam-viruses.txt).

Key features:
- Robustly parses accessions from list lines like "NC 001806.1 Human herpesvirus 1, complete genome"
  by normalizing to NC_001806.1 (underscore instead of space).
- Matches FASTA headers that may contain forms like ">NC_001806.1 ..." or ">ref|NC_001806.1| ...".
- Optional version-insensitive matching (default: on) so NC_001806.1 matches NC_001806.2.

Usage:
  python filter_fasta_by_list.py --fasta path/to/viruses.fasta --list path/to/sam-viruses.txt --out filtered.fasta

Exit code is 0 on success. Writes summary to stderr.
"""

from __future__ import annotations

import argparse
import sys
import re
from typing import Iterable, Set, Tuple, Optional


ACCESSION_RE = re.compile(r"[A-Z]{1,3}_\d+\.\d+")


def parse_targets(list_path: str) -> Tuple[Set[str], Set[str]]:
    """Parse the target accession list file.

    Returns two sets:
    - with_version: e.g., {"NC_001806.1", ...}
    - no_version: e.g., {"NC_001806", ...}
    """
    with_version: Set[str] = set()
    no_version: Set[str] = set()

    # Helper to normalize to underscore form and split version
    def _add(acc: str) -> None:
        acc = acc.strip()
        if not acc:
            return
        # Normalize any accidental space between prefix and numeric part
        acc = acc.replace(" ", "_")
        # Remove surrounding bars like "ref|NC_...|" if present
        if "|" in acc:
            # try to find an accession pattern inside
            m = ACCESSION_RE.search(acc)
            if m:
                acc = m.group(0)
        # After normalization, keep only the standard accession form if embedded
        m = ACCESSION_RE.search(acc)
        if m:
            acc = m.group(0)
        # Accept as-is if it already looks good
        if "." in acc:
            with_version.add(acc)
            no_version.add(acc.split(".")[0])
        else:
            # If version is missing, keep only no-version form
            no_version.add(acc)

    with open(list_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # Strategy 1: If the line contains a clear accession pattern, use it
            m = ACCESSION_RE.search(line.replace(" ", "_"))
            if m:
                _add(m.group(0))
                continue
            # Strategy 2: Join first two tokens like "NC" + "001806.1" -> NC_001806.1
            toks = line.split()
            if len(toks) >= 2 and toks[0].isalpha() and re.match(r"^\d+\.\d+$", toks[1]):
                _add(f"{toks[0]}_{toks[1]}")
                continue
            # Strategy 3: If first token already looks like an accession
            if toks:
                _add(toks[0])

    return with_version, no_version


def fasta_records(stream: Iterable[str]) -> Iterable[Tuple[str, str]]:
    """Yield (header, sequence) tuples from a FASTA stream.

    header includes the leading '>' and no trailing newline.
    sequence is concatenated without newlines.
    """
    header: Optional[str] = None
    seq_parts = []
    for line in stream:
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_parts)
            header = line.rstrip("\n")
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        yield header, "".join(seq_parts)


def extract_header_accession(header: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract accession with and without version from a FASTA header line.

    Returns (with_version, no_version). May return (None, None) if not found.
    """
    # Remove leading '>' if present
    h = header[1:] if header.startswith(">") else header
    # Look for standard accession anywhere in the header
    m = ACCESSION_RE.search(h)
    if m:
        acc = m.group(0)
        return acc, acc.split(".")[0]
    # Fallback: take the first whitespace token and try to clean it
    tok = h.split()[0] if h.split() else ""
    tok = tok.strip("|")
    tok = tok.replace(" ", "_")
    m2 = ACCESSION_RE.search(tok)
    if m2:
        acc = m2.group(0)
        return acc, acc.split(".")[0]
    return None, None


def filter_fasta(
    fasta_path: str,
    list_path: str,
    out_path: Optional[str] = None,
    ignore_version: bool = True,
) -> Tuple[int, int]:
    """Filter the FASTA at fasta_path, writing only records that match list_path.

    Returns (kept, total) counts.
    """
    targets_with_version, targets_no_version = parse_targets(list_path)

    out_fh = sys.stdout if not out_path or out_path == "-" else open(out_path, "w", encoding="utf-8")
    kept = 0
    total = 0
    try:
        with open(fasta_path, "r", encoding="utf-8") as fh:
            for header, seq in fasta_records(fh):
                total += 1
                acc_wv, acc_nv = extract_header_accession(header)
                match = False
                if acc_wv and acc_wv in targets_with_version:
                    match = True
                elif ignore_version and acc_nv and acc_nv in targets_no_version:
                    match = True

                if match:
                    kept += 1
                    # Write in wrapped format (60 chars per line) for readability
                    print(header, file=out_fh)
                    for i in range(0, len(seq), 60):
                        print(seq[i : i + 60], file=out_fh)
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()

    return kept, total


def main(argv: Optional[Iterable[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Filter a FASTA by a list of accessions")
    ap.add_argument("--fasta", required=True, help="Path to input FASTA (e.g., viruses.fasta)")
    ap.add_argument("--list", required=True, help="Path to list file (e.g., sam-viruses.txt)")
    ap.add_argument(
        "--out",
        default="-",
        help="Path to output FASTA (default: '-' for stdout)",
    )
    ap.add_argument(
        "--no-ignore-version",
        dest="ignore_version",
        action="store_false",
        help="Require exact version match (default matches by accession regardless of version)",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    kept, total = filter_fasta(
        fasta_path=args.fasta,
        list_path=args.list,
        out_path=args.out,
        ignore_version=args.ignore_version,
    )
    print(f"[filter_fasta_by_list] Kept {kept} / {total} records", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
