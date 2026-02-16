#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import List, Optional

import pandas as pd
from Bio import AlignIO
from Bio.Align import AlignInfo


def _infer_ogs_id(filename: str) -> str:
    name = filename
    name = re.sub(r"\.(fasta|fas|fa)$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"_cat.*$", "", name)
    return name


def _gc_percent_from_alignment(alignment) -> float:
    gc = 0
    atgc = 0
    for rec in alignment:
        seq = str(rec.seq).upper()
        for ch in seq:
            if ch in {"A", "C", "G", "T"}:
                atgc += 1
                if ch in {"C", "G"}:
                    gc += 1
    if atgc == 0:
        return float("nan")
    return (gc / atgc) * 100.0


def compute_alignment_stats(aln_path: Path) -> dict:
    alignment = AlignIO.read(str(aln_path), "fasta")

    align_len = int(alignment.get_alignment_length())
    no_seq = int(len(alignment))

    if no_seq == 0 or align_len == 0:
        return {
            "OGS": _infer_ogs_id(aln_path.name),
            "sub_rate": float("nan"),
            "Avg_record": float("nan"),
            "align_len": align_len,
            "No_seq": no_seq,
            "meanGC": float("nan"),
        }

    total_gaps = 0
    for rec in alignment:
        total_gaps += str(rec.seq).count("-")
    avg_record = total_gaps / no_seq

    summary = AlignInfo.SummaryInfo(alignment)
    consensus = summary.dumb_consensus(ambiguous="X")
    consensus_str = str(consensus).upper()
    sub_rate = consensus_str.count("X") / len(consensus_str) if consensus_str else float("nan")

    mean_gc = _gc_percent_from_alignment(alignment)

    return {
        "OGS": _infer_ogs_id(aln_path.name),
        "sub_rate": float(sub_rate),
        "Avg_record": float(avg_record),
        "align_len": align_len,
        "No_seq": no_seq,
        "meanGC": float(mean_gc),
    }


def _iter_alignment_files(input_dir: Path, list_file: Optional[Path]) -> List[Path]:
    if list_file is not None:
        files: List[Path] = []
        for line in list_file.read_text().splitlines():
            name = line.strip()
            if not name:
                continue
            files.append(input_dir / name)
        return files

    patterns = ["*.fa", "*.fas", "*.fasta"]
    out: List[Path] = []
    for pat in patterns:
        out.extend(sorted(input_dir.glob(pat)))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute per-alignment statistics used for ML: sub_rate, Avg_record, align_len, No_seq, meanGC."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing alignment FASTA files.",
    )
    parser.add_argument(
        "--list-file",
        type=Path,
        default=None,
        help="Optional text file listing alignment filenames (one per line) inside --input-dir.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output TSV path.",
    )

    args = parser.parse_args()

    input_dir: Path = args.input_dir
    if not input_dir.exists():
        raise FileNotFoundError(f"Input dir not found: {input_dir}")

    files = _iter_alignment_files(input_dir, args.list_file)
    rows = []
    for p in files:
        if not p.exists():
            continue
        if p.is_file() and p.stat().st_size == 0:
            rows.append(
                {
                    "OGS": _infer_ogs_id(p.name),
                    "sub_rate": float("nan"),
                    "Avg_record": float("nan"),
                    "align_len": 0,
                    "No_seq": 0,
                    "meanGC": float("nan"),
                }
            )
            continue

        try:
            rows.append(compute_alignment_stats(p))
        except Exception:
            rows.append(
                {
                    "OGS": _infer_ogs_id(p.name),
                    "sub_rate": float("nan"),
                    "Avg_record": float("nan"),
                    "align_len": float("nan"),
                    "No_seq": float("nan"),
                    "meanGC": float("nan"),
                }
            )

    df = pd.DataFrame(rows, columns=["OGS", "sub_rate", "Avg_record", "align_len", "No_seq", "meanGC"])
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
