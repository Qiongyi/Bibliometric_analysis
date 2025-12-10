#!/usr/bin/env python3
import re
import pandas as pd
from pathlib import Path

INFILE = "RNA-seq_records.tsv"
OUT_KEEP = "RNA-seq_records.filtered.tsv"
OUT_DROP = "RNA-seq_records.removed_sc_only.tsv"

# ---------- Patterns ----------
# RNA-seq (bulk) keywords
RNA_SEQ_PAT = re.compile(
    r"\b("
    r"rna[\s\-]?seq|"                 # RNA-seq, RNA seq, RNAseq, RNA-Seq
    r"rna[\s\-]?sequencing"           # RNA sequencing / RNA-sequencing
    r")\b",
    re.IGNORECASE
)

# scRNA-seq / single-cell RNA-seq variants to strip first
SCRNA_SEQ_PAT = re.compile(
    r"\b("
    r"scrna[\s\-]?seq|"               # scRNA-seq, scRNAseq, scRNA seq, scRNA-Seq
    r"single[\s\-]?cell[\s\-]?rna[\s\-]?seq|"          # single-cell RNA-seq etc
    r"single[\s\-]?cell[\s\-]?rna[\s\-]?sequencing|"   # single cell RNA sequencing
    r"single[\s\-]?cell[\s\-]?rna[\s\-]?sequencing"    # (dup ok)
    r")\b",
    re.IGNORECASE
)

def normalise_text(s: str) -> str:
    if pd.isna(s):
        return ""
    return str(s)

def is_sc_only(title: str, abstract: str) -> bool:
    text = f"{normalise_text(title)} {normalise_text(abstract)}"
    if not SCRNA_SEQ_PAT.search(text):
        return False

    # Remove scRNA-seq terms, then re-check for bulk RNA-seq terms
    stripped = SCRNA_SEQ_PAT.sub(" ", text)
    stripped = re.sub(r"\s+", " ", stripped).strip()

    # If bulk RNA-seq terms are absent after stripping, it's sc-only
    return RNA_SEQ_PAT.search(stripped) is None

def main():
    infile = Path(INFILE)
    if not infile.exists():
        raise FileNotFoundError(f"Cannot find input file: {INFILE}")

    df = pd.read_csv(INFILE, sep="\t", dtype=str, keep_default_na=False)

    for col in ["title", "abstract"]:
        if col not in df.columns:
            raise ValueError(f"Missing column '{col}' in {INFILE}. Found columns: {list(df.columns)}")

    mask_sc_only = df.apply(lambda r: is_sc_only(r.get("title", ""), r.get("abstract", "")), axis=1)

    df_drop = df[mask_sc_only].copy()
    df_keep = df[~mask_sc_only].copy()

    df_keep.to_csv(OUT_KEEP, sep="\t", index=False)
    df_drop.to_csv(OUT_DROP, sep="\t", index=False)

    print("=== RNA-seq sc-only filter ===")
    print(f"Input:   {len(df)}")
    print(f"Kept:    {len(df_keep)}  -> {OUT_KEEP}")
    print(f"Removed: {len(df_drop)}  -> {OUT_DROP}")

if __name__ == "__main__":
    main()
