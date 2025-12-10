#!/usr/bin/env python3
import re
import pandas as pd
from pathlib import Path

INPUT_FILES = [
    "EST_or_microarray_records.tsv",
    "RNA-seq_records.filtered.tsv",
    "scRNA-seq_records.tsv",
    "spatial_or_multimodal_integration_records.tsv",
]

OUTFILE = "AI_records.tsv"

# AI keyword patterns
# 'deep learning' OR 'deep neural network*' OR
# 'graph neural network*' OR 'graph convolutional network*' OR
# 'deep generative model*' OR 'generative adversarial network*' OR
# 'variational autoencoder*' OR 'transformer model*' OR 'transformer-based'
AI_PATTERNS = [
    r"\bdeep\s+learning\b",
    r"\bdeep\s+neural\s+networks?\b",
    r"\bgraph\s+neural\s+networks?\b",
    r"\bgraph\s+convolutional\s+networks?\b",
    r"\bdeep\s+generative\s+models?\b",
    r"\bgenerative\s+adversarial\s+networks?\b",
    r"\bvariational\s+autoencoders?\b",
    r"\btransformer\s+models?\b",
    r"\btransformer[-\s]?based\b",
]

AI_REGEX = re.compile(r"(?:%s)" % "|".join(AI_PATTERNS), flags=re.IGNORECASE)

def file_exists_or_die(p: str):
    if not Path(p).exists():
        raise FileNotFoundError(f"Cannot find file: {p}")

def load_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    # Ensure columns exist
    for col in ("title", "abstract"):
        if col not in df.columns:
            raise ValueError(f"Missing column '{col}' in {path}. Columns: {list(df.columns)}")
    return df

def ai_hit(df: pd.DataFrame) -> pd.Series:
    text = (df["title"].fillna("") + " " + df["abstract"].fillna("")).astype(str)
    return text.str.contains(AI_REGEX, na=False)

def main():
    all_hits = []
    per_file_stats = []

    for f in INPUT_FILES:
        file_exists_or_die(f)
        df = load_tsv(f)
        mask = ai_hit(df)
        hits = df.loc[mask].copy()
        hits["source_file"] = f  # keep provenance (optional but useful)
        all_hits.append(hits)
        per_file_stats.append((f, len(df), int(mask.sum())))

    if not all_hits:
        print("No inputs processed.")
        return

    out = pd.concat(all_hits, ignore_index=True)

    if "pmid" in out.columns:
        out = out.drop_duplicates(subset=["pmid"], keep="first")
    else:
        out = out.drop_duplicates(subset=[c for c in ["journal", "year", "title"] if c in out.columns], keep="first")

    # Optional: sort by year then journal
    if "year" in out.columns:
        out["year_num"] = pd.to_numeric(out["year"], errors="coerce")
        out = out.sort_values(["year_num", "journal"], na_position="last").drop(columns=["year_num"])

    out.to_csv(OUTFILE, sep="\t", index=False)

    print("=== AI record extraction (title+abstract) ===")
    for f, n_total, n_hit in per_file_stats:
        print(f"{f:40s}  total={n_total:6d}  AI_hits={n_hit:5d}")
    print(f"\nSaved total unique AI hits: {len(out)} -> {OUTFILE}")

if __name__ == "__main__":
    main()
