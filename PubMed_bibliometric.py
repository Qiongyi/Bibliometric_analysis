from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import Entrez, Medline


REQUIRED_KEYS = {
    "ENTREZ_EMAIL",
    "JOURNALS",
    "TOPICS",
    "SLEEP_TIME",
    "START_YEAR",
    "END_YEAR",
}


def load_py_config(path: str) -> dict[str, Any]:
    cfg_path = Path(path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")

    code = cfg_path.read_text(encoding="utf-8")

    safe_builtins = {
        "True": True,
        "False": False,
        "None": None,
        "dict": dict,
        "list": list,
        "set": set,
        "tuple": tuple,
        "str": str,
        "int": int,
        "float": float,
        "len": len,
        "range": range,
    }

    env: dict[str, Any] = {}
    exec(compile(code, str(cfg_path), "exec"), {"__builtins__": safe_builtins}, env)

    missing = [k for k in REQUIRED_KEYS if k not in env]
    if missing:
        raise ValueError(f"Missing required config keys: {missing}")

    if not isinstance(env["ENTREZ_EMAIL"], str) or not env["ENTREZ_EMAIL"].strip():
        raise ValueError("ENTREZ_EMAIL must be a non-empty string")

    if not isinstance(env["JOURNALS"], dict) or not env["JOURNALS"]:
        raise ValueError("JOURNALS must be a non-empty dict")

    if not isinstance(env["TOPICS"], dict) or not env["TOPICS"]:
        raise ValueError("TOPICS must be a non-empty dict")

    if not isinstance(env["SLEEP_TIME"], (int, float)) or env["SLEEP_TIME"] < 0:
        raise ValueError("SLEEP_TIME must be a non-negative number")

    if not isinstance(env["START_YEAR"], int) or not isinstance(env["END_YEAR"], int):
        raise ValueError("START_YEAR and END_YEAR must be integers")

    if env["START_YEAR"] > env["END_YEAR"]:
        raise ValueError("START_YEAR must be <= END_YEAR")

    env.setdefault("ENTREZ_API_KEY", "")
    if not isinstance(env["ENTREZ_API_KEY"], str):
        raise ValueError("ENTREZ_API_KEY must be a string")

    return env


def build_query(topic_query: str, journal_name: str, year: int) -> str:
    yy = f'"{year}"[Date - Publication]'
    jj = f'"{journal_name}"[Journal]'
    return f"({topic_query}) AND ({jj}) AND ({yy})"


def get_pmids_for_query(query: str, retmax: int = 100000) -> list[str]:
    delays = [0.0, 0.5, 1.0]
    for attempt, delay in enumerate(delays):
        if delay > 0:
            time.sleep(delay)
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
            record = Entrez.read(handle)
            handle.close()
            return record.get("IdList", [])
        except Exception as e:
            print(f"  [Retry {attempt + 1}/3] esearch failed: {e}")

    print("  !! esearch failed after 3 attempts, returning empty list")
    return []


def fetch_records(pmids: list[str], sleep_time: float) -> list[dict[str, str]]:
    records_data: list[dict[str, str]] = []
    batch_size = 200

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i : i + batch_size]
        ids_str = ",".join(batch)

        handle = Entrez.efetch(db="pubmed", id=ids_str, rettype="medline", retmode="text")
        medline_records = Medline.parse(handle)

        for rec in medline_records:
            pmid = rec.get("PMID", "")
            title = rec.get("TI", "")
            abstract = rec.get("AB", "")
            journal = rec.get("JT", "")
            year = rec.get("DP", "")[:4]

            records_data.append(
                {
                    "pmid": pmid,
                    "year": year,
                    "journal": journal,
                    "title": title,
                    "abstract": abstract,
                }
            )

        handle.close()
        time.sleep(sleep_time)

    return records_data


def collect_and_save_by_topic(
    topics: dict[str, str],
    journals: dict[str, str],
    start_year: int,
    end_year: int,
    sleep_time: float,
) -> None:
    for topic_name, topic_query in topics.items():
        print(f"\n=== Topic: {topic_name} ===")
        all_pmids: set[str] = set()

        for journal_label, journal_field in journals.items():
            for year in range(start_year, end_year + 1):
                query = build_query(topic_query, journal_field, year)
                try:
                    pmids = get_pmids_for_query(query)
                except Exception as e:
                    print(f"Error in esearch for {topic_name} / {journal_label} / {year}: {e}")
                    pmids = []
                all_pmids.update(pmids)
                time.sleep(sleep_time)

        print(f"  Total PMIDs for {topic_name}: {len(all_pmids)}")
        if not all_pmids:
            continue

        records = fetch_records(list(all_pmids), sleep_time=sleep_time)
        df = pd.DataFrame(records)

        safe_name = topic_name.replace(" ", "_").replace("/", "_")
        out_file = f"{safe_name}_records.tsv"
        df.to_csv(out_file, sep="\t", index=False)
        print(f"  Saved {len(df)} records to {out_file}")

        yearly_counts = (
            df.groupby(["journal", "year"])
            .size()
            .reset_index(name="n_articles")
            .sort_values(["journal", "year"])
        )
        yearly_counts_file = f"{safe_name}_yearly_counts.tsv"
        yearly_counts.to_csv(yearly_counts_file, sep="\t", index=False)
        print(f"  Saved yearly counts to {yearly_counts_file}")


def get_total_count_for_journal_year(journal_field: str, year: int) -> int:
    year_str = f'"{year}"[Date - Publication]'
    journal_clause = f'"{journal_field}"[Journal]'
    query = f"({journal_clause}) AND ({year_str})"

    delays = [0.0, 0.5, 1.0]
    for attempt, delay in enumerate(delays):
        if delay > 0:
            time.sleep(delay)
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
            record = Entrez.read(handle)
            handle.close()
            return int(record.get("Count", 0))
        except Exception as e:
            print(f"  [Retry {attempt + 1}/3] esearch (count) failed for {journal_field} {year}: {e}")

    print(f"  !! esearch (count) failed after 3 attempts for {journal_field} {year}, returning 0")
    return 0


def collect_journal_yearly_totals(
    journals: dict[str, str],
    start_year: int,
    end_year: int,
    sleep_time: float,
) -> None:
    rows: list[dict[str, Any]] = []

    for journal_label, journal_field in journals.items():
        print(f"\n=== Journal totals: {journal_label} ===")
        for year in range(start_year, end_year + 1):
            try:
                count = get_total_count_for_journal_year(journal_field, year)
                print(f"  {year}: {count} articles")
            except Exception as e:
                print(f"Error in esearch for totals {journal_label} / {year}: {e}")
                count = 0

            rows.append(
                {
                    "journal_label": journal_label,
                    "journal_query_name": journal_field,
                    "year": year,
                    "n_articles": count,
                }
            )
            time.sleep(sleep_time)

    df_totals = pd.DataFrame(rows)
    df_totals.to_csv("journal_yearly_totals.tsv", sep="\t", index=False)
    print("\nSaved journal yearly totals to journal_yearly_totals.tsv")


def main() -> int:
    ap = argparse.ArgumentParser(description="PubMed bibliometric analysis")
    ap.add_argument("config", help="Path to config file, e.g. config_tiab.txt")
    args = ap.parse_args()

    cfg = load_py_config(args.config)

    Entrez.email = cfg["ENTREZ_EMAIL"]
    if cfg.get("ENTREZ_API_KEY"):
        Entrez.api_key = cfg["ENTREZ_API_KEY"]

    collect_and_save_by_topic(
        topics=cfg["TOPICS"],
        journals=cfg["JOURNALS"],
        start_year=cfg["START_YEAR"],
        end_year=cfg["END_YEAR"],
        sleep_time=float(cfg["SLEEP_TIME"]),
    )

    collect_journal_yearly_totals(
        journals=cfg["JOURNALS"],
        start_year=cfg["START_YEAR"],
        end_year=cfg["END_YEAR"],
        sleep_time=float(cfg["SLEEP_TIME"]),
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
