import argparse
import time
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import Entrez


REQUIRED_KEYS = {
    "ENTREZ_EMAIL",
    "ENTREZ_API_KEY",
    "START_YEAR",
    "END_YEAR",
    "SOFTWARE_LIST",
    "SLEEP_TIME",
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

    if not isinstance(env["ENTREZ_API_KEY"], str):
        raise ValueError("ENTREZ_API_KEY must be a string (can be empty)")

    if not isinstance(env["START_YEAR"], int) or not isinstance(env["END_YEAR"], int):
        raise ValueError("START_YEAR and END_YEAR must be integers")

    if env["START_YEAR"] > env["END_YEAR"]:
        raise ValueError("START_YEAR must be <= END_YEAR")

    if not isinstance(env["SOFTWARE_LIST"], list) or not env["SOFTWARE_LIST"]:
        raise ValueError("SOFTWARE_LIST must be a non-empty list")

    if not all(isinstance(x, str) and x.strip() for x in env["SOFTWARE_LIST"]):
        raise ValueError("SOFTWARE_LIST must contain non-empty strings")

    if not isinstance(env["SLEEP_TIME"], (int, float)) or env["SLEEP_TIME"] < 0:
        raise ValueError("SLEEP_TIME must be a non-negative number")

    env.setdefault("REQUEST_TIMEOUT", 30)
    if not isinstance(env["REQUEST_TIMEOUT"], (int, float)) or env["REQUEST_TIMEOUT"] <= 0:
        raise ValueError("REQUEST_TIMEOUT must be a positive number")

    return env


def build_pmc_query(software_name: str, year: int, open_access_only: bool = True) -> str:
    year_str = f'"{year}"[pdat]'
    oa_filter = ' AND "open access"[filter]' if open_access_only else ""
    return f'"{software_name}" AND {year_str}{oa_filter}'


def get_article_count(query: str, request_timeout: float) -> int:
    delays = [0.0, 0.5, 1.0]

    for attempt, delay in enumerate(delays):
        if delay > 0:
            time.sleep(delay)
        try:
            handle = Entrez.esearch(
                db="pmc",
                term=query,
                retmax=0,
                timeout=request_timeout,
            )
            record = Entrez.read(handle)
            handle.close()
            return int(record.get("Count", 0))
        except Exception as e:
            print(f"  [Retry {attempt+1}/3] esearch failed: {e}")

    print("  !! esearch failed after 3 attempts")
    return 0


def collect_software_usage_statistics(
    software_list: list[str],
    start_year: int,
    end_year: int,
    sleep_time: float,
    request_timeout: float,
) -> pd.DataFrame:
    print("\n" + "=" * 70)
    print("PMC Software Usage Statistics")
    print("=" * 70)
    print(f"Tools: {len(software_list)} | Years: {start_year}-{end_year}")
    total_queries = len(software_list) * (end_year - start_year + 1)
    print(f"Total queries: {total_queries}")
    print("=" * 70)
    print("\nTip: data is saved after each tool, so you can Ctrl+C and rerun.")
    print("=" * 70)

    all_results: list[dict[str, Any]] = []
    current_query = 0

    detail_file = "pmc_software_usage_by_year.tsv"
    summary_file = "pmc_software_usage_summary.tsv"

    for idx, software in enumerate(software_list, start=1):
        print(f"\n[{idx}/{len(software_list)}] {software}")
        print("-" * 70)

        for year in range(start_year, end_year + 1):
            current_query += 1
            progress = (current_query / total_queries) * 100
            print(f"  {year} ", end="", flush=True)

            query = build_pmc_query(software, year, open_access_only=True)

            try:
                count = get_article_count(query, request_timeout=request_timeout)
                print(f"-> {count:4d} articles  [{progress:5.1f}%]")
            except Exception as e:
                print(f"-> ERROR: {e}")
                count = 0

            all_results.append({"software": software, "year": year, "n_articles": count})
            time.sleep(sleep_time)

        pd.DataFrame(all_results).to_csv(detail_file, sep="\t", index=False)
        print(f"  ✓ Saved: {detail_file}")

    df_results = pd.DataFrame(all_results)

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE - Summary")
    print("=" * 70)

    summary_total = df_results.groupby("software")["n_articles"].sum().sort_values(ascending=False)

    recent_df = df_results[df_results["year"].between(2020, 2024)]
    summary_recent = recent_df.groupby("software")["n_articles"].sum().sort_values(ascending=False)

    latest_df = df_results[df_results["year"].between(2022, 2024)]
    summary_latest = latest_df.groupby("software")["n_articles"].sum().sort_values(ascending=False)

    summary_all = pd.DataFrame(
        {
            "software": summary_total.index,
            "total_2000_2025": summary_total.values,
            "recent_2020_2024": [summary_recent.get(s, 0) for s in summary_total.index],
            "latest_2022_2024": [summary_latest.get(s, 0) for s in summary_total.index],
        }
    )
    summary_all.to_csv(summary_file, sep="\t", index=False)

    print("\nSaved files:")
    print(f"  1) {detail_file}")
    print(f"  2) {summary_file}")

    return df_results


def main() -> int:
    ap = argparse.ArgumentParser(description="PMC software mention counter (config-driven).")
    ap.add_argument("config", help="Path to config file, e.g. config_tools.txt")
    args = ap.parse_args()

    cfg = load_py_config(args.config)

    Entrez.email = cfg["ENTREZ_EMAIL"]
    if cfg["ENTREZ_API_KEY"]:
        Entrez.api_key = cfg["ENTREZ_API_KEY"]

    print("\nPMC Software Usage Analysis")
    print("=" * 70)
    print("This script will:")
    print("  - Search PMC for software mentions (open access filter on)")
    print("  - Count unique articles (Count from esearch)")
    print(f"  - Years: {cfg['START_YEAR']}-{cfg['END_YEAR']}")
    print("=" * 70)

    collect_software_usage_statistics(
        software_list=cfg["SOFTWARE_LIST"],
        start_year=cfg["START_YEAR"],
        end_year=cfg["END_YEAR"],
        sleep_time=float(cfg["SLEEP_TIME"]),
        request_timeout=float(cfg["REQUEST_TIMEOUT"]),
    )
    print("\nAnalysis complete!\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
