# Bibliometric Analysis Workflow (2000–2025)

This repository contains the scripts and intermediate outputs used to reproduce our bibliometric analyses and figures. We performed two main analyses:

1. **PubMed journal-level trends** (Figure 2a)
2. **Tool adoption in PMC full text** (Figure 2b)

---

## Overview

### PubMed Journal-Level Trends (Figure 2a)
We queried PubMed for articles published between 2000 and 2025 in six leading bioinformatics journals:
- Briefings in Bioinformatics
- Bioinformatics
- Nature Biotechnology
- Genome Biology
- Genome Research
- Nucleic Acids Research

### Tool Adoption in PMC Full Text (Figure 2b)
We quantified tool adoption using **7,078,831** PubMed Central (PMC) open-access full-text articles published between 2000 and 2025.

---

## Analysis 1: PubMed Trends in Six Journals (Figure 2a)

### Step 1.1: Run PubMed Queries (by Topic/Category)
```bash
python PubMed_bibliometric.py config_tiab.txt
```

This script searches PubMed using the topics and journals defined in `config_tiab.txt`.

**Before running:**
- Edit `config_tiab.txt` to match your use case:
  - `ENTREZ_EMAIL`: Your email (required by NCBI Entrez)
  - `ENTREZ_API_KEY`: Optional but recommended (faster and more stable)
  - `START_YEAR`, `END_YEAR`: Time range for analysis
  - `JOURNALS`: Journal display names and PubMed journal field names
  - `TOPICS`: PubMed query strings (Title/Abstract-focused in our setup)
  - `SLEEP_TIME`: Delay between requests to avoid rate limits

**Outputs:**
- For each topic/category:
  - `*_records.tsv`: PMID-level records (year/journal/title/abstract, etc.)
  - `*_yearly_counts.tsv`: Yearly counts grouped by journal
- `journal_yearly_totals.tsv`: Total number of PubMed articles per journal per year
  - Useful for normalizing topic counts by overall journal output if needed

### Step 1.2: Remove "scRNA-seq-only" Papers from RNA-seq Category
```bash
python filter_rnaseq_remove_scrnaseq_only.py
```

This step cleans the RNA-seq category by removing papers that are only about scRNA-seq (owing to overlapping terminology).

**Output:**
- `RNA-seq_records.filtered.tsv`

### Step 1.3: Identify AI-Related Papers (for AI Proportion)
```bash
python search_AI_records.py
```

This step searches for AI-related papers and generates records/counts that are later used to compute the AI proportion shown as the dashed green line in Figure 2a.

### Step 1.4: Plot Figure 2a
```bash
Rscript plot_Fig2a_tech.R
```

This script reads the "records.tsv" tables and generates the stacked area plot with the AI proportion overlay.

---

## Analysis 2: Tool Adoption Using PMC Open-Access Full Text (Figure 2b)

### Step 2.1: Run PMC Full-Text Tool Counting
```bash
python PMC_tools.py config_fulltext_tools.txt
```

This script queries PMC and counts (per year) how many open-access full-text articles mention each tool listed in the config.

**Before running:**
- Edit `config_fulltext_tools.txt` to set:
  - `ENTREZ_EMAIL`, `ENTREZ_API_KEY`
  - `START_YEAR`, `END_YEAR`
  - `SOFTWARE_LIST`: List of tools to search for
  - `SLEEP_TIME` (optional)
  - `REQUEST_TIMEOUT`

**Output:**
- `pmc_software_usage_by_year.tsv`: A year-by-year table with the number of PMC open-access articles that mention each tool

### Step 2.2: Plot Figure 2b
```bash
Rscript plot_Fig2b_tools.R
```

This script generates the tool adoption panel, including:
- Annual counts per tool shown as heatmap
- Publication-year markers (white circles in heatmap)
- Highlighting of the top tools in recent five years
- Inset trend plot for the top tools

---

## Requirements

### Python Dependencies
```bash
pip install biopython pandas
```

### R Dependencies
```R
install.packages(c("tidyverse", "scales"))
```

---

## Citation

If you use this workflow in your research, please cite our paper:
```
*Coming soon.*
```

---

## License

[MIT]

---

## Contact

For questions or issues, please open an issue on GitHub or contact [q.zhao@uq.edu.au].