
### Set your working directory:
#setwd("XXXXXX")

library(tidyverse)
library(scales)

# input files
files <- tribble(
  ~topic_file, ~topic_name,
  "EST_or_microarray_records.tsv", "EST & microarray",
  "RNA-seq_records.filtered.tsv", "RNA-seq",
  "scRNA-seq_records.tsv", "scRNA-seq",
  "spatial_or_multimodal_integration_records.tsv", "spatial+multimodal",
  "AI_records.tsv", "AI"
)

# categories
topic_to_category <- c(
  "EST & microarray"      = "EST & microarray",
  "RNA-seq"               = "RNA-seq",
  "scRNA-seq"             = "scRNA-seq",
  "spatial+multimodal"    = "spatial+multimodal",
  "AI"                    = "AI"
)

# Journals to KEEP (remove Nat Methods + PLOS Comp Biol)
keep_journals <- paste(
  c("Briefings in Bioinformatics",
    "^Bioinformatics$",
    "Nature Biotechnology",
    "Genome Research",
    "Genome Biology",
    "Nucleic Acids Research"),
  collapse = "|"
)

# load all records
read_one <- function(path, topic_name){
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(topic = topic_name) %>%
    mutate(year = suppressWarnings(as.integer(year))) %>%
    filter(!is.na(year), year >= 2000, year <= 2025)
}

df <- files %>%
  mutate(data = map2(topic_file, topic_name, read_one)) %>%
  select(-topic_file) %>%
  unnest(data)

# filter journals
df_filt <- df %>%
  filter(str_detect(journal, regex(keep_journals, ignore_case = TRUE)))

# annual counts by 4 big categories + AI
df_counts <- df_filt %>%
  mutate(category = topic_to_category[topic]) %>%
  count(year, category, name = "n") %>%
  complete(year = 2000:2025, category, fill = list(n = 0))

tech4 <- df_counts %>%
  filter(category != "AI")

ai <- df_counts %>%
  filter(category == "AI") %>%
  select(year, n_ai = n)

tech_total <- tech4 %>%
  group_by(year) %>%
  summarise(n_tech = sum(n), .groups = "drop")


### plot stacked area + AI% dashed line (right axis)
# AI proportion among tech papers
df_prop <- tech_total %>%
  left_join(ai, by = "year") %>%
  mutate(n_ai = replace_na(n_ai, 0),
         ai_share = ifelse(n_tech == 0, NA_real_, n_ai / n_tech))

max_left <- max(tech_total$n_tech, na.rm = TRUE)
scale_factor_pct <- max_left*6   # ai_share in [0,20%], so multiply by max_left to overlay

cat_cols <- c(
  "EST & microarray"   = "#039bbc",  # Teal
  "RNA-seq"            = "#ce1575",  # Pink
  "scRNA-seq"          = "#6c3e98",  # Purple
  "spatial+multimodal" = "#e68025"   # Orange
)

p1 <- ggplot() +
  geom_area(
    data = tech4,
    aes(x = year, y = n, fill = category),
    alpha = 0.85, colour = NA
  ) +
  scale_fill_manual(values = cat_cols) +
  # AI% dashed line overlay (scaled onto left axis)
  geom_line(
    data = df_prop,
    aes(x = year, y = ai_share * scale_factor_pct),
    linewidth = 1.1,
    linetype = "dashed",
    colour = "#a2be3b"
  ) +
  scale_y_continuous(
    labels = comma,
    name = "Articles per year",
    sec.axis = sec_axis(
      ~ . / scale_factor_pct,
      labels = percent_format(accuracy = 1),
      name = "AI proportion"
    )
  ) +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  labs(x = NULL, fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title.position = "plot"
  )

ggsave("Fig2a_stacked_area_chart.pdf", p1, width = 11, height = 6)

