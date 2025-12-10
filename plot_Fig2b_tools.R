
### Set your working directory:
#setwd("XXXXXX")

library(tidyverse)
library(scales)

file <- "pmc_software_usage_by_year.tsv"

df <- readr::read_tsv(file, show_col_types = FALSE) %>%
  mutate(
    year = as.integer(year),
    software = as.character(software),
    n_articles = as.integer(n_articles)
  )


remove_tools <- c("STAR", "TRIAGE", "Salmon", "GLUE", "MAGIC", "Harmony", "FLAIR")
df_clean <- df %>% filter(!software %in% remove_tools)


release_map <- tribble(
  ~software,        ~category,               ~release_year,
  "CAP3",           "EST & microarray",      1999,
  "Limma",          "EST & microarray",      2004,
  "TopHat",         "RNA-seq",               2009, 
  "edgeR",          "RNA-seq",               2010, 
  "Cufflinks",      "RNA-seq",               2010,
  "DESeq",          "RNA-seq",               2010,
  "RSEM",           "RNA-seq",               2011,
  "Monocle",        "scRNA-seq",             2014,
  "featureCounts",  "RNA-seq",               2014,
  "HTSeq",          "RNA-seq",               2015,
  "StringTie",      "RNA-seq",               2015,
  "Ballgown",       "RNA-seq",               2015,
  "HISAT",          "RNA-seq",               2015,
  "Seurat",         "scRNA-seq",             2015,
  "kallisto",       "RNA-seq",               2016,
  "minimap",        "RNA-seq",               2016,
  

  "Cell Ranger",    "scRNA-seq",             2017,
  "SCENIC",         "scRNA-seq",             2017,
  "Scanpy",         "scRNA-seq",             2018,
  
  "MOFA",           "Spatial & multimodal",  2018,
  
  
  "NicheNet",       "Spatial & multimodal",  2019,
  "SCTransform",    "scRNA-seq",             2019,
  "CellPhoneDB",    "Spatial & multimodal",  2020,
  "totalVI",        "Spatial & multimodal",  2021,
  "CellChat",       "Spatial & multimodal",  2021,
  
  "CellOracle",     "AI",                    2023,
  "Geneformer",     "AI",                    2023,
  "CellBender",     "AI",                    2023,
  "scGPT",          "AI",                    2024,
  "scFoundation",   "AI",                    2024
)

#software_order <- rev(release_map$software)
software_order <- release_map$software

# 把 release_map join 到数据，并过滤只画你这份清单里的软件
df_plot <- df_clean %>%
  inner_join(release_map, by = "software") %>%
  mutate(
    software = factor(software, levels = software_order),
    # release_year 之前留空（NA），其余用 log1p(count)
    fill_value = if_else(year < release_year, NA_real_, log1p(n_articles))
  )

# 标注点：每个软件在 release_year 那一格打点
df_mark <- release_map %>%
  mutate(software = factor(software, levels = software_order))

# =========================
# Heatmap：按你顺序 + 空白 + release_year 标注
# =========================
p_heat_release <- ggplot(df_plot, aes(x = year, y = software, fill = fill_value)) +
  geom_tile() +
  geom_point(
    data = df_mark,
    aes(x = release_year, y = software),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.9,
    stroke = 0.25,
    fill = "white"
  ) +
  scale_x_continuous(breaks = seq(min(df_plot$year), max(df_plot$year), by = 5)) +
  scale_fill_viridis_c(name = "log1p(count)", na.value = "transparent") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey85")
  ) +
  labs(
    x = "Year",
    y = ""
  )
p_heat_release
ggsave("Figure2b_heatmap.pdf", p_heat_release, width = 10, height = 8, dpi = 300)



#################################################################################
### List top 3 tools per year (2021–2025)
top3_2021_2025 <- df_clean %>%
  filter(year >= 2021, year <= 2025) %>%
  group_by(year) %>%
  arrange(desc(n_articles), software) %>%
  slice_head(n = 3) %>%
  ungroup()

top3_2021_2025 %>%
  arrange(year, desc(n_articles), software) %>%
  select(year, software, n_articles) %>%
  print(n = Inf)

top3_2021_2025 %>%
  arrange(year, desc(n_articles), software) %>%
  group_by(year) %>%
  summarise(
    top3 = paste0(software, " (", n_articles, ")", collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(line = paste0(year, ": ", top3)) %>%
  pull(line) %>%
  cat(sep = "\n")
cat("\n")

# 2021: Limma (6517) | edgeR (4780) | HTSeq (2710)
# 2022: Limma (8760) | edgeR (5210) | Seurat (3918)
# 2023: Limma (7850) | Seurat (5432) | edgeR (4941)
# 2024: Limma (8358) | Seurat (7492) | edgeR (5191)
# 2025: Seurat (9357) | Limma (9075) | edgeR (5214)



#################################################################################
### plot top 3 tools in each of 2021-25 for the line plot (inset)
tools <- c("Limma", "edgeR", "Seurat", "HTSeq")
df_line <- df_clean %>%
  filter(year >= 2021, year <= 2025, software %in% tools) %>%
  mutate(
    software = factor(software, levels = tools),
    year = as.integer(year),
    n_articles = as.integer(n_articles)
  )

p_line <- ggplot(df_line, aes(x = year, y = n_articles, colour = software, group = software)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2021:2025) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    x = "",
    y = "Number of PMC articles",
    colour = "Tool"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )


ggsave("Figure2b_top_tools_2021_2025_line.pdf", p_line, width = 7.5, height = 4.8)
#ggsave("Figure2b_top_tools_2021_2025_line.png", p_line, width = 7.5, height = 4.8, dpi = 300)
