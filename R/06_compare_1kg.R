# 06_compare_1kg.R — 같은 PCA, 다른 해상도: 1000 Genomes와 PBMC의 PC1–PC2를 나란히
# 입력: work/pbmc_annot.rds (04_annotate.R 결과)
#       ../human-diversity-1kg/work/B_pca.eigenvec, ../human-diversity-1kg/work/pops.txt
#       (1000 Genomes 프로젝트 폴더가 이 폴더 옆에 있어야 한다. 그 프로젝트의 07_pca.R 결과)
# 출력: 그림 11 (results/figures/fig11_same_pca_two_scales.png)
#
# 왜 하나:
#   두 프로젝트를 잇는 메시지 "같은 PCA를 사람 집단에 쓰면 조상이, 세포 집단에 쓰면 세포 종류가 보인다"를
#   그림 한 장으로 보여 준다. README 맨 위에 넣는다.
#   두 그림 모두 PCA(PC1–PC2)로 그린다. UMAP을 쓰면 "같은 방법"이 아니게 되기 때문이다.
#
# 실행 방법: Session → Restart R 후, 처음부터 끝까지 한 번에 실행 (판단할 값 없음)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

kg_dir <- "../human-diversity-1kg"
stopifnot("1000 Genomes 폴더를 찾을 수 없습니다" = dir.exists(kg_dir))

# ── 1. 1000 Genomes PCA 점수 (사람 2,491명) ──────────────────
pops <- read_tsv(file.path(kg_dir, "work/pops.txt"), col_types = cols(.default = "c")) |>
  rename(IID = `#IID`)
kg <- read_tsv(file.path(kg_dir, "work/B_pca.eigenvec"), show_col_types = FALSE) |>
  rename(IID = `#IID`) |>
  mutate(IID = as.character(IID)) |>
  inner_join(pops, by = "IID")
stopifnot(nrow(kg) == 2491)

super_labels  <- c(AFR = "아프리카", EUR = "유럽", SAS = "남아시아", EAS = "동아시아", AMR = "아메리카")
super_colours <- c(AFR = "#2a78d6", EUR = "#eb6834", SAS = "#1baf7a", EAS = "#eda100", AMR = "#e87ba4")
kg <- kg |> mutate(group = factor(SUPER, levels = names(super_labels)))

# ── 2. PBMC PCA 점수 (세포 2,511개) ─────────────────────────
pbmc <- readRDS("work/pbmc_annot.rds")
pb <- as_tibble(Embeddings(pbmc, "pca")[, 1:2], rownames = "cell") |>
  rename(PC1 = PC_1, PC2 = PC_2) |>
  mutate(group = pbmc$cell_type)

# ── 3. PC 부호 맞추기 ──────────────────────────────────────
# PC의 부호(+/−)는 임의라서, 뒤집어도 의미는 같다.
# 보기 쉽게 "가장 먼저 갈라지는 쪽"(아프리카, 단핵구)이 두 그림 모두 왼쪽에 오도록 맞춘다.
if (mean(kg$PC1[kg$SUPER == "AFR"]) > 0) kg$PC1 <- -kg$PC1
mono <- pb$group %in% c("CD14+ Mono", "FCGR3A+ Mono")
if (mean(pb$PC1[mono]) > 0) pb$PC1 <- -pb$PC1

# ── 4. 그림 11 ─────────────────────────────────────────────
base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.title = element_blank(), plot.subtitle = element_text(colour = "grey30"))

p_kg <- ggplot(kg, aes(PC1, PC2, colour = group)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_colour_manual(values = super_colours, labels = super_labels) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 2)) +
  labs(title = "사람 집단: 1000 Genomes",
       subtitle = "점 하나 = 사람 1명 (2,491명), 변수 = SNP 8,551개\nPC1: 아프리카 vs 나머지") +
  base_theme

p_pb <- ggplot(pb, aes(PC1, PC2, colour = group)) +
  geom_point(size = 0.5, alpha = 0.6) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 3)) +
  labs(title = "세포 집단: PBMC 단일세포",
       subtitle = "점 하나 = 세포 1개 (2,511개), 변수 = 유전자 2,000개\nPC1: 단핵구 vs 림프구") +
  base_theme

p <- (p_kg | p_pb) +
  plot_annotation(title = "같은 PCA, 다른 해상도",
                  caption = "PC의 부호는 임의라서, 가장 먼저 갈라지는 쪽(아프리카, 단핵구)을 왼쪽에 두었다.",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
ggsave("results/figures/fig11_same_pca_two_scales.png", p, width = 12, height = 6.5, dpi = 200, bg = "white")
cat("저장: results/figures/fig11_same_pca_two_scales.png\n")
