# 01_qc.R — 세포 품질 관리(QC)
# 하는 일:
#   1. 10x 행렬을 Seurat 객체로 만든다 (너무 드문 유전자, 유전자가 너무 적은 바코드는 처음부터 뺀다)
#   2. 세포마다 QC 지표 3개를 계산한다
#        nFeature_RNA : 검출된 유전자 수 (0보다 큰 유전자 개수)
#        nCount_RNA   : 전체 UMI 수 (= 잡힌 RNA 분자 수)
#        percent.mt   : 미토콘드리아 유전자(MT-)에서 나온 UMI 비율
#   3. 분포를 그림으로 보고 기준을 정해 세포를 거른다
#   4. 거른 객체를 work/pbmc_qc.rds로 저장한다 (다음 단계의 입력)
#
# 왜 이 지표인가:
#   - 죽거나 터진 세포: 세포질 RNA는 새어 나가고 미토콘드리아 안의 RNA만 남음 → percent.mt 높음, 유전자 수 적음
#   - 빈 방울(세포 없이 주변 RNA만 잡힘): 유전자 수가 매우 적음
#   - doublet(세포 두 개가 한 방울에): 유전자 수와 UMI가 비정상적으로 많음

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

dir.create("work", showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ── 1. Seurat 객체 만들기 ───────────────────────────────────
counts <- Read10X("data/filtered_gene_bc_matrices/hg19/")

# min.cells = 3    : 세포 3개 미만에서만 검출된 유전자는 뺀다 (정보가 거의 없음)
# min.features = 200 : 유전자가 200개 미만인 바코드는 뺀다 (빈 방울일 가능성)
pbmc <- CreateSeuratObject(counts, project = "pbmc3k", min.cells = 3, min.features = 200)
cat("처음:", nrow(counts), "유전자 x", ncol(counts), "세포\n")
cat("객체 생성 후:", nrow(pbmc), "유전자 x", ncol(pbmc), "세포\n\n")

# ── 2. QC 지표 ──────────────────────────────────────────────
# hg19 유전자 이름에서 미토콘드리아 유전자는 "MT-"로 시작한다 (MT-CO1, MT-ND1 등 13개 단백질 + rRNA/tRNA)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
cat("MT- 유전자 수:", sum(grepl("^MT-", rownames(pbmc))), "\n\n")

qc <- pbmc@meta.data |> rownames_to_column("barcode")

# 분포 요약: 기준을 정하기 전에 숫자로 먼저 본다
summ <- qc |>
  summarise(across(c(nFeature_RNA, nCount_RNA, percent.mt),
                   list(min = min, q01 = ~quantile(.x, 0.01), median = median,
                        q99 = ~quantile(.x, 0.99), max = max))) |>
  pivot_longer(everything(), names_to = c("metric", "stat"), names_pattern = "(.*)_(min|q01|median|q99|max)") |>
  pivot_wider(names_from = stat, values_from = value)
print(summ, width = Inf)
cat("\n")

# ── 3. 기준 정하기 ──────────────────────────────────────────
# 기준 (판단 근거는 notes.md에 기록)
#   nFeature_RNA > 200  : 빈 방울·터진 세포 제외 (객체 생성 때 이미 적용)
#   nFeature_RNA < 2500 : doublet 의심 제외
#   percent.mt   < 5    : 죽어 가는 세포 제외
min_feat <- 200
max_feat <- 2500
max_mt   <- 5

qc <- qc |>
  mutate(fail = case_when(
    nFeature_RNA >= max_feat ~ "유전자 많음 (doublet 의심)",
    percent.mt   >= max_mt   ~ "미토콘드리아 비율 높음",
    TRUE ~ "통과"
  ))

fail_tab <- qc |> count(fail) |> mutate(pct = round(100 * n / sum(n), 1))
print(fail_tab)
write_csv(fail_tab, "results/tables/table1_qc_filter.csv")

# 기준 하나만 바꿨을 때 몇 개가 빠지는지 (민감도 확인)
sens <- expand_grid(max_feat = c(2000, 2500, 3000), max_mt = c(5, 10, 20)) |>
  rowwise() |>
  mutate(kept = sum(qc$nFeature_RNA < max_feat & qc$percent.mt < max_mt)) |>
  ungroup()
cat("\n기준별 남는 세포 수\n")
print(pivot_wider(sens, names_from = max_mt, values_from = kept, names_prefix = "mt<"))
write_csv(sens, "results/tables/table2_qc_sensitivity.csv")

# ── 그림 ────────────────────────────────────────────────────
fill_cols <- c("통과" = "grey60", "유전자 많음 (doublet 의심)" = "#D55E00",
               "미토콘드리아 비율 높음" = "#0072B2")

violin <- function(var, cut, lab) {
  ggplot(qc, aes(x = "", y = .data[[var]])) +
    geom_violin(fill = "grey90", colour = "grey50") +
    geom_jitter(aes(colour = fail), width = 0.25, size = 0.3, alpha = 0.6) +
    geom_hline(yintercept = cut, linetype = "dashed", colour = "red") +
    scale_colour_manual(values = fill_cols, guide = "none") +
    labs(x = NULL, y = lab) +
    theme_bw()
}
p1 <- violin("nFeature_RNA", c(min_feat, max_feat), "검출 유전자 수")
p2 <- violin("nCount_RNA", NULL, "UMI 수")
p3 <- violin("percent.mt", max_mt, "미토콘드리아 비율 (%)")
ggsave("results/figures/fig1_qc_violin.png", p1 + p2 + p3, width = 9, height = 4, dpi = 150)

# UMI 수와 유전자 수는 함께 움직여야 정상. 벗어나는 점과 MT 비율의 관계를 본다
s1 <- ggplot(qc, aes(nCount_RNA, nFeature_RNA, colour = fail)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = max_feat, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = fill_cols, name = NULL) +
  labs(x = "UMI 수", y = "검출 유전자 수") + theme_bw()
s2 <- ggplot(qc, aes(nCount_RNA, percent.mt, colour = fail)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = max_mt, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = fill_cols, name = NULL) +
  labs(x = "UMI 수", y = "미토콘드리아 비율 (%)") + theme_bw()
ggsave("results/figures/fig2_qc_scatter.png", s1 + s2 + plot_layout(guides = "collect") &
         theme(legend.position = "bottom"), width = 9, height = 4.5, dpi = 150)

# ── 4. 거르고 저장 ──────────────────────────────────────────
pbmc <- subset(pbmc, subset = nFeature_RNA > min_feat & nFeature_RNA < max_feat & percent.mt < max_mt)
cat("\nQC 후:", nrow(pbmc), "유전자 x", ncol(pbmc), "세포\n")
saveRDS(pbmc, "work/pbmc_qc.rds")
