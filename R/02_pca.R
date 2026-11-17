# 02_pca.R — 정규화 → 변이 큰 유전자 → PCA
# 입력: work/pbmc_qc.rds (01_qc.R 결과)
# 출력: work/pbmc_pca.rds, 그림 3–4, 표 3
#
# 실행 방법: RStudio에서 위에서부터 한 덩어리씩 (Ctrl+Enter) 실행한다.
#   Part A를 끝까지 돌리고 그림 4(Elbow plot)를 본 뒤,
#   Part B의 "직접 정할 값"을 채우고 Part B를 실행한다.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

pbmc <- readRDS("work/pbmc_qc.rds")

# ═══ Part A. PCA까지 ═══════════════════════════════════════

# ── 1. 정규화 ───────────────────────────────────────────────
# 세포마다 잡힌 RNA 총량이 다르다 (UMI 중앙값 2,196, 최대 1만 이상).
# 총량이 많은 세포는 모든 유전자 값이 커 보이므로, 세포마다 총합을 10,000으로 맞추고 log(1 + x)로 바꾼다.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ── 2. 변이가 큰 유전자 고르기 ──────────────────────────────
# 13,714개 유전자 중 대부분은 모든 세포에서 비슷하거나 거의 0 → 세포를 구분하는 정보가 없다.
# 평균 발현량에 비해 세포 간 차이가 큰 유전자 2,000개만 고른다.
# (1000 Genomes에서 MAF가 낮은 SNP를 뺀 것과 비슷한 이유: 정보가 있는 변수만 PCA에 넣는다)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
cat("변이가 가장 큰 유전자 10개:", paste(top10, collapse = ", "), "\n\n")

p_var <- LabelPoints(VariableFeaturePlot(pbmc), points = top10, repel = TRUE) +
  labs(x = "평균 발현량", y = "표준화한 분산") + theme_bw()
ggsave("results/figures/fig3_variable_genes.png", p_var, width = 7, height = 5, dpi = 150)

# ── 3. 스케일링 → PCA ───────────────────────────────────────
# 유전자마다 평균 0, 분산 1로 맞춘다 → 원래 발현량이 큰 유전자가 PCA를 독차지하지 않게
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 50, verbose = FALSE)

# ── 4. 몇 번째 PC까지 쓸지 판단할 자료 ─────────────────────
# (a) Elbow plot: PC마다 설명하는 분산. 꺾이는 지점 이후는 잡음에 가깝다
sdev    <- pbmc[["pca"]]@stdev
var_pct <- 100 * sdev^2 / sum(sdev^2)   # 50개 PC 안에서의 비율
elbow <- tibble(PC = 1:50, var_pct = var_pct, cum_pct = cumsum(var_pct))
write_csv(elbow, "results/tables/table3_pca_variance.csv")

p_elbow <- ggplot(elbow, aes(PC, var_pct)) +
  geom_line(colour = "grey50") + geom_point(size = 1.5) +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  labs(x = "주성분(PC)", y = "설명 분산 (%, PC 50개 안에서)") + theme_bw()
ggsave("results/figures/fig4a_elbow.png", p_elbow, width = 7, height = 4, dpi = 150)

# (b) PC별 로딩 상위 유전자: 각 PC가 어떤 세포 차이를 잡는지
#     알려진 마커(예: LYZ=단핵구, MS4A1=B세포, NKG7=NK세포)가 보이면 생물학적 신호, 아니면 잡음일 수 있다
print(pbmc[["pca"]], dims = 1:15, nfeatures = 5)

# (c) PC1–PC2 산점도 (1000 Genomes 그림 3과 같은 그림)
p_pca <- DimPlot(pbmc, reduction = "pca", dims = c(1, 2)) + NoLegend() + theme_bw()
ggsave("results/figures/fig4b_pca_pc1_pc2.png", p_pca, width = 5, height = 5, dpi = 150)

# (d) PC별 히트맵: 세포를 PC 점수 순서로 놓고 상위 유전자 발현을 본다
#     뒤쪽 PC로 갈수록 무늬가 흐려지면 그 PC는 잡음에 가깝다
png("results/figures/fig4c_pc_heatmaps.png", width = 1400, height = 1800, res = 120)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE, ncol = 3)
dev.off()

# ═══ Part B. 직접 정할 값 ══════════════════════════════════
# 그림 4a–4c를 보고 정한다. 이유는 notes.md "Day 2 판단"에 적는다.

# ▶ 직접 정할 값
n_pcs <- 9   # 다음 단계(군집화, UMAP)에 쓸 PC 개수

stopifnot("n_pcs를 먼저 정하세요" = !is.na(n_pcs))
pbmc@misc$n_pcs <- n_pcs
saveRDS(pbmc, "work/pbmc_pca.rds")
cat("저장: work/pbmc_pca.rds (PC", n_pcs, "개 사용)\n")
