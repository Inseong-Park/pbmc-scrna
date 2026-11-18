# 03_cluster.R — 군집화 + UMAP
# 입력: work/pbmc_pca.rds (02_pca.R 결과, PC 개수가 함께 저장되어 있음)
# 출력: 그림 5–6, 표 4–5, work/pbmc_cluster.rds
#
# 실행 방법: RStudio에서 한 덩어리씩 (Ctrl+Enter) 실행한다.
#   Part A를 끝까지 돌리고 그림 5a–5c와 표 4를 본 뒤,
#   Part B의 "직접 정할 값"(resolution)을 채우고 Part B를 실행한다.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

pbmc  <- readRDS("work/pbmc_pca.rds")
n_pcs <- pbmc@misc$n_pcs
cat("세포:", ncol(pbmc), "/ 사용할 PC:", n_pcs, "개\n")

# ═══ Part A. 군집화와 UMAP ═════════════════════════════════

# ── 1. 이웃 그래프 만들기 ───────────────────────────────────
# PC 1~n_pcs 공간에서 세포마다 가장 가까운 세포 20개(이웃)를 찾고,
# 이웃을 많이 공유하는 세포끼리 강하게 연결한 그래프(SNN)를 만든다.
# → "비슷한 세포끼리 손을 잡은 네트워크"
pbmc <- FindNeighbors(pbmc, dims = 1:n_pcs, verbose = FALSE)

# ── 2. 여러 resolution으로 군집 나누기 ─────────────────────
# 그래프에서 서로 촘촘히 연결된 무리를 찾는다 (Louvain 알고리즘).
# resolution: 무리를 얼마나 잘게 나눌지 정하는 값
#   작게 (0.3) → 큰 무리 몇 개 / 크게 (1.2) → 작은 무리 여러 개
# 한 값만 쓰지 않고 네 값으로 나눠 보고 비교한다 (1000 Genomes 전처리 비교와 같은 생각)
res_try <- c(0.3, 0.5, 0.8, 1.2)
for (r in res_try) {
  pbmc <- FindClusters(pbmc, resolution = r, random.seed = 42, verbose = FALSE)
}
res_cols <- paste0("RNA_snn_res.", res_try)

# 표 4: resolution별 군집 수와 가장 작은 군집의 세포 수
res_tab <- tibble(resolution = res_try,
                  군집수 = map_int(res_cols, ~ n_distinct(pbmc@meta.data[[.x]])),
                  가장_작은_군집 = map_int(res_cols, ~ min(table(pbmc@meta.data[[.x]]))),
                  가장_큰_군집 = map_int(res_cols, ~ max(table(pbmc@meta.data[[.x]]))))
print(res_tab)
write_csv(res_tab, "results/tables/table4_resolution.csv")

# ── 3. UMAP ────────────────────────────────────────────────
# PC 1~n_pcs를 2차원 그림으로 펼친다. 가까운 세포끼리 가깝게 놓는 데 집중한다.
# 주의: UMAP에서 군집 "사이의 거리"와 군집의 "크기"는 그대로 해석하면 안 된다.
#       같은 군집 안에서 누가 이웃인지는 믿을 만하지만, 멀리 떨어진 두 군집이 얼마나 다른지는 PCA로 봐야 한다.
pbmc <- RunUMAP(pbmc, dims = 1:n_pcs, seed.use = 42, verbose = FALSE)

# 그림 5a: resolution 네 개의 UMAP을 나란히
p_res <- map2(res_cols, res_try, function(col, r) {
  DimPlot(pbmc, reduction = "umap", group.by = col, label = TRUE, label.size = 3) +
    NoLegend() + labs(title = paste0("resolution ", r, " (", n_distinct(pbmc@meta.data[[col]]), "개)")) +
    theme_bw()
})
ggsave("results/figures/fig5a_umap_resolutions.png", wrap_plots(p_res, ncol = 2),
       width = 10, height = 10, dpi = 150)

# 그림 5b: resolution이 커질 때 군집이 어떻게 쪼개지나 (0.3 → 0.5 → 0.8 → 1.2)
#   각 칸: 낮은 resolution의 군집(행)에 속한 세포가 높은 resolution에서 어느 군집(열)으로 갔나
#   한 행이 한 칸에 몰려 있으면 그대로 유지, 여러 칸으로 나뉘면 쪼개진 것
split_plot <- function(a, b) {
  tab <- as.data.frame(table(from = pbmc@meta.data[[a]], to = pbmc@meta.data[[b]])) |> filter(Freq > 0)
  ggplot(tab, aes(to, from, fill = Freq)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = Freq), size = 2.5) +
    scale_fill_gradient(low = "grey95", high = "#2a78d6", name = "세포 수") +
    labs(x = paste("resolution", sub("RNA_snn_res.", "", b), "군집"),
         y = paste("resolution", sub("RNA_snn_res.", "", a), "군집")) +
    theme_bw()
}
p_split <- split_plot(res_cols[1], res_cols[2]) / split_plot(res_cols[2], res_cols[3]) / split_plot(res_cols[3], res_cols[4])
ggsave("results/figures/fig5b_cluster_splits.png", p_split, width = 9, height = 14, dpi = 150)

# 그림 5c: UMAP 위에 QC 지표와 대표 마커를 칠해 본다
#   - QC 지표가 한 군집에 몰려 있으면 "세포 종류"가 아니라 "품질" 때문에 생긴 군집일 수 있다
#   - 마커는 Day 4에서 자세히 본다. 여기서는 큰 덩어리가 무엇인지 감만 잡는다
p_qc <- FeaturePlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) &
  theme_bw()
p_mk <- FeaturePlot(pbmc, features = c("CD3E", "MS4A1", "LYZ", "NKG7", "PPBP", "FCER1A"), ncol = 3) &
  theme_bw()
ggsave("results/figures/fig5c_umap_qc.png", p_qc, width = 12, height = 4, dpi = 150)
ggsave("results/figures/fig5d_umap_markers.png", p_mk, width = 12, height = 8, dpi = 150)

# 그림 6: 같은 세포를 PCA(PC1–PC2)와 UMAP으로 나란히
#   PCA는 세포 사이의 큰 차이(축)를, UMAP은 이웃 관계를 보여 준다
p_pca  <- DimPlot(pbmc, reduction = "pca",  group.by = res_cols[2]) + NoLegend() +
  labs(title = "PCA (PC1–PC2)") + theme_bw()
p_umap <- DimPlot(pbmc, reduction = "umap", group.by = res_cols[2], label = TRUE) + NoLegend() +
  labs(title = paste0("UMAP (PC 1–", n_pcs, ")")) + theme_bw()
ggsave("results/figures/fig6_pca_vs_umap.png", p_pca + p_umap, width = 10, height = 5, dpi = 150)

# ═══ Part B. 직접 정할 값 ══════════════════════════════════
# 그림 5a–5d와 표 4를 보고 정한다. 이유는 notes.md "Day 3 판단"에 적는다.
# res_try 중 하나(0.3, 0.5, 0.8, 1.2)를 고른다.

# ▶ 직접 정할 값
resolution <- 0.5

stopifnot("resolution을 먼저 정하세요 (0.3, 0.5, 0.8, 1.2 중 하나)" = resolution %in% res_try)
res_col <- paste0("RNA_snn_res.", resolution)
Idents(pbmc) <- res_col
pbmc$cluster <- Idents(pbmc)
pbmc@misc$resolution <- resolution

# 표 5: 고른 resolution의 군집별 세포 수와 QC 지표 중앙값
clu_tab <- pbmc@meta.data |>
  group_by(cluster) |>
  summarise(세포수 = n(),
            비율 = round(100 * n() / ncol(pbmc), 1),
            UMI_중앙값 = median(nCount_RNA),
            유전자_중앙값 = median(nFeature_RNA),
            MT_중앙값 = round(median(percent.mt), 2))
print(clu_tab)
write_csv(clu_tab, "results/tables/table5_clusters.csv")

saveRDS(pbmc, "work/pbmc_cluster.rds")
cat("저장: work/pbmc_cluster.rds (resolution", resolution, "→", nlevels(Idents(pbmc)), "개 군집)\n")
