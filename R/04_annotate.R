# 04_annotate.R — 마커 유전자 찾기 + 세포 유형 이름 붙이기
# 입력: work/pbmc_cluster.rds (03_cluster.R 결과, resolution 0.5 → 9개 군집)
# 출력: 그림 7–8, 표 6–7, work/pbmc_annot.rds
#
# 실행 방법: RStudio에서 한 덩어리씩 (Ctrl+Enter) 실행한다.
#   Part A를 끝까지 돌리고 표 6, 그림 7a–7d를 본 뒤,
#   Part B의 "직접 정할 값"(군집별 세포 이름)을 채우고 Part B를 실행한다.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

pbmc <- readRDS("work/pbmc_cluster.rds")
cat("세포:", ncol(pbmc), "/ resolution:", pbmc@misc$resolution, "/ 군집:", nlevels(Idents(pbmc)), "개\n")

# ═══ Part A. 마커 유전자 찾기 ══════════════════════════════

# ── 1. 군집마다 마커 찾기 (데이터가 알려 주는 마커) ─────────
# 군집 하나 vs 나머지 모든 세포를 비교해서, 그 군집에서 특히 높게 켜진 유전자를 찾는다.
#   only.pos = TRUE      : 그 군집에서 "높은" 유전자만 (낮은 유전자는 빼고)
#   min.pct = 0.25       : 그 군집 세포의 25% 이상에서 검출된 유전자만 → 몇 세포에서만 우연히 잡힌 유전자 제외
#   logfc.threshold = 0.25 : 나머지 세포보다 어느 정도 이상 높은 유전자만
# 결과 열:
#   avg_log2FC : 그 군집이 나머지보다 몇 배 높은가 (log2: 1이면 2배, 2이면 4배)
#   pct.1      : 그 군집 세포 중 이 유전자가 검출된 비율
#   pct.2      : 나머지 세포 중 이 유전자가 검출된 비율
#   p_val_adj  : 유전자 1만여 개를 동시에 검정했으므로 보정한 p-value (Bonferroni)
# 1~3분 걸린다.
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write_csv(markers, "results/tables/table6_all_markers.csv")

# 표 6: 군집마다 avg_log2FC가 큰 상위 8개 유전자
top_markers <- markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 8, with_ties = FALSE) |>
  ungroup()
top_markers |>
  group_by(cluster) |>
  summarise(상위_마커 = paste(gene, collapse = ", ")) |>
  print(n = Inf, width = Inf)
write_csv(top_markers, "results/tables/table6_top_markers.csv")

# 그림 7a: 군집별 상위 마커 히트맵 (세포를 군집 순서로 늘어놓음)
#   대각선 블록이 선명하면 군집마다 자기만의 마커가 있다는 뜻
top5 <- markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 5, with_ties = FALSE) |> pull(gene) |> unique()
p_heat <- DoHeatmap(pbmc, features = top5, size = 3) + NoLegend()
ggsave("results/figures/fig7a_marker_heatmap.png", p_heat, width = 10, height = 10, dpi = 150)

# ── 2. 알려진 마커로 확인 (교과서가 알려 주는 마커) ────────
# PBMC 세포 종류마다 널리 쓰이는 마커. 데이터가 찾은 마커(표 6)와 교차 확인한다.
known <- c(
  "CD3E",                 # T세포 전체
  "IL7R", "CCR7",         # Naive CD4 T (IL7R + CCR7)
  "S100A4",               # Memory CD4 T (IL7R + S100A4)
  "CD8A",                 # CD8 T
  "GNLY", "NKG7",         # NK (+ 세포독성 T)
  "MS4A1", "CD79A",       # B
  "CD14", "LYZ",          # CD14+ 단핵구
  "FCGR3A", "MS4A7",      # FCGR3A+ 단핵구
  "FCER1A", "CST3",       # 수지상세포 (CST3는 단핵구에도)
  "PPBP"                  # 혈소판
)

# 그림 7b: Dot plot
#   점 크기 = 그 군집에서 이 유전자가 검출된 세포 비율 (pct)
#   점 색   = 그 군집의 평균 발현 (진할수록 높음, 유전자마다 스케일링)
p_dot <- DotPlot(pbmc, features = known) + RotatedAxis() +
  labs(x = NULL, y = "군집") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/figures/fig7b_known_markers_dot.png", p_dot, width = 10, height = 5, dpi = 150)

# 그림 7c: 같은 마커를 바이올린 그림으로 (군집별 발현 분포)
p_vln <- VlnPlot(pbmc, features = known, stack = TRUE, flip = TRUE, fill.by = "ident") + NoLegend()
ggsave("results/figures/fig7c_known_markers_violin.png", p_vln, width = 8, height = 10, dpi = 150)

# ── 3. 앞에서 세운 가설 확인 ───────────────────────────────
# (1) Day 1+ 리보솜 비율: 림프구는 높고, 단핵구는 낮고, 혈소판은 아주 낮을 것
# (2) Day 3 5번 군집: RNA가 많은데 FCGR3A+ 단핵구인가, doublet인가
#     → 자기만의 마커(FCGR3A, MS4A7)가 있으면 진짜 세포, T세포 마커까지 섞여 있으면 doublet 의심
# (3) Day 2 PC9 인터페론 반응: 특정 군집에 몰려 있나, 여러 군집에 흩어져 있나
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
ifn_genes <- c("ISG15", "IFIT1", "IFIT3", "MX1", "IFI6")
pbmc <- AddModuleScore(pbmc, features = list(ifn_genes), name = "IFN_score", seed = 42)

hyp_tab <- pbmc@meta.data |>
  group_by(cluster) |>
  summarise(세포수 = n(),
            리보솜_중앙값 = round(median(percent.ribo), 1),
            인터페론점수_중앙값 = round(median(IFN_score1), 3),
            인터페론점수_상위5퍼센트_세포수 = sum(IFN_score1 > quantile(pbmc$IFN_score1, 0.95)))
print(hyp_tab)
write_csv(hyp_tab, "results/tables/table7_hypotheses.csv")

# 그림 7d: 가설 확인용 — 리보솜 비율(군집별), 인터페론 점수(UMAP)
p_ribo <- VlnPlot(pbmc, features = "percent.ribo", pt.size = 0) + NoLegend() +
  labs(x = "군집", y = "리보솜 유전자 비율 (%)", title = "리보솜 비율")
p_ifn  <- FeaturePlot(pbmc, features = "IFN_score1") + labs(title = "인터페론 반응 점수") + theme_bw()
ggsave("results/figures/fig7d_hypotheses.png", wrap_plots(p_ribo, p_ifn, ncol = 2), width = 12, height = 5, dpi = 150)

# ═══ Part B. 직접 정할 값 ══════════════════════════════════
# 표 6, 그림 7a–7d를 보고 군집마다 세포 이름을 정한다. 이유는 notes.md "Day 4"에 적는다.
# 쓸 수 있는 이름 예: "Naive CD4 T", "Memory CD4 T", "CD8 T", "NK", "B",
#                     "CD14+ Mono", "FCGR3A+ Mono", "DC", "Platelet"
# 두 군집이 같은 세포라고 판단하면 같은 이름을 넣으면 된다 (자동으로 합쳐진다).
# 판단이 어려우면 "Unknown"이라고 넣어도 된다.

# ▶ 직접 정할 값
cell_names <- c(
  "0" ="Naive CD4 T",
  "1" = "CD14+ Mono",
  "2" = "Memory CD4 T",
  "3" = "B",
  "4" = "CD8 T",
  "5" = "FCGR3A+ Mono",
  "6" = "NK",
  "7" = "DC",
  "8" = "Platelet"
)

stopifnot("cell_names의 NA를 모두 이름으로 바꾸세요" = !any(is.na(cell_names)))
pbmc <- RenameIdents(pbmc, cell_names)
pbmc$cell_type <- Idents(pbmc)

# 그림 8: 이름 붙인 UMAP
p_named <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4) +
  NoLegend() + labs(title = "PBMC 세포 유형") + theme_bw()
ggsave("results/figures/fig8_umap_celltypes.png", p_named, width = 7, height = 6, dpi = 150)

# 표 8: 세포 유형별 세포 수와 비율 → 건강한 사람 PBMC의 대략적 비율과 비교
#   T세포 약 45–70%, 단핵구 약 10–20%, NK 약 5–15%, B 약 5–10%, 수지상세포 약 1–2%
type_tab <- pbmc@meta.data |>
  dplyr::count(cell_type, name = "세포수") |>
  mutate(비율 = round(100 * 세포수 / sum(세포수), 1)) |>
  arrange(desc(세포수))
print(type_tab)
write_csv(type_tab, "results/tables/table8_cell_types.csv")

saveRDS(pbmc, "work/pbmc_annot.rds")
cat("저장: work/pbmc_annot.rds (", nlevels(Idents(pbmc)), "개 세포 유형)\n")
