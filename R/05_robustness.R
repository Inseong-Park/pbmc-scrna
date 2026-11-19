# 05_robustness.R — 기준을 바꿔도 세포 유형 결과가 유지되나? (민감도 분석)
# 입력: 원본 행렬, work/doublet_scores.csv (01b), work/pbmc_annot.rds (04, 기준 결과)
# 출력: 그림 9a–9d, 10a–10c, 표 9–12 (9b, 9c, 12b, 12c 포함), work/robustness_runs.rds
#
# 왜 하나:
#   Day 1–4에서 MT 5%, doublet 제거, PC 9개, resolution 0.5를 "직접" 정했다.
#   정한 값이 결론(9개 세포 유형, 비율)을 좌우한다면 결론이 약하고,
#   값을 바꿔도 같은 세포 유형이 나오면 결론이 튼튼하다.
#   1000 Genomes에서 MAF·LD pruning·친척 제외를 바꿔 가며 PCA를 다시 본 것과 같은 생각이다.
#
# 바꿔 볼 기준 (한 번에 하나만 바꾼다 → 차이가 나면 그 기준 때문이라고 말할 수 있다)
#   mt10     : MT 5% → 10%      (notes Day 1: 빠진 세포 중 정상 세포가 섞였을 수 있다)
#   dbl_keep : doublet 제거 → 남김 (notes Day 1+: 남기면 가짜 중간 군집이 생길까?)
#   pc7/pc10 : PC 9 → 7 / 10     (notes Day 2: PC8(DC), PC9(인터페론)를 빼면?)
#   res08    : resolution 0.5 → 0.8 (notes Day 3–4: Naive CD4 T가 세 덩어리로 쪼개짐)
#
# 비교 방법:
#   새로 나온 군집에는 이름이 없으므로, 기준 결과에도 있는 세포(공통 세포)의 이름을 빌린다.
#   새 군집마다 "공통 세포 중 가장 많은 기준 이름"을 붙이고(다수결), 그 비율을 순도로 본다.
#   기준에는 없던 세포(MT 5–10% 세포, doublet)는 자기가 들어간 군집의 이름을 받는다.
#
# 실행 방법: Session → Restart R (Ctrl+Shift+F10) 후 시작.
#   Part A(5–10분 걸림)를 돌리면 work/robustness_runs.rds가 저장된다.
#   다시 열 때는 Part A 대신 맨 아래 "다시 불러오기" 한 줄만 실행하고 Part B부터 진행한다.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

annot <- readRDS("work/pbmc_annot.rds")
base_lab <- setNames(as.character(annot$cell_type), colnames(annot))   # 바코드 → 기준 세포 이름
type_levels <- levels(annot$cell_type)
cat("기준 결과:", length(base_lab), "세포,", length(type_levels), "개 세포 유형\n")

# ═══ Part A. 기준을 하나씩 바꿔서 처음부터 다시 돌리기 ══════

# ── 1. 원본에서 다시 시작 (01_qc.R와 같은 방식) ─────────────
counts <- Read10X("data/filtered_gene_bc_matrices/hg19/")
raw    <- CreateSeuratObject(counts, project = "pbmc3k", min.cells = 3, min.features = 200)
raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^MT-")
dbl <- read_csv("work/doublet_scores.csv", show_col_types = FALSE)
raw$dbl_class <- dbl$dbl_class[match(colnames(raw), dbl$barcode)]

# ── 2. 바꿔 볼 설정 ────────────────────────────────────────
# base는 Day 1–4와 똑같은 설정이다. 다시 돌려서 기준 결과와 똑같이 나오는지 먼저 확인한다
# (똑같이 안 나오면 나머지 비교도 믿을 수 없다).
settings <- tribble(
  ~name,      ~label,              ~max_mt, ~remove_dbl, ~n_pcs, ~res,
  "base",     "기준",                    5,  TRUE,         9,  0.5,
  "mt10",     "MT 10%",                 10,  TRUE,         9,  0.5,
  "dbl_keep", "doublet 남김",            5,  FALSE,        9,  0.5,
  "pc7",      "PC 7개",                  5,  TRUE,         7,  0.5,
  "pc10",     "PC 10개",                 5,  TRUE,        10,  0.5,
  "res08",    "resolution 0.8",          5,  TRUE,         9,  0.8
)

# ── 3. QC → 정규화 → PCA → 군집 → UMAP을 한 번에 ───────────
# 02_pca.R, 03_cluster.R과 같은 함수·값·seed를 쓴다.
# ScaleData는 변이 큰 유전자 2,000개만 한다. PCA에는 이 유전자만 쓰이므로 결과는 같고 더 빠르다.
run_pipeline <- function(max_mt, remove_dbl, n_pcs, res) {
  obj <- subset(raw, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < max_mt)
  if (remove_dbl) obj <- subset(obj, subset = dbl_class == "singlet")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:n_pcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res, random.seed = 42, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:n_pcs, seed.use = 42, verbose = FALSE)
  obj
}

runs <- list()
for (i in seq_len(nrow(settings))) {
  s <- settings[i, ]
  cat("실행 중:", s$label, "...\n")
  runs[[s$name]] <- run_pipeline(s$max_mt, s$remove_dbl, s$n_pcs, s$res)
}
saveRDS(runs, "work/robustness_runs.rds")
cat("저장: work/robustness_runs.rds\n")

# ═══ Part B. 기준 결과와 비교 ══════════════════════════════

# ── 4. 새 군집에 이름 빌려 오기 (다수결) ────────────────────
transfer_labels <- function(obj) {
  md <- obj@meta.data |>
    rownames_to_column("barcode") |>
    transmute(barcode, cluster = as.character(seurat_clusters), base = base_lab[barcode],
              dbl_class, percent.mt)
  clu <- md |>
    filter(!is.na(base)) |>
    dplyr::count(cluster, base) |>
    group_by(cluster) |>
    mutate(순도 = n / sum(n)) |>
    slice_max(n, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(cluster, label = base, 순도 = round(순도, 3))
  md |> left_join(clu, by = "cluster") |>
    mutate(label = replace_na(label, "기준에 없는 군집"))
}

# ARI (Adjusted Rand Index): 두 가지 나누기가 얼마나 같은지 나타내는 값
#   세포 두 개를 아무렇게나 골랐을 때 "두 방법 모두 같은 무리" 또는 "두 방법 모두 다른 무리"로
#   일치하는 정도를 세고, 우연히 일치할 만큼을 뺀 값이다. 1 = 완전히 같음, 0 = 우연 수준.
#   주의: 군집이 더 잘게 쪼개지면(res 0.8) 같은 세포 유형이 둘로 나뉜 것도 "불일치"로 세서 낮아진다.
#   그래서 "이름 일치율"(공통 세포 중 빌려 온 이름 = 원래 이름인 비율)도 함께 본다.
ari <- function(a, b) {
  tab <- table(a, b)
  c2 <- function(x) sum(choose(x, 2))
  s_ab <- c2(tab); s_a <- c2(rowSums(tab)); s_b <- c2(colSums(tab))
  expected <- s_a * s_b / choose(sum(tab), 2)
  (s_ab - expected) / ((s_a + s_b) / 2 - expected)
}

res_md <- imap(runs, ~ transfer_labels(.x) |> mutate(setting = .y)) |> list_rbind()

# 표 9: 설정별 요약
summary_tab <- res_md |>
  group_by(setting) |>
  summarise(세포수 = n(),
            공통_세포 = sum(!is.na(base)),
            기준에_없던_세포 = sum(is.na(base)),
            군집수 = n_distinct(cluster),
            ARI = round(ari(cluster[!is.na(base)], base[!is.na(base)]), 3),
            이름_일치율 = round(mean(label[!is.na(base)] == base[!is.na(base)]), 3),
            가장_낮은_순도 = min(순도, na.rm = TRUE)) |>
  left_join(select(settings, setting = name, 설정 = label), by = "setting") |>
  relocate(설정) |>
  arrange(match(setting, settings$name))
print(summary_tab, width = Inf)
write_csv(summary_tab, "results/tables/table9_robustness_summary.csv")
cat("\n기준 재현 확인 (base의 ARI와 이름 일치율이 1이어야 한다):",
    summary_tab$ARI[summary_tab$setting == "base"], "/",
    summary_tab$이름_일치율[summary_tab$setting == "base"], "\n\n")

# 표 9b: 설정별 군집 → 이름 대응표 (군집마다 세포 수, 빌려 온 이름, 순도)
cluster_tab <- res_md |>
  group_by(setting, cluster, label, 순도) |>
  summarise(세포수 = n(), 기준에_없던_세포 = sum(is.na(base)),
            doublet = sum(dbl_class == "doublet"), .groups = "drop") |>
  arrange(match(setting, settings$name), as.integer(cluster))
print(cluster_tab, n = Inf, width = Inf)
write_csv(cluster_tab, "results/tables/table9b_robustness_clusters.csv")

# 표 9c: 설정별 세포 유형 비율 (%) — 비율이 기준에 따라 얼마나 흔들리나
prop_tab <- res_md |>
  dplyr::count(setting, label) |>
  group_by(setting) |>
  mutate(비율 = round(100 * n / sum(n), 1)) |>
  ungroup() |>
  select(-n) |>
  pivot_wider(names_from = setting, values_from = 비율, values_fill = 0) |>
  select(세포유형 = label, any_of(settings$name)) |>
  arrange(match(세포유형, type_levels))
print(prop_tab, width = Inf)
write_csv(prop_tab, "results/tables/table9c_robustness_proportions.csv")

# ── 5. 기준에 없던 세포는 어디로 갔나 ───────────────────────
# mt10: MT 5–10% 세포 → 죽어 가는 세포라면 한쪽에 따로 뭉치거나 마커가 흐릿할 것이고,
#       정상 세포라면 기존 세포 유형 군집에 고르게(또는 대사가 활발한 단핵구 쪽에) 들어갈 것이다
# dbl_keep: doublet → 가짜 중간 군집을 만드나, 기존 군집에 흩어지나
extra_tab <- res_md |>
  filter(setting %in% c("mt10", "dbl_keep"), is.na(base)) |>
  dplyr::count(setting, label, name = "세포수") |>
  group_by(setting) |>
  mutate(비율 = round(100 * 세포수 / sum(세포수), 1)) |>
  ungroup() |>
  arrange(setting, desc(세포수))
print(extra_tab, n = Inf)
write_csv(extra_tab, "results/tables/table10_extra_cells.csv")

# ── 그림 ───────────────────────────────────────────────────
type_cols <- setNames(scales::hue_pal()(length(type_levels)), type_levels)
type_cols["기준에 없는 군집"] <- "black"

# 그림 9a: 설정별 UMAP, 빌려 온 이름으로 색칠 (UMAP 모양은 설정마다 다르게 그려진다)
p_umaps <- imap(runs, function(obj, nm) {
  md <- filter(res_md, setting == nm)
  obj$label <- md$label[match(colnames(obj), md$barcode)]
  DimPlot(obj, reduction = "umap", group.by = "label", cols = type_cols, pt.size = 0.3,
          label = TRUE, repel = TRUE, label.size = 3) +
    labs(title = settings$label[settings$name == nm]) + theme_bw() + NoLegend()
})
ggsave("results/figures/fig9a_robustness_umaps.png", wrap_plots(p_umaps, ncol = 3),
       width = 15, height = 10, dpi = 150)

# 그림 9b: 기준 세포 유형(행)의 세포가 새 군집(열)에 어떻게 나뉘었나 (행마다 합 = 100%)
#   한 행에 진한 칸이 하나뿐이면 그 세포 유형이 그대로 유지된 것
#   한 행이 여러 칸으로 나뉘면 쪼개진 것, 한 열에 여러 행이 진하면 합쳐진 것
map_long <- res_md |>
  filter(!is.na(base)) |>
  dplyr::count(setting, base, cluster) |>
  group_by(setting, base) |>
  mutate(비율 = n / sum(n)) |>
  ungroup() |>
  mutate(base = factor(base, levels = rev(type_levels)),
         cluster = factor(cluster, levels = as.character(0:30)),
         setting = factor(settings$label[match(setting, settings$name)], levels = settings$label))
p_map <- ggplot(map_long, aes(cluster, base, fill = 비율)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(비율 >= 0.1, round(100 * 비율), "")), size = 2.5) +
  scale_fill_gradient(low = "grey95", high = "#0072B2", labels = scales::percent) +
  facet_wrap(~ setting, scales = "free_x", ncol = 3) +
  labs(x = "새 군집 번호", y = "기준 세포 유형") + theme_bw()
ggsave("results/figures/fig9b_type_vs_cluster.png", p_map, width = 15, height = 8, dpi = 150)

# 그림 9c: 기준에 없던 세포의 위치 (왼쪽: MT 5–10% 세포, 오른쪽: doublet)
extra_bc <- function(nm) res_md |> filter(setting == nm, is.na(base)) |> pull(barcode)
p_mt  <- DimPlot(runs$mt10, cells.highlight = extra_bc("mt10"), sizes.highlight = 0.6,
                 cols.highlight = "#D55E00", pt.size = 0.3) +
  labs(title = "MT 10%: 새로 들어온 MT 5–10% 세포") + theme_bw() + NoLegend()
p_dbl <- DimPlot(runs$dbl_keep, cells.highlight = extra_bc("dbl_keep"), sizes.highlight = 0.6,
                 cols.highlight = "#D55E00", pt.size = 0.3) +
  labs(title = "doublet 남김: scDblFinder doublet") + theme_bw() + NoLegend()
ggsave("results/figures/fig9c_extra_cells.png", p_mt + p_dbl, width = 12, height = 5.5, dpi = 150)

# ═══ Part C. resolution 0.8: T세포가 어떻게 쪼개지나 ════════
# Day 4에서 확신이 가장 낮았던 것: Naive CD4 T(0번)와 Memory CD4 T(2번)의 경계.
# res 0.8에서 T세포 군집이 몇 개로 나뉘는지, 나뉜 군집마다 자기만의 마커가 있는지 본다.
#   - 자기만의 마커가 뚜렷하면(예: CCR7·SELL vs S100A4·GZMK) 실제 하위 종류일 가능성
#   - 마커 차이가 작고 같은 유전자의 "정도 차이"뿐이면 연속적인 변화를 잘게 자른 것일 가능성
r08 <- runs$res08
md08 <- filter(res_md, setting == "res08")
r08$label <- md08$label[match(colnames(r08), md08$barcode)]
t_types <- c("Naive CD4 T", "Memory CD4 T", "CD8 T")
t_clusters <- md08 |> filter(label %in% t_types) |> distinct(cluster, label) |> arrange(label, as.integer(cluster))
print(t_clusters)

# 그림 9d: T세포 쪽 군집의 T세포 마커 dot plot
t_markers <- c("CD3E", "CCR7", "SELL", "LEF1", "TCF7",       # Naive
               "IL7R", "S100A4", "CD40LG", "GPR183",         # Memory CD4
               "CD8A", "CD8B", "GZMK", "CCL5", "NKG7",       # CD8 / 세포독성
               "ISG15", "IFIT3")                              # 인터페론 상태
r08_t <- subset(r08, idents = t_clusters$cluster)
p_dot <- DotPlot(r08_t, features = t_markers) + RotatedAxis() +
  labs(x = NULL, y = "res 0.8 군집") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/figures/fig9d_res08_tcells.png", p_dot, width = 10, height = 4.5, dpi = 150)

# 표 11: T세포 쪽 군집끼리만 비교한 마커 (다른 T세포 군집보다 높은 유전자)
t_mk <- FindAllMarkers(r08_t, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
t_top <- t_mk |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 8, with_ties = FALSE) |>
  ungroup() |>
  left_join(mutate(t_clusters, cluster = factor(cluster, levels = levels(t_mk$cluster))), by = "cluster")
t_top |>
  group_by(cluster, label) |>
  summarise(상위_마커 = paste(gene, collapse = ", "), .groups = "drop") |>
  print(n = Inf, width = Inf)
write_csv(t_top, "results/tables/table11_res08_tcell_markers.csv")

# ═══ Part D. res 0.8의 군집 6은 진짜 다른 세포인가 ═════════
# Part C에서 본 것: 군집 1과 5는 같은 Naive CD4 T로 보이고, 군집 6은 S100A4가 높아 Memory 쪽에 가까워 보였다.
# 하지만 군집 6은 순도가 낮고(CD8 T 16% 섞임), 표 11의 마커 차이가 작았다.
# 최종 결과는 res 0.5로 유지하고, 여기서는 군집 6이 무엇인지만 네 가지로 확인한다.
# (Part B, C를 먼저 실행해 r08, md08이 있어야 한다)

c6     <- "6"          # 그림 9b·9d에서 본 군집 번호 (seed가 같아서 다시 돌려도 같은 번호)
naive  <- c("1", "5")  # 같다고 판단한 Naive CD4 T 군집
memory <- "2"
cd8    <- "4"
r08$base <- md08$base[match(colnames(r08), md08$barcode)]   # Day 4에서 붙인 이름

# ── ① 군집 6 vs Naive, 군집 6 vs Memory 직접 비교 ───────────
# 차이 나는 유전자(p_val_adj < 0.05)가 몇 개인지 센다.
# 기준점으로 "같다고 판단한 1 vs 5"도 같이 센다 → 6이 진짜 다르면 1 vs 5보다 훨씬 많아야 한다.
de_compare <- function(a, b, nm) {
  FindMarkers(r08, ident.1 = a, ident.2 = b, logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE) |>
    rownames_to_column("gene") |>
    mutate(비교 = nm)
}
de <- bind_rows(
  de_compare(c6, naive,  "6 vs Naive(1+5)"),
  de_compare(c6, memory, "6 vs Memory(2)"),
  de_compare("1", "5",   "1 vs 5 (기준점: 같다고 본 두 군집)")
)
de_count <- de |>
  group_by(비교) |>
  summarise(유의한_유전자 = sum(p_val_adj < 0.05),
            앞_군집에서_높음 = sum(p_val_adj < 0.05 & avg_log2FC > 0),
            뒤_군집에서_높음 = sum(p_val_adj < 0.05 & avg_log2FC < 0),
            `2배_이상_차이` = sum(p_val_adj < 0.05 & abs(avg_log2FC) >= 1))
print(de_count, width = Inf)
write_csv(de_count, "results/tables/table12_c6_de_counts.csv")

de_top <- de |>
  filter(p_val_adj < 0.05) |>
  mutate(방향 = if_else(avg_log2FC > 0, "앞 군집에서 높음", "뒤 군집에서 높음")) |>
  group_by(비교, 방향) |>
  slice_max(abs(avg_log2FC), n = 10, with_ties = FALSE) |>
  summarise(유전자 = paste(gene, collapse = ", "), .groups = "drop")
print(de_top, n = Inf, width = Inf)
write_csv(de, "results/tables/table12_c6_de_all.csv")

# ── ② CD8 T가 섞인 세포를 빼고 다시 보기 ───────────────────
# 군집 6 세포를 Day 4 이름으로 나눈다: 기준에서 CD8 T였던 세포 vs 나머지.
# CD8 T를 빼도 S100A4·CCL5 등이 높으면 → 섞여 든 CD8 T 때문이 아니다.
# 메타데이터 열 이름은 영어로 쓴다 (Windows에서 한글 열 이름은 group.by가 찾지 못하는 경우가 있다)
r08$grp <- case_when(
  r08$seurat_clusters %in% naive ~ "Naive (1+5)",
  r08$seurat_clusters == c6 & r08$base == "CD8 T" ~ "6 중 기준 CD8 T",
  r08$seurat_clusters == c6 ~ "6 (CD8 T 제외)",
  r08$seurat_clusters == memory ~ "Memory CD4 (2)",
  r08$seurat_clusters == cd8 ~ "CD8 T (4)",
  TRUE ~ NA_character_
)
grp_levels <- c("Naive (1+5)", "6 (CD8 T 제외)", "6 중 기준 CD8 T", "Memory CD4 (2)", "CD8 T (4)")
r08_g <- subset(r08, cells = colnames(r08)[!is.na(r08$grp)])
r08_g$grp <- factor(r08_g$grp, levels = grp_levels)

cat("\n군집 6 세포의 Day 4 이름:\n")
print(table(r08$base[r08$seurat_clusters == c6]))
cat("\n그룹별 세포 수:\n")
print(table(r08_g$grp))

chk_genes <- c("CCR7", "SELL", "LEF1", "TCF7", "IL7R", "S100A4", "CCL5", "NKG7", "GZMK", "CD8A", "CD8B")
# 그림 10a: 그룹별 마커 (점 크기 = 검출 비율, 색 = 평균 발현)
p_grp <- DotPlot(r08_g, features = chk_genes, group.by = "grp") + RotatedAxis() +
  labs(x = NULL, y = NULL) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("results/figures/fig10a_c6_groups.png", p_grp, width = 9, height = 4, dpi = 150)

# 표 12b: 같은 내용을 숫자로 (검출 비율 %)
grp_pct <- FetchData(r08_g, vars = chk_genes, layer = "counts") |>
  mutate(그룹 = r08_g$grp) |>
  group_by(그룹) |>
  summarise(세포수 = n(), across(all_of(chk_genes), ~ round(100 * mean(.x > 0))))
print(grp_pct, width = Inf)
write_csv(grp_pct, "results/tables/table12b_c6_marker_pct.csv")

# ── ③ QC 지표 비교 ─────────────────────────────────────────
# 표 11에 MT-ND5, MT-ND1이 있었다. 군집 6이 "세포 상태"가 아니라 "품질" 차이로 생긴 것인지 본다.
# 품질 문제라면 MT 비율이 높거나, UMI·유전자 수가 낮을 것이다.
r08$percent.ribo <- PercentageFeatureSet(r08, pattern = "^RP[SL]")
t_cl <- c(naive, c6, memory, cd8)
qc_t <- r08@meta.data |>
  filter(seurat_clusters %in% t_cl) |>
  group_by(군집 = seurat_clusters) |>
  summarise(세포수 = n(),
            UMI_중앙값 = median(nCount_RNA),
            유전자_중앙값 = median(nFeature_RNA),
            MT_중앙값 = round(median(percent.mt), 2),
            리보솜_중앙값 = round(median(percent.ribo), 1))
print(qc_t)
write_csv(qc_t, "results/tables/table12c_c6_qc.csv")

r08_t2 <- subset(r08, idents = t_cl)
p_qc <- VlnPlot(r08_t2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                pt.size = 0, ncol = 4) &
  labs(x = "res 0.8 군집") & theme_bw() & NoLegend()
ggsave("results/figures/fig10b_c6_qc.png", p_qc, width = 14, height = 4, dpi = 150)

# ── ④ UMAP 위치 ───────────────────────────────────────────
# 왼쪽: T세포 군집만 색칠 (나머지 회색) → 6이 Naive와 Memory 사이에 있나, 따로 떨어져 있나
# 오른쪽: 같은 UMAP을 Day 4 이름으로 색칠 → 6 자리에 원래 어떤 이름의 세포가 있었나
r08$t_group <- if_else(r08$seurat_clusters %in% t_cl, as.character(r08$seurat_clusters), "기타")
t_cols <- c("1" = "#56B4E9", "5" = "#0072B2", "6" = "#D55E00", "2" = "#009E73", "4" = "#CC79A7", "기타" = "grey85")
p_u1 <- DimPlot(r08, group.by = "t_group", cols = t_cols, pt.size = 0.4, label = TRUE, label.size = 5,
                order = rev(names(t_cols))) +
  labs(title = "res 0.8 T세포 군집") + theme_bw() + NoLegend()
p_u2 <- DimPlot(r08, group.by = "base", cols = type_cols, pt.size = 0.4, label = TRUE, repel = TRUE, label.size = 3) +
  labs(title = "같은 UMAP, Day 4 이름") + theme_bw() + NoLegend()
ggsave("results/figures/fig10c_c6_umap.png", p_u1 + p_u2, width = 13, height = 6, dpi = 150)

# ═══ 직접 판단할 것 (notes.md "Day 5"에 적는다) ════════════
# 1. 표 9: base가 1 / 1로 재현되었나? 다른 설정의 ARI와 이름 일치율은 어느 정도인가?
# 2. 그림 9b·표 9b: 어떤 설정에서 세포 유형이 쪼개지거나 합쳐졌나? 순도가 낮은 군집이 있나?
#    PC 7개에서 DC(PC8)가 따로 군집을 이루나? PC 10개에서 새 군집이 생기나?
# 3. 표 9c: 세포 유형 비율이 설정에 따라 몇 %p나 흔들리나? Day 4 표 8의 결론이 유지되나?
# 4. 표 10·그림 9c: MT 5–10% 세포와 doublet은 어디로 갔나? 따로 뭉쳤나, 흩어졌나?
#    → Day 1의 MT 5%, Day 1+의 doublet 제거 판단을 유지할 것인가? (Day 5에서 둘 다 유지로 결정)
# 5. 그림 9d·표 11: res 0.8에서 쪼개진 T세포 군집은 실제 하위 종류인가, 연속적인 변화를 자른 것인가?
#    → 최종 결과를 res 0.5(9개)로 유지할 것인가? (Day 5에서 0.5 유지로 결정)
# 6. Part D: 군집 6은 진짜 다른 세포인가?
#    ① 표 12: 6 vs Naive, 6 vs Memory의 차이 유전자 수가 기준점(1 vs 5)보다 훨씬 많은가?
#    ② 그림 10a·표 12b: CD8 T를 빼도 S100A4·CCL5가 높은가?
#    ③ 표 12c·그림 10b: MT 비율, UMI, 유전자 수가 다른 T세포 군집과 비슷한가?
#    ④ 그림 10c: UMAP에서 Naive와 Memory 사이에 있나, 따로 떨어져 있나?
#    → Day 5 결론: 별개의 종류가 아니라 Naive–Memory 중간 상태 세포 + 경계의 CD8 T

# 다시 불러오기 (Part A를 다시 돌리지 않을 때)
# runs <- readRDS("work/robustness_runs.rds")
