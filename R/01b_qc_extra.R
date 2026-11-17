# 01b_qc_extra.R — 추가 QC: 리보솜 유전자 비율 + doublet 전용 도구(scDblFinder)
# 입력: 원본 행렬 (01_qc.R와 같은 2,700개 세포에서 시작), work/pbmc_qc.rds
# 출력: 그림 2b–2d, 표 2b, work/doublet_scores.csv
#       (Part C에서 제거하기로 정하면 work/pbmc_qc.rds를 갱신)
#
# 왜 하나:
#   01_qc.R의 유전자 수 상한(2,500)은 "아주 큰" doublet만 잡는다.
#   작은 세포 두 개가 합쳐진 doublet은 유전자 수가 보통 세포와 비슷해서 기준선에 안 걸린다.
#   scDblFinder는 데이터 안에서 세포 두 개를 인공적으로 섞어 "가짜 doublet"을 만들고,
#   진짜 세포 중 가짜 doublet과 닮은 세포를 찾는다 → 유전자 수가 아니라 "발현 패턴"으로 판단
#
# 실행 방법: Part A → 그림 보기 → Part B → 그림 보기 → Part C의 "직접 정할 값" 채우고 실행

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(scDblFinder)
  library(SingleCellExperiment)
})

counts <- Read10X("data/filtered_gene_bc_matrices/hg19/")
pbmc   <- CreateSeuratObject(counts, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]]   <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 01_qc.R에서 정한 기준 (같은 값으로 맞춘다)
max_feat <- 2500
max_mt   <- 5

# ═══ Part A. 리보솜 유전자 비율 ═══════════════════════════
# 리보솜 단백질 유전자(RPS*, RPL*)는 단백질 공장의 부품 → 거의 모든 세포에서 많이 켜진다.
#   - 비율이 매우 낮은 세포: RNA 품질이 나쁘거나, 리보솜이 적은 특수한 세포(예: 혈소판)
#   - 세포 종류마다 원래 다르다: 쉬고 있는 T세포는 높고, 단핵구는 낮은 편
# → 보통 거르는 기준으로 쓰지 않고, "이상한 세포 집단이 있는지" 확인하는 용도
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
cat("RPS/RPL 유전자 수:", sum(grepl("^RP[SL]", rownames(pbmc))), "\n")

qc <- pbmc@meta.data |> rownames_to_column("barcode")
cat("리보솜 비율 요약:\n")
print(summary(qc$percent.ribo))

r1 <- ggplot(qc, aes(x = "", y = percent.ribo)) +
  geom_violin(fill = "grey90", colour = "grey50") +
  geom_jitter(width = 0.25, size = 0.3, alpha = 0.4) +
  labs(x = NULL, y = "리보솜 유전자 비율 (%)") + theme_bw()
r2 <- ggplot(qc, aes(nFeature_RNA, percent.ribo, colour = percent.mt >= max_mt)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#0072B2"),
                      labels = c("MT 정상", paste0("MT ", max_mt, "% 이상")), name = NULL) +
  labs(x = "검출 유전자 수", y = "리보솜 유전자 비율 (%)") + theme_bw()
r3 <- ggplot(qc, aes(percent.mt, percent.ribo)) +
  geom_point(size = 0.5, alpha = 0.6, colour = "grey40") +
  geom_vline(xintercept = max_mt, linetype = "dashed", colour = "red") +
  labs(x = "미토콘드리아 비율 (%)", y = "리보솜 유전자 비율 (%)") + theme_bw()
ggsave("results/figures/fig2b_ribo.png", r1 + r2 + r3 + plot_layout(widths = c(1, 2, 2)),
       width = 12, height = 4, dpi = 150)

# ═══ Part B. scDblFinder ═══════════════════════════════════
# 1. Seurat 객체 → SingleCellExperiment(Bioconductor 형식)로 바꿔서 넣는다
# 2. 가짜 doublet을 만들 때 무작위로 세포를 고르므로 set.seed로 결과를 고정한다 (재현성)
# 3. 결과: 세포마다 doublet 점수(0~1)와 판정(singlet / doublet)
#    - dbr: 예상 doublet 비율. 10x는 세포 1,000개당 약 0.8% → 2,700개면 약 2%
#      지정하지 않으면 scDblFinder가 세포 수로 자동 계산한다
sce <- SingleCellExperiment(assays = list(counts = LayerData(pbmc, layer = "counts")))
set.seed(42)
sce <- scDblFinder(sce)

qc <- qc |>
  mutate(dbl_score = sce$scDblFinder.score,
         dbl_class = sce$scDblFinder.class,
         over_feat = nFeature_RNA >= max_feat)

cat("\nscDblFinder 판정:\n")
print(table(scDblFinder = qc$dbl_class, `유전자 2,500 이상` = qc$over_feat))

# doublet 판정 세포는 보통 세포와 무엇이 다른가?
dbl_summary <- qc |>
  group_by(dbl_class) |>
  summarise(세포수 = n(),
            UMI_중앙값 = median(nCount_RNA),
            유전자_중앙값 = median(nFeature_RNA),
            MT_중앙값 = round(median(percent.mt), 2),
            MT_5이상 = sum(percent.mt >= max_mt))
print(dbl_summary)
write_csv(dbl_summary, "results/tables/table2b_doublet_summary.csv")
write_csv(select(qc, barcode, nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, dbl_score, dbl_class),
          "work/doublet_scores.csv")

# 그림 2c: 01_qc.R 그림 2 왼쪽에 scDblFinder 판정을 겹쳐 그린다
d1 <- ggplot(arrange(qc, dbl_class == "doublet"), aes(nCount_RNA, nFeature_RNA, colour = dbl_class)) +
  geom_point(size = 0.6, alpha = 0.7) +
  geom_hline(yintercept = max_feat, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c(singlet = "grey70", doublet = "#D55E00"), name = "scDblFinder") +
  labs(x = "UMI 수", y = "검출 유전자 수") + theme_bw()
d2 <- ggplot(qc, aes(dbl_score)) +
  geom_histogram(bins = 50, fill = "grey60") +
  labs(x = "doublet 점수 (0 = 단일 세포, 1 = doublet)", y = "세포 수") + theme_bw()
ggsave("results/figures/fig2c_doublets.png", d1 + d2 + plot_layout(widths = c(3, 2)),
       width = 11, height = 4.5, dpi = 150)

# 그림 2d: doublet 판정 세포는 어떤 마커를 함께 갖고 있나?
#   서로 다른 세포 종류의 마커가 한 세포에서 같이 켜져 있으면 doublet일 가능성이 높다
#   (예: T세포 마커 CD3E + 단핵구 마커 LYZ)
markers <- c(T세포 = "CD3E", B세포 = "MS4A1", 단핵구 = "LYZ", NK세포 = "NKG7", 혈소판 = "PPBP")
mk <- as.matrix(counts[markers, qc$barcode] > 0)
rownames(mk) <- names(markers)
mk_long <- as_tibble(t(mk)) |>
  mutate(dbl_class = qc$dbl_class) |>
  pivot_longer(-dbl_class, names_to = "마커", values_to = "검출") |>
  group_by(dbl_class, 마커) |>
  summarise(비율 = mean(검출), .groups = "drop")
n_types <- colSums(mk[c("T세포", "B세포", "단핵구"), ])   # T/B/단핵구 마커 중 몇 개가 동시에 켜졌나
cat("\nT·B·단핵구 마커가 2개 이상 동시에 검출된 세포 비율\n")
print(tapply(n_types >= 2, qc$dbl_class, mean) |> round(3))

p_mk <- ggplot(mk_long, aes(마커, 비율, fill = dbl_class)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(singlet = "grey60", doublet = "#D55E00"), name = "scDblFinder") +
  labs(x = NULL, y = "마커가 검출된 세포 비율") + theme_bw()
ggsave("results/figures/fig2d_doublet_markers.png", p_mk, width = 7, height = 4, dpi = 150)

# ═══ Part C. 직접 정할 값 ══════════════════════════════════
# 그림 2c, 2d와 위 표를 보고 정한다. 이유는 notes.md "추가 QC"에 적는다.
#   TRUE  : scDblFinder가 doublet으로 판정한 세포를 QC 결과에서 추가로 뺀다
#   FALSE : 빼지 않는다 (01_qc.R 결과 그대로 사용)
# 주의: TRUE로 저장한 뒤 FALSE로 바꾸려면 01_qc.R의 마지막 세 줄(114–116줄)을 먼저 다시 실행해
#       원래 2,638개 파일로 되돌린 뒤 이 Part C를 실행한다. 01_qc.R을 다시 돌렸을 때도 이 Part C를 다시 실행한다.

# ▶ 직접 정할 값
remove_dbl <- TRUE

stopifnot("remove_dbl을 먼저 정하세요 (TRUE 또는 FALSE)" = !is.na(remove_dbl))
pbmc_qc <- readRDS("work/pbmc_qc.rds")
dbl_bc  <- qc$barcode[qc$dbl_class == "doublet"]
pbmc_qc$dbl_class <- qc$dbl_class[match(colnames(pbmc_qc), qc$barcode)]   # 나중에 확인용으로 판정을 붙여 둔다
if (remove_dbl) {
  before  <- ncol(pbmc_qc)
  pbmc_qc <- subset(pbmc_qc, cells = setdiff(colnames(pbmc_qc), dbl_bc))
  cat("\ndoublet 추가 제거:", before, "→", ncol(pbmc_qc), "세포\n")
} else {
  cat("\ndoublet 판정은 표시만 하고 제거하지 않음:", ncol(pbmc_qc), "세포\n")
}
saveRDS(pbmc_qc, "work/pbmc_qc.rds")
