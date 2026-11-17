# 00_setup.R — 패키지 설치와 PBMC 3k 데이터 다운로드
# 프로젝트 폴더(pbmc-scrna)를 작업 디렉터리로 두고 실행한다.

options(repos = c(CRAN = "https://cloud.r-project.org"))

# 1. 패키지 -----------------------------------------------------------------
# Seurat: 단일세포 분석 (QC → 정규화 → PCA → 군집화 → UMAP → 마커)
# biomaRt(lncRNA 단계)는 나중에 Bioconductor로 설치한다:
#   install.packages("BiocManager"); BiocManager::install("biomaRt")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

# 2. 데이터 -----------------------------------------------------------------
# 10x Genomics PBMC 3k: 건강한 공여자 1명의 말초혈액 단핵세포 2,700개 (hg19)
url  <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
dest <- "data/pbmc3k_filtered_gene_bc_matrices.tar.gz"
dir.create("data", showWarnings = FALSE)
if (!file.exists(dest)) download.file(url, dest, mode = "wb")
if (!dir.exists("data/filtered_gene_bc_matrices")) untar(dest, exdir = "data")

# 3. 확인 -------------------------------------------------------------------
library(Seurat)
counts <- Read10X("data/filtered_gene_bc_matrices/hg19/")
cat("유전자 x 세포:", dim(counts), "\n")
cat("세포당 검출 유전자 수 중앙값:", median(colSums(counts > 0)), "\n")
cat("세포당 UMI 수 중앙값:", median(colSums(counts)), "\n")
cat("Seurat 버전:", as.character(packageVersion("Seurat")), "\n")
