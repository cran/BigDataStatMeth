## ----setup, include = FALSE---------------------------------------------------
library(knitr)
library(BiocStyle)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "",
  cache = FALSE,
  message = FALSE,
  warning = FALSE,
  width = 80,
  crop = NULL
)

## ----load-package-------------------------------------------------------------
library(BigDataStatMeth)

## ----options-view-------------------------------------------------------------
old_opts <- hdf5matrix_options()
old_opts

## ----options-set--------------------------------------------------------------
hdf5matrix_options(
  paral = TRUE,
  threads = 2L,
  block_size = 512L,
  compression = 6L
)

hdf5matrix_options()

## ----create-matrix------------------------------------------------------------
h5file <- tempfile(fileext = ".h5")

set.seed(123)
X <- matrix(rnorm(500 * 100), nrow = 500, ncol = 100)

X_h5 <- hdf5_create_matrix(
  filename = h5file,
  dataset = "data/X",
  data = X,
  overwrite = TRUE
)

X_h5
dim(X_h5)
nrow(X_h5)
ncol(X_h5)

## ----open-matrix--------------------------------------------------------------
list_datasets(h5file)

X_h5_reopened <- hdf5_matrix(
  filename = h5file,
  path = "data/X"
)

dim(X_h5_reopened)

## ----import-tabular-data------------------------------------------------------
csv_file <- system.file(
  "extdata", "colesterol.csv",
  package = "BigDataStatMeth"
)

if (!nzchar(csv_file)) {
  csv_file <- system.file(
    "data", "colesterol.csv",
    package = "BigDataStatMeth"
  )
}

if (!nzchar(csv_file) && file.exists("colesterol.csv")) {
  csv_file <- "colesterol.csv"
}

stopifnot(nzchar(csv_file))

h5_csv <- tempfile(fileext = ".h5")

cholesterol_h5 <- hdf5_import(
  source = csv_file,
  filename = h5_csv,
  dataset = "cholesterol/data",
  sep = ",",
  header = TRUE,
  overwrite = TRUE
)

cholesterol_h5
dim(cholesterol_h5)
cholesterol_h5[1:5, 1:min(6, ncol(cholesterol_h5))]

## ----subsetting---------------------------------------------------------------
X_h5[1:5, 1:6]

sub_X <- X_h5[1:20, 1:10]
dim(sub_X)

## ----assignment---------------------------------------------------------------
X_edit <- hdf5_create_matrix(
  h5file,
  "data/X_edit",
  data = X[1:10, 1:6],
  overwrite = TRUE
)

X_edit[1, 1] <- 999
X_edit[1:3, 1:3]

## ----dimnames-example---------------------------------------------------------
DN_h5 <- hdf5_create_matrix(
  h5file,
  "data/dimnames_example",
  data = matrix(seq_len(12), nrow = 4, ncol = 3),
  overwrite = TRUE
)

rownames(DN_h5) <- paste0("sample_", seq_len(nrow(DN_h5)))
colnames(DN_h5) <- paste0("feature_", seq_len(ncol(DN_h5)))

rownames(DN_h5)
colnames(DN_h5)
dimnames(DN_h5)

## ----convert-to-memory--------------------------------------------------------
X_small <- as.matrix(X_h5[1:10, 1:5])
X_small

## ----arithmetic---------------------------------------------------------------
set.seed(1)

A <- matrix(rnorm(300 * 80), nrow = 300, ncol = 80)
B <- matrix(rnorm(300 * 80), nrow = 300, ncol = 80)
C <- matrix(rnorm(80 * 40), nrow = 80, ncol = 40)

A_h5 <- hdf5_create_matrix(
  h5file, "data/A", data = A,
  overwrite = TRUE
)

B_h5 <- hdf5_create_matrix(
  h5file, "data/B", data = B,
  overwrite = TRUE
)

C_h5 <- hdf5_create_matrix(
  h5file, "data/C", data = C,
  overwrite = TRUE
)

S_h5 <- A_h5 + B_h5
D_h5 <- A_h5 - B_h5

S_h5
dim(S_h5)

all.equal(as.matrix(S_h5), A + B)
all.equal(as.matrix(D_h5), A - B)

## ----arithmetic-output-location-----------------------------------------------
list_datasets(h5file, group = "OUTPUT", recursive = TRUE)

## ----multiplication-----------------------------------------------------------
M_h5 <- A_h5 %*% C_h5

M_h5
dim(M_h5)
all.equal(as.matrix(M_h5), A %*% C)

## ----crossproducts------------------------------------------------------------
XtX_h5 <- crossprod(
  A_h5,
  outgroup = "RESULTS",
  outdataset = "A_crossprod"
)

XXt_h5 <- tcrossprod(
  A_h5,
  outgroup = "RESULTS",
  outdataset = "A_tcrossprod"
)

XtX_h5
list_datasets(h5file, group = "RESULTS", recursive = TRUE)

all.equal(as.matrix(XtX_h5), crossprod(A))
all.equal(as.matrix(XXt_h5), tcrossprod(A))

## ----bind-example-------------------------------------------------------------
A1_h5 <- hdf5_create_matrix(
  h5file,
  "bind/A1",
  data = A[1:50, 1:10],
  overwrite = TRUE
)

A2_h5 <- hdf5_create_matrix(
  h5file,
  "bind/A2",
  data = A[1:50, 11:20],
  overwrite = TRUE
)

Cbind_h5 <- cbind(
  A1_h5, A2_h5,
  out_group = "BIND_RESULTS",
  out_dataset = "A1_A2_cbind",
  overwrite = TRUE
)

Rbind_h5 <- rbind(
  A1_h5, A1_h5,
  out_group = "BIND_RESULTS",
  out_dataset = "A1_A1_rbind",
  overwrite = TRUE
)

dim(Cbind_h5)
dim(Rbind_h5)

all.equal(as.matrix(Cbind_h5), cbind(A[1:50, 1:10], A[1:50, 11:20]))
all.equal(as.matrix(Rbind_h5), rbind(A[1:50, 1:10], A[1:50, 1:10]))

## ----aggregations-------------------------------------------------------------
all.equal(colMeans(A_h5), colMeans(A))
all.equal(rowSums(A_h5), rowSums(A))

all.equal(colVars(A_h5), apply(A, 2, var))
all.equal(rowSds(A_h5), apply(A, 1, sd))

## ----scaling------------------------------------------------------------------
A_scaled_h5 <- scale(A_h5)
A_scaled <- scale(A)

all.equal(
  as.matrix(A_scaled_h5),
  A_scaled,
  check.attributes = FALSE
)

## ----correlation--------------------------------------------------------------
Cor_h5 <- cor(A_h5)

all.equal(
  as.matrix(Cor_h5),
  cor(A),
  tolerance = 1e-6
)

## ----sweep-example------------------------------------------------------------
col_means_h5 <- hdf5_create_matrix(
  h5file,
  "stats/A_col_means",
  data = matrix(colMeans(A), nrow = 1),
  overwrite = TRUE
)

A_centered_h5 <- sweep(A_h5, MARGIN = 2, STATS = col_means_h5, FUN = "-")

all.equal(
  as.matrix(A_centered_h5),
  sweep(A, MARGIN = 2, STATS = colMeans(A), FUN = "-"),
  check.attributes = FALSE
)

## ----svd-example--------------------------------------------------------------
set.seed(123)
X_svd <- matrix(rnorm(120 * 300), nrow = 120, ncol = 300)

X_svd_h5 <- hdf5_create_matrix(
  h5file,
  "decomp/X_svd",
  data = X_svd,
  overwrite = TRUE
)

svd_h5 <- svd(
  X_svd_h5,
  nu = 10,
  nv = 10,
  center = TRUE,
  scale = TRUE,
  overwrite = TRUE
)

head(svd_h5$d)
dim(svd_h5$u)
dim(svd_h5$v)

## ----svd-validation-----------------------------------------------------------
svd_r <- svd(scale(X_svd), nu = 10, nv = 10)

all.equal(
  svd_h5$d[1:10],
  svd_r$d[1:10],
  tolerance = 1e-6
)

## ----block-svd-example--------------------------------------------------------
svd_blk_h5 <- svd(
  X_svd_h5,
  nu = 5,
  nv = 5,
  k = 4,
  q = 1,
  threads = 2,
  overwrite = TRUE
)

head(svd_blk_h5$d)
dim(svd_blk_h5$u)
dim(svd_blk_h5$v)

## ----pca-example--------------------------------------------------------------
set.seed(124)

n <- 180
p <- 40
group <- rep(c("Group 1", "Group 2", "Group 3"), each = n / 3)

X_pca <- matrix(rnorm(n * p), nrow = n, ncol = p)
X_pca[group == "Group 2", 1:8] <- X_pca[group == "Group 2", 1:8] + 1.5
X_pca[group == "Group 3", 9:16] <- X_pca[group == "Group 3", 9:16] - 1.5

X_pca_h5 <- hdf5_create_matrix(
  h5file,
  "decomp/X_pca",
  data = X_pca,
  overwrite = TRUE
)

pca_h5 <- prcomp(
  X_pca_h5,
  center = TRUE,
  scale. = TRUE,
  ncomponents = 5,
  overwrite = TRUE
)

pca_h5
head(pca_h5$sdev)

## ----pca-plot, fig.width = 7, fig.height = 5, fig.cap = "PCA scores computed from an HDF5-backed matrix."----
class(pca_h5$x)
dim(pca_h5$x)

pca_scores <- as.matrix(pca_h5$x[, 1:2])
pca_df <- data.frame(
  PC1 = pca_scores[, 1],
  PC2 = pca_scores[, 2],
  group = group
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, colour = group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.85) +
    ggplot2::stat_ellipse(linewidth = 0.6, alpha = 0.7) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "PCA of an HDF5-backed matrix",
      subtitle = "Scores returned by prcomp.HDF5Matrix()",
      x = "PC1",
      y = "PC2",
      colour = "Synthetic group"
    )
} else {
  plot(
    pca_df$PC1,
    pca_df$PC2,
    pch = 19,
    xlab = "PC1",
    ylab = "PC2",
    main = "PCA of an HDF5-backed matrix"
  )
}

## ----qr-example---------------------------------------------------------------
QR_h5 <- qr(A_h5, thin = TRUE, overwrite = TRUE)

Q <- as.matrix(QR_h5$Q)
R <- as.matrix(QR_h5$R)

all.equal(Q %*% R, A, tolerance = 1e-6)
all.equal(crossprod(Q), diag(ncol(Q)), tolerance = 1e-6)

## ----tsqr-example-------------------------------------------------------------
X_tsqr <- matrix(rnorm(600 * 30), nrow = 600, ncol = 30)
X_tsqr_h5 <- hdf5_create_matrix(
  h5file,
  "decomp/X_tsqr",
  data = X_tsqr,
  overwrite = TRUE
)

QR_tsqr_h5 <- qr(
  X_tsqr_h5,
  thin = TRUE,
  method = "tsqr",
  threads = 2L,
  overwrite = TRUE
)

dim(QR_tsqr_h5$Q)
dim(QR_tsqr_h5$R)

## ----chol-solve-example-------------------------------------------------------
set.seed(321)
Z <- matrix(rnorm(200 * 50), nrow = 200, ncol = 50)
SPD <- crossprod(Z) + diag(1e-6, 50)

SPD_h5 <- hdf5_create_matrix(
  h5file,
  "decomp/SPD",
  data = SPD,
  overwrite = TRUE
)

Ch_h5 <- chol(SPD_h5, overwrite = TRUE)
Inv_h5 <- solve(SPD_h5, overwrite = TRUE)

Ch <- as.matrix(Ch_h5)
chol_error <- min(
  max(abs(crossprod(Ch) - SPD)),
  max(abs(tcrossprod(Ch) - SPD))
)
chol_error < 1e-6

all.equal(as.matrix(Inv_h5), solve(SPD), tolerance = 1e-6)

## ----eigen-pseudoinverse------------------------------------------------------
Eig_h5 <- eigen(
  SPD_h5,
  symmetric = TRUE,
  k = 5L,
  overwrite = TRUE
)

head(Eig_h5$values)
dim(Eig_h5$vectors)

Pinv_h5 <- pseudoinverse(
  A_h5,
  overwrite = TRUE
)

dim(Pinv_h5)

## ----compression-example------------------------------------------------------
set.seed(123)
X_cmp <- round(matrix(rnorm(2500 * 250), nrow = 2500, ncol = 250), 2)

f0 <- tempfile(fileext = ".h5")
f6 <- tempfile(fileext = ".h5")

t0 <- system.time(
  hdf5_create_matrix(
    f0, "data/X",
    data = X_cmp,
    compression = 0,
    overwrite = TRUE
  )
)

t6 <- system.time(
  hdf5_create_matrix(
    f6, "data/X",
    data = X_cmp,
    compression = 6,
    overwrite = TRUE
  )
)

data.frame(
  compression = c(0, 6),
  elapsed = c(t0[["elapsed"]], t6[["elapsed"]]),
  file_size_MB = round(file.info(c(f0, f6))$size / 1024^2, 3)
)

## ----close-single-object------------------------------------------------------
close(X_h5_reopened)

## ----restore-options-and-close-all--------------------------------------------
hdf5matrix_options(
  paral = old_opts$paral,
  block_size = old_opts$block_size,
  threads = old_opts$threads,
  compression = old_opts$compression
)

hdf5_close_all()
gc()

## ----session-info-------------------------------------------------------------
sessionInfo()

