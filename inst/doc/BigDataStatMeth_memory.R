## ----setup, include = FALSE---------------------------------------------------
library(knitr)
library(rmarkdown)
library(BiocStyle)

options(tinytex.verbose = TRUE)

knitr::opts_chunk$set(collapse = TRUE, comment = "", cache=FALSE, message = FALSE, width = 180, crop = NULL)

## ----install_required, eval=FALSE---------------------------------------------
#  # Install BiocManager (if not previously installed)
#  install.packages("BiocManager")
#  
#  # Install required packages
#  BiocManager::install(c("Matrix", "RcppEigen", "RSpectra", "DelayedArray",
#                         "HDF5Array", "rhdf5"))

## ---- install, eval=FALSE-----------------------------------------------------
#  # Install devtools and load library (if not previously installed)
#  install.packages("devtools")
#  library(devtools)
#  
#  # Install BigDataStatMeth
#  install_github("isglobal-brge/BigDataStatMeth")

## ---- load--------------------------------------------------------------------
library(Matrix)
library(DelayedArray)
library(BigDataStatMeth)
library(ggplot2)

## ---- load2-------------------------------------------------------------------
library(microbenchmark)

## ----mat_sim------------------------------------------------------------------

# Define small matrix A and B
set.seed(123456)
n <- 500
p <- 750
A <- matrix(rnorm(n*n), nrow=n, ncol=n)
B <- matrix(rnorm(n*p), nrow=n, ncol=p)


# Define big matrix Abig and Bbig
n <- 1000
p <- 10000
Abig <- matrix(rnorm(n*n), nrow=n, ncol=n)
Bbig <- matrix(rnorm(n*p), nrow=n, ncol=p)

## ----blockmult----------------------------------------------------------------
# Use 10x10 blocks
AxB <- bdblockmult(A, B, block_size = 10)

## ----check_equal_x------------------------------------------------------------
all.equal(AxB, A%*%B)

## ----blockmultparal-----------------------------------------------------------
AxB <- bdblockmult(A, B, block_size = 10, paral = TRUE)

all.equal(AxB,A%*%B)

## ----blockmultbm1-------------------------------------------------------------
# We want to force it to run in memory
AxBNOBig <- bdblockmult(Abig, Bbig, onmemory = TRUE) 

# Run matrix multiplication with data on memory using submatrices of 256x256
AxBBig3000 <- bdblockmult(Abig, Bbig, block_size = 256 , onmemory = TRUE)

## ----bench2, cache=FALSE------------------------------------------------------
bench1 <- microbenchmark( 
  # Parallel block size = 128
  Paral128Mem = bdblockmult(Abig, Bbig, paral = TRUE), 
  # On disk parallel block size = 256
  Paral256Disk = bdblockmult(Abig, Bbig, block_size=256, paral=TRUE), 
  Paral256Mem = bdblockmult(Abig, Bbig, block_size=256, 
                            paral=TRUE, onmemory=TRUE),
  Paral1024Mem = bdblockmult(Abig, Bbig, block_size=1024, 
                             paral=TRUE, onmemory=TRUE), times = 3 )

bench1

## ---- plotbench2, cache=FALSE-------------------------------------------------
ggplot2::autoplot(bench1)

## ----sparsematmult------------------------------------------------------------
k <- 1e3

# Generate 2 sparse matrix x_sparse and y_sparse
set.seed(1)
x_sparse <- sparseMatrix(
   i = sample(x = k, size = k),
   j = sample(x = k, size = k),
   x = rnorm(n = k)
)

set.seed(2)
y_sparse <- sparseMatrix(
   i = sample(x = k, size = k),
   j = sample(x = k, size = k),
   x = rnorm(n = k)
)

d <- bdblockmult_sparse(x_sparse, y_sparse)
f <- x_sparse%*%y_sparse

all.equal(d,f)


## ----bench3, cache=FALSE------------------------------------------------------
res <- microbenchmark( 
  sparse_mult = bdblockmult_sparse(x_sparse, y_sparse),
  matrix_mult = bdblockmult(as.matrix(x_sparse), as.matrix(y_sparse)), 
  RSparse_mult = x_sparse%*% y_sparse, 
  times = 3 )

res

## ---- plotbench3, cache=FALSE-------------------------------------------------
ggplot2::autoplot(res)

## ----crossprod----------------------------------------------------------------
n <- 500
m <- 250
A <- matrix(rnorm(n*m), nrow=n, ncol=m)

# Cross Product of a standard R matrix
cpA <- bdCrossprod(A)

## ----check_cp-----------------------------------------------------------------
all.equal(cpA, crossprod(A))

## ----nocrossprod--------------------------------------------------------------
# Transposed Cross Product R matrices
tcpA <- bdtCrossprod(A)

## ----check_tcp----------------------------------------------------------------
all.equal(tcpA, tcrossprod(A))

## ----benchmark_bdcrossprod----------------------------------------------------
res <- microbenchmark(
  bdcrossp_tr = bdtCrossprod(A),
  rcrossp_tr = tcrossprod(A),
  bdcrossp = bdCrossprod(A),
  rcrossp = crossprod(A),
  times = 3)
res

## ---- plotbenchmark_bdcrossprod-----------------------------------------------
ggplot2::autoplot(res)

## ----wcrossprod---------------------------------------------------------------
n <- 250
X <- matrix(rnorm(n*n), nrow=n, ncol=n)
u <- runif(n)
w <- u * (1 - u)

wcpX <- bdwproduct(X, w,"xwxt")

wcpX[1:5,1:5]

## ----check_wcp----------------------------------------------------------------
all.equal( wcpX, X%*%diag(w)%*%t(X) )

## ----wtcrossprod--------------------------------------------------------------
wtcpX <- bdwproduct(X, w,"xtwx")

wtcpX[1:5,1:5]

## ----check_wtcp---------------------------------------------------------------
all.equal(wtcpX, t(X)%*%diag(w)%*%X)

## ----invChols-----------------------------------------------------------------
# Generate a positive definite matrix
Posdef <- function (n, ev = runif(n, 0, 10))
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

A <- Posdef(n = 500, ev = 1:500)

invchol <- bdInvCholesky(A)

round(invchol[1:5,1:5],8)

## ----invCholsequal------------------------------------------------------------
all.equal(invchol, solve(A))

## ----benchmark_invChol--------------------------------------------------------
res <- microbenchmark(invchol = bdInvCholesky(A),
                      invcholR = solve(A),
                      times = 3)
res

## ---- plotbenchmark_invChol---------------------------------------------------
ggplot2::autoplot(res)

## ----pseudoinv----------------------------------------------------------------

m <- 1300
n <- 1200

A <- matrix(rnorm(n*m), nrow=n, ncol=m)

pseudoinv <- bdpseudoinv(A)
zapsmall(pseudoinv)[1:5,1:5]

## ----QRdecomp-----------------------------------------------------------------


QR_A <- bdQR(A)
QR_R <- qr(A)

# Show results for Q
zapsmall(QR_A$Q[1:5,1:5])

# Show results for R
zapsmall(QR_A$R[1:5,1:5])

# Test Q 
all.equal(QR_A$Q, qr.Q(QR_R), check.attributes=FALSE)


## ----matrix_matrixEQ----------------------------------------------------------

# Simulate data
m <- 1000
n <- 1000

A <- matrix(runif(n*m), nrow = n, ncol = m)
B <- matrix(runif(n*2), nrow = n)

# Solve matrix equation
X <- bdSolve(A, B)

# Show results
X[1:5,]


## ----check_matrix_matrixEQ----------------------------------------------------

testB <- bdblockmult(A,X)

B[1:5,]
testB[1:5,]

all.equal(B, testB)


## ----svd_default--------------------------------------------------------------
# Matrix simulation
set.seed(413)
n <- 500
A <- matrix(rnorm(n*n), nrow=n, ncol=n)

# SVD
bsvd <- bdSVD(A)

# Singular values, and right and left singular vectors
bsvd$d[1:5]
bsvd$u[1:5,1:5]
bsvd$v[1:5,1:5]

## ----check_svd----------------------------------------------------------------
all.equal( sqrt( svd( tcrossprod( scale(A) ) )$d[1:10] ), bsvd$d[1:10] )

## ----svd_nonorm---------------------------------------------------------------
bsvd <- bdSVD(A, bcenter = FALSE, bscale = FALSE)


bsvd$d[1:5]
bsvd$u[1:5,1:5]
bsvd$v[1:5,1:5]

all.equal( sqrt(svd(tcrossprod(A))$d[1:10]), bsvd$d[1:10] )


## ----mirNA--------------------------------------------------------------------
data(miRNA)
data(cancer)
dim(miRNA)

## ----tab----------------------------------------------------------------------
table(cancer)

## ----pca----------------------------------------------------------------------
pc <- prcomp(miRNA)

## ----plot_pca-----------------------------------------------------------------
plot(pc$x[, 1], pc$x[, 2],
     main = "miRNA data on tumor samples",
     xlab = "PC1", ylab = "PC2", type="n")
abline(h=0, v=0, lty=2)
points(pc$x[, 1], pc$x[, 2], col = cancer,
       pch=16, cex=1.2)
legend("topright", levels(cancer), pch=16, col=1:3)

## ----cia_da-------------------------------------------------------------------
miRNA.c <- sweep(miRNA, 2, colMeans(miRNA), "-")
svd.da <- bdSVD(miRNA.c, bcenter = FALSE, bscale = FALSE)

## ----plot_svd_da--------------------------------------------------------------
plot(svd.da$u[, 1], svd.da$u[, 2],
     main = "miRNA data on tumor samples",
     xlab = "PC1", ylab = "PC2", type="n")
abline(h=0, v=0, lty=2)
points(svd.da$u[, 1], svd.da$u[, 2], col = cancer,
       pch=16, cex=1.2)
legend("topright", levels(cancer), pch=16, col=1:3)

## ----sesinfo------------------------------------------------------------------
sessionInfo()

