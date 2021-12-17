## ----setup, include = FALSE---------------------------------------------------
library(knitr)
library(BiocStyle)

knitr::opts_chunk$set(collapse = TRUE, comment = "", cache = TRUE, message = FALSE, width = 180)

## ---- cleanup, echo=FALSE, include=FALSE--------------------------------------
if( isTRUE(file.exists('delayed.hdf5'))) {
    file.remove('delayed.hdf5')
}
if( isTRUE(file.exists('robject.hdf5'))){
    file.remove('robject.hdf5')
}
if( isTRUE(file.exists('colesterol_file.hdf5'))){
    file.remove('colesterol_file.hdf5')
}
if( isTRUE(file.exists('tmp_blockmult.hdf5'))){
    file.remove('tmp_blockmult.hdf5')
}

## ----install_required, eval=FALSE---------------------------------------------
#  # Install BiocManager (if not previously installed)
#  install.packages("BiocManager")
#  
#  # Install required packages
#  BiocManager::install(c("Matrix", "RcppEigen", "RSpectra",
#                         "beachmat", "DelayedArray",
#                         "HDF5Array", "rhdf5"))

## ---- install, eval=FALSE-----------------------------------------------------
#  # Install devtools and load library (if not previously installed)
#  install.packages("devtools")
#  library(devtools)
#  
#  # Install BigDataStatMeth
#  install_github("isglobal-brge/BigDataStatMeth")

## ---- load, cache=FALSE-------------------------------------------------------
library(rhdf5)
library(BigDataStatMeth)

## ----hdf5Img, out.width = '100%', fig.align = 'center', fig.cap = "HDF5 hierarchical structure", echo=FALSE----
knitr::include_graphics("imgs/hdf5_squema.jpg")

## ----hdf5Create---------------------------------------------------------------
library(rhdf5)

set.seed(5234)
n <- 500
m <- 600
A <- matrix(rnorm(n*m,mean=0,sd=1), n,m)

# We also can create a dataset from R matrix object
bdCreate_hdf5_matrix_file(filename = "robject.hdf5", 
                          object = A, 
                          group = "INPUT", 
                          dataset = "A")

## ----ls-----------------------------------------------------------------------
dir()

## ----hdf5AddDataset-----------------------------------------------------------
set.seed(5234)
n <- 100
m <- 10000
A <- matrix(rnorm(n*m,mean=0,sd=1), n, m)

set.seed(5234)
n <- 50
m <- 12000
B <- matrix(rnorm(n*m,mean=3,sd=0.5), n, m)

# We create another data file (delayed.hdf5) with a matrix A.
# The group is called INPUT
bdCreate_hdf5_matrix_file(filename = "delayed.hdf5", 
                        object = A, 
                        group = "INPUT", 
                        dataset = "A")

# And them, we add another matrix B to the same group
bdAdd_hdf5_matrix(object = B, 
                filename = "delayed.hdf5", 
                group = "INPUT", 
                dataset = "B")

## ----hdf5AddtransposedDataset-------------------------------------------------
# Create dataset BTransposed with B transposed matrix in delayed.hdf5 file at INPUT group
bdAdd_hdf5_matrix(B, "delayed.hdf5", "INPUT", "BTransposed", transp = TRUE);

## ----hdf5Show-----------------------------------------------------------------
# Examine hierarchy before open file
h5ls("delayed.hdf5")

## ----hdf5Open, cache=FALSE----------------------------------------------------
# Open file
h5fdelay <- H5Fopen("delayed.hdf5")
# Show hdf5 hierarchy (groups)
h5fdelay

## ----hdf5Dataset--------------------------------------------------------------
Bdata <- h5fdelay$INPUT$B
Bdata[1:3,1:5]

## ----hdf5DatasetTransposed----------------------------------------------------
BdataTrans <- h5fdelay$INPUT$BTransposed
BdataTrans[1:5,1:3]

## ----get_pvals----------------------------------------------------------------
out <- sapply(1:100, function(i){
  mean(h5fdelay$INPUT$A[,i])
  })
head(out)

## ----hdf5Close----------------------------------------------------------------
# Close delayed.hdf5 file
H5Fclose(h5fdelay)

# Open 2 files and close all
h5fdelay <- H5Fopen("delayed.hdf5")
h5fr <- H5Fopen("robject.hdf5")

h5closeAll()

## ----convert_HDF5, cache=FALSE------------------------------------------------
bdImport_text_to_hdf5(filename = "colesterol.csv", 
                      sep=',', 
                      outputfile = "colesterol_file.hdf5", 
                      outGroup = "COLESTEROL", 
                      outDataset = "COLESTEROLDATA", 
                      header = TRUE)


## ----read_data_col_HDF5-------------------------------------------------------

h5ls("colesterol_file.hdf5")

# We can open the file and have access to the data
h5Colesterol <- H5Fopen("colesterol_file.hdf5")

# Show hdf5 content dataset
h5Colesterol$COLESTEROL$COLESTEROLDATA[1:5, 1:6]

# Show colnames 
head(h5Colesterol$COLESTEROL$.COLESTEROLDATA_dimnames$`2`)

H5Fclose(h5Colesterol)

## ----blockmultbm1-------------------------------------------------------------
n <- 1000
p <- 10000
Abig <- matrix(rnorm(n*n), nrow=n, ncol=n)
Bbig <- matrix(rnorm(n*p), nrow=n, ncol=p)


# We want to force it to run in memory
AxBNOBig <- bdblockmult(Abig, Bbig, bigmatrix = 100000)

# We consider a big matrix if number of rows or columns are > 500
AxBBig3000 <- bdblockmult(Abig, Bbig, bigmatrix = 500)

## ----blockmultresmat----------------------------------------------------------
class(AxBNOBig)
AxBNOBig[1:5,1:5]

## ----blockmultresfile---------------------------------------------------------
h5data <- H5Fopen(AxBBig3000$file)

# We can get where data is stored
AxBBig3000$dataset

# We observe that the data is in folder OUTPUT dataset C
reshdf5 <- h5data$OUTPUT$C

reshdf5[1:5,1:5]

all.equal(reshdf5, AxBNOBig)

## ----blockmultresfileclose----------------------------------------------------
# Close file
H5Fclose(h5data)

## ----blockmult_hdf5_exec------------------------------------------------------
res <- bdblockmult_hdf5(filename = 'tmp_blockmult.hdf5', group="INPUT",
                        a="A", b="B", outgroup = 'HDF5_RES')

# We show the hdf5 content
 h5ls(res$file)

## ----blockmult_hdf5_res-------------------------------------------------------
# Get content
h5data <- H5Fopen(res$file)

# We can get where data is stored
res$dataset

# We get the results with bdblockmult (previous step)
resblockmult <- h5data$OUTPUT$C

# and the results obtained from bdblockmult_hdf5
resblockmult_hdf5 <- h5data$HDF5_RES$A_x_B
resblockmult_hdf5[1:5,1:5]

# Close the file
h5closeAll()

all.equal(resblockmult, resblockmult_hdf5)

## ----crossprod----------------------------------------------------------------
res <- bdCrossprod_hdf5(filename = 'robject.hdf5', 
                        group="INPUT", 
                        A="A", 
                        outgroup = "RESULTS")

## ----check_cp-----------------------------------------------------------------

# Get content
h5data <- H5Fopen(res$file)

# We can get wher data is stored
res$dataset

# We get the Crossprod Results and the original matrix
resCProd <- h5data$RESULTS$CrossProd_AxA
A <- h5data$INPUT$A

# Close the file
h5closeAll()

# Show results
resCProd[1:5,1:5]

all.equal(resCProd, crossprod(A))

## ----crossprod_ab-------------------------------------------------------------

set.seed(5234)
n <- 500
m <- 600
B <- matrix(rnorm(n*m,mean=0,sd=1), n,m)

bdAdd_hdf5_matrix(B, filename = 'robject.hdf5', group="INPUT", dataset = "B2")

# Get Crossprod with two matrix
res <- bdCrossprod_hdf5(filename = 'robject.hdf5', group="INPUT", A="A", 
                        groupB = "INPUT", B = "B2", outgroup = "RESULTS")


## ----check_cp_ab--------------------------------------------------------------

# Get content
h5data <- H5Fopen(res$file)

# We can get wher data is stored
res$dataset

# We get the Crossprod Results and the original matrix
resCProd <- h5data$RESULTS$CrossProd_AxB2
A <- h5data$INPUT$A
B <- h5data$INPUT$B2

# Close the file
h5closeAll()

# Show results
resCProd[1:5,1:5]

all.equal(resCProd, crossprod(A,B))

## ----tcrossprod---------------------------------------------------------------
res <- bdtCrossprod_hdf5(filename = 'robject.hdf5', group="INPUT", A="A", outgroup = "RESULTS")

## ----check_tcp----------------------------------------------------------------

# Get content
h5data <- H5Fopen(res$file)

# We can get wher data is stored
res$dataset

# We get the Crossprod Results and the original matrix
restCProd <- h5data$RESULTS$tCrossProd_AxA
A <- h5data$INPUT$A

# Close the file
h5closeAll()

# Show results
restCProd[1:5,1:5]

all.equal(restCProd, tcrossprod(A))

## ----tcrossprod_ab------------------------------------------------------------

set.seed(5234)
n <- 500
m <- 600
B <- matrix(rnorm(n*m,mean=0,sd=1), n,m)

bdAdd_hdf5_matrix(B, filename = 'robject.hdf5', group="INPUT", dataset = "B3")

# Get Crossprod with two matrix
res <- bdtCrossprod_hdf5(filename = 'robject.hdf5', group="INPUT", A="A", 
                         groupB = "INPUT", B = "B3", outgroup = "RESULTS")


## ----check_tcp_ab-------------------------------------------------------------

# Get content
h5data <- H5Fopen(res$file)

# We can get wher data is stored
res$dataset

# We get the Crossprod Results and the original matrix
restCProd <- h5data$RESULTS$tCrossProd_AxB3
A <- h5data$INPUT$A
B <- h5data$INPUT$B3

# Close the file
h5closeAll()

# Show results
restCProd[1:5,1:5]

all.equal(restCProd, tcrossprod(A,B))

## ----BSVDImg, out.width = '100%', fig.align = 'center', fig.cap = "Flowchart for a two-level hierarchical Block SVD algorithm", echo=FALSE----
knitr::include_graphics("imgs/blocksvd.png")

## ----BlockSVDNorm-------------------------------------------------------------
# Create dataframe data with 'odata' matrix in delayed hdf5 file at OMIC group
set.seed(5234)
n <- 150000
m <- 50
odata <- matrix(rnorm(n*m, mean=0, sd=1), n,m)

bdAdd_hdf5_matrix(odata, "delayed.hdf5", "OMICS", "data")

# Direct from hdf5 data file
svdh5 <- bdSVD_hdf5("delayed.hdf5", "OMICS", "data")

# with R implementation from data in memory
test <- H5Fopen("delayed.hdf5")
# get results svd (d)
svd_hdf5_d <- test$SVD$data$d[1:7]
# get data
omdata <- test$OMICS$data
h5closeAll()

# Results in hdf5 file for d
svd_hdf5_d[1:7]

svd <- svd(scale(omdata))
svd$d[1:7]

## ----BlockSVDNotNorm----------------------------------------------------------
# Direct from hdf5 data file
svdh5 <- bdSVD_hdf5("delayed.hdf5", "OMICS", "data", 
                    bcenter = FALSE, bscale = FALSE)

# get results svd (d)
test <- H5Fopen("delayed.hdf5")
svd_hdf5_d <- test$SVD$data$d[1:7] 
h5closeAll()

# SVD (d) from file - data not normalized
svd_hdf5_d

# with R implementation from data in memory
svd <- svd(omdata)
svd$d[1:7]

## ----BlockSVDk4---------------------------------------------------------------
# Block decomposition with 1 level and 4 local SVDs at each level
svdh5 <- bdSVD_hdf5("delayed.hdf5", "OMICS", "data", q=1, k=4)

# get results svd (d)
fprova <- H5Fopen("delayed.hdf5")
  svd_hdf5_d <- fprova$SVD$data$d[1:7]
h5closeAll()

# SVD (d) from file - data not normalized
svd_hdf5_d

# with R implementation from data in memory
svd <- svd(scale(omdata))
svd$d[1:7]

## ----snpgenotype--------------------------------------------------------------
set.seed(108432)
geno.sim <- matrix(sample(0:3, 10000, replace = TRUE), byrow = TRUE, ncol = 10)
bdAdd_hdf5_matrix(geno.sim, "delayed.hdf5", "OMICS", "geno")

# Get data and show the first 5 rows
h5fsvd = H5Fopen("delayed.hdf5")
geno <- h5fsvd$OMICS$geno
h5closeAll()

geno[1:5,1:10]


## ----impute-------------------------------------------------------------------
bdImpute_snps_hdf5("delayed.hdf5", group="OMICS", dataset="geno",
                outgroup="OMICS", outdataset="imputedgeno")

## ----impute2------------------------------------------------------------------
# Get imputed data and show the first 5 rows
h5fsvd = H5Fopen("delayed.hdf5")
imputedgeno <- h5fsvd$OMICS$imputedgeno
h5closeAll()

imputedgeno[1:5,1:10]

## ----normalize----------------------------------------------------------------

bdNormalize_hdf5("delayed.hdf5", group="OMICS", dataset="imputedgeno")

## ----normalizeres-------------------------------------------------------------
# Geno after normalization
h5fsvd = H5Fopen("delayed.hdf5")
genonormalized <- h5fsvd$NORMALIZED$OMICS$geno
h5closeAll()

genonormalized[1:5,1:10]

## ----removelow----------------------------------------------------------------
bdRemovelowdata("delayed.hdf5", group="OMICS", dataset="geno",
                outgroup="OMICS", outdataset="removedgeno",
                bycols = TRUE, pcent = 0.39)

## ----removelow2---------------------------------------------------------------
# Get imputed data and show the first 5 rows
h5fsvd = H5Fopen("delayed.hdf5")
removedgeno <- h5fsvd$OMICS$removedgeno
h5closeAll()

removedgeno[1:5,1:10]

## ----sesinfo------------------------------------------------------------------
sessionInfo()

