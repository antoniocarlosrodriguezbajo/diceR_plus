#' Gene expression data for High Grade Serous Carcinoma from TCGA
#'
#' There are 489 samples measured on 321 genes. Sample IDs are in the row names
#' and gene names are in the column names. This data set is used for clustering
#' HGSC into subtypes with prognostic significance. The cluster assignments
#' obtained by TCGA are indicated by the last six characters of each row name in
#' `hgsc`: `MES.C1`, `IMM.C2`, `DIF.C4`, and `PRO.C5`
#'
#' @format A data frame with 489 rows and 321 columns.
"hgsc"

#' Near-infrared spectra of homogenized meat samples
#'
#' This dataset contains spectral measurements (400nm-2049nm) of 231 meat samples
#' from five species: Chicken, Beef, Turkey, Lamb, and Pork. The first 1050
#' columns record spectral data, while the last column indicates the sample type.
#'
#' @references
#' McElhinney, J., Downey, G., Fearn, T. (1999). Chemometric processing of visible
#' and near infrared reflectance spectra for species identification in selected raw homogenised meats.
#' *Journal of Near Infrared Spectroscopy, 7*, 145–154.
#'
#' @format
#' A data frame with 231 rows and 1051 columns.
"meats"

#' Gene Expression Data for Acute Leukemia
#'
#' This dataset contains gene expression profiles from patients with
#' Acute Lymphoblastic Leukemia (ALL) and Acute Myeloid Leukemia (AML).
#' The data includes expression levels for 7129 genes across 72 patient samples.
#'
#' Each row represents a patient, while each column (except the last) corresponds
#' to a gene. The final column `Diagnosis` classifies the patients as either `ALL` or `AML`.
#'
#' @format A data frame with 72 rows and 7130 columns:
#' \describe{
#'   \item{X1, X2, ..., X7129}{Gene expression values (numeric).}
#'   \item{class}{Factor indicating the leukemia type (`ALL` or `AML`).}
#' }
#'
#' @references
#' Fodor, S.P.A. (1997). Massively parallel genomics. *Science*, **277**(5324), 393.
#'
"ALLAML"

#'
#' Gene Expression Data for Prostate Cancer
#'
#' This dataset contains gene expression profiles from patients with
#' prostate cancer. The data includes expression levels for 5966 genes across 102 patient samples.
#'
#' Each row represents a patient, while each column (except the last) corresponds
#' to a gene. The final column `Diagnosis` classifies the patients as either `Tumor` or `Normal`.
#'
#' @format A data frame with 102 rows and 5967 columns:
#' \describe{
#'   \item{X1, X2, ..., X5966}{Gene expression values (numeric).}
#'   \item{class}{Factor indicating the sample type (`Tumor` or `Normal`).}
#' }
#'
#' @references
#' Singh, D.; Febbo, P.G.; Ross, K.; Jackson, D.G.; Manola, J.; Ladd, C.; Tamayo, P.; Renshaw, A.A.; DAmico, A.V.; Richie, J.P.; et al. (2002). Gene expression correlates of clinical prostate cancer behavior. *Cancer Cell*, **1**(2), 203–209.
#'
"prostate_ge"

#'
#' COIL-20 Object Image Dataset
#'
#' The COIL-20 dataset consists of grayscale images of 20 distinct objects. Each object
#' was placed on a turntable and imaged at pose intervals of 5 degrees, resulting in 72 images per object (totaling 1,440 images).
#'
#' Each row represents an image sample, while each column (except the last) corresponds to a pixel value.
#' The final column `Object` indicates the object's identity (from 1 to 20).
#'
#' @format A data frame with 1,440 rows and 1,025 columns:
#' \describe{
#'   \item{X1, X2, ..., X1024}{Grayscale pixel values (numeric).}
#'   \item{class}{Factor indicating the object label (1–20).}
#' }
#'
#' @references
#' Nene, S.A., Nayar, S.K., & Murase, H. (1996). Columbia Object Image Library (COIL-20). Technical Report CUCS-005-96, Department of Computer Science, Columbia University.
#'
"COIL20"
#
#'
#' warpAR10P Face Image Dataset
#'
#' The warpAR10P dataset consists of grayscale face images from 10 distinct individuals. Each individual's face was imaged under varying conditions, with each image warped for alignment, resulting in high-dimensional data.
#'
#' Each row represents an image sample, while each column (except the last) corresponds to a pixel value.
#' The final column `Person` indicates the subject's identity (from 1 to 10).
#'
#' @format A data frame with 130 rows and 2,401 columns:
#' \describe{
#'   \item{X1, X2, ..., X2400}{Grayscale pixel values (numeric).}
#'   \item{class}{Factor indicating the subject label (1–10).}
#' }
#'
#' @references
#' Li, J., Wang, K., et al. (scikit-feature). WarpAR10P available at: https://jundongl.github.io/scikit-feature/datasets.html
#'
"warpAR10P"
#
#'
#' warpPIE10P Face Image Dataset
#'
#' The warpPIE10P dataset consists of grayscale face images from 10 distinct individuals, sourced from the PIE face database. Each subject's face was captured under varying conditions, and all images were geometrically warped for alignment, resulting in high-dimensional feature data.
#'
#' Each row represents an image sample, while each column (except the last) corresponds to a pixel value.
#' The final column `Person` indicates the subject's identity (from 1 to 10).
#'
#' @format A data frame with 210 rows and 2,421 columns:
#' \describe{
#'   \item{X1, X2, ..., X2420}{Grayscale pixel values (numeric).}
#'   \item{class}{Factor indicating the subject label (1–10).}
#' }
#'
#' @references
#' Li, J., Wang, K., et al. (scikit-feature). WarpPIE10P datasets available at: https://jundongl.github.io/scikit-feature/datasets.html
#'
"warpPIE10P"
