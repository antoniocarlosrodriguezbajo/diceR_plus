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
#' The data includes expression levels for 7,129 genes across 72 patient samples.
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
#' Gene Expression Data for Lymphoma
#'
#' This dataset contains gene expression profiles from patients with
#' lymphoma. The data includes expression levels for 4026 genes across 62 patient samples.
#'
#' Each row represents a patient, while each column (except the last) corresponds
#' to a gene. The final column `class` classifies the patients as either `DLBCL`
#' (Diffuse large B-cell lymphoma), `FL` (Follicular lymphoma), or `CLL` (Chronic lymphocytic leukemia).
#'
#' @format A data frame with 62 rows and 4027 columns:
#' \describe{
#'   \item{X1, X2, ..., X4026}{Gene expression values (numeric).}
#'   \item{class}{Factor indicating the sample type (`DLBCL`, `FL`, or `CLL`).}
#' }
#'
#' @references
#' Alizadeh, A. A.; Eisen, M. B.; Davis, R. E.; Ma, C.; Lossos, I. S.; Rosenwald, A.; Boldrick, J. C.; Sabet, H.; Tran, T.; Yu, X.; et al. (2000). Distinct types of diffuse large B-cell lymphoma identified by gene expression profiling. *Nature*, **403**(6769), 503–511.
#'
"lymphoma"

#'
#' Gene Expression Data for Leukemia
#'
#' This dataset contains gene expression profiles from patients with
#' leukemia. The data includes expression levels for 999 genes across 38 patient samples.
#'
#' Each row represents a patient, while each column (except the last) corresponds
#' to a gene. The final column `class` classifies the patients as either `AML`
#' (acute myeloid leukemia), `ALL-B` (B-lineage acute lymphoblastic leukemia), or `ALL-T`
#' (T-lineage acute lymphoblastic leukemia).
#'
#' @format A data frame with 38 rows and 7,130 columns:
#' \describe{
#'   \item{AFFX.BioB.5_at, AFFX.BioB.M_at, ..., Z78285_f_at}{Gene expression values (numeric).}
#'   \item{class}{Factor indicating the leukemia subtype (`AML`, `ALL-B`, or `ALL-T`).}
#' }
#'
#' @references
#' Golub, T. R.; Slonim, D. K.; Tamayo, P.; Huard, C.; Gaasenbeek, M.; Mesirov, J. P.; Coller, H.; Loh, M. L.; Downing, J. R.; Caligiuri, M. A.; Bloomfield, C. D.; Lander, E. S. (1999). Molecular classification of cancer: class discovery and class prediction by gene expression monitoring. *Science*, **286**(5439), 531–537.
#'
"leukemia"

#' Gene Expression Data for Lung Cancer
#'
#' This dataset contains gene expression profiles from human lung tissue samples,
#' including both tumor and normal tissues. The data encompasses expression levels
#' for 3,313 genes across 203 samples. The samples represent four known classes:
#' Adenocarcinoma (AD), Squamous Cell Carcinoma (SQ), Carcinoid (COID), and
#' Normal Lung (NL). The AD class is known to be highly heterogeneous, and
#' substructure is believed to exist within this group.
#'
#' Each row represents a patient sample, while each column (except the last) corresponds
#' to a gene. The final column `class` indicates the sample type.
#'
#' @format A data frame with 203 rows and 3,313 columns:
#' \describe{
#'   \item{X1, X2, ..., X3313}{Gene expression values (numeric).}
#'   \item{class}{Factor indicating the tissue type (`AD`, `SQ`, `COID`, or `NL`).}
#' }
#'
#' @details
#' The samples comprise 139 adenocarcinomas, 21 squamous cell carcinomas, 20 carcinoids,
#' and 17 normal lung tissues. Gene selection was performed to enrich for class-distinguishing
#' genes prior to analysis.
#'
#' @references
#' Bhattacharjee, A., et al. (2001). Classification of human lung carcinomas by mRNA expression profiling reveals
#' distinct adenocarcinoma sub-classes. *Proceedings of the National Academy of Sciences*, **98**(24), 13790–13795.
#'
"lung_cancer"

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
#' Li, J., Wang, K., et al. (scikit-feature).
#' Dataset available at: https://jundongl.github.io/scikit-feature/datasets.html
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
#' Li, J., Wang, K., et al. (scikit-feature).
#' Dataset available at: https://jundongl.github.io/scikit-feature/datasets.html
#'
"warpPIE10P"
#
#' Meat Data
#'
#' This is the near-infrared spectroscopic meat data used in Murphy, Dean and Raftery (2009) <doi:10.1214/09-AOAS279> and originally collected by McElhinney, Downey and Fearn (1999) <doi:10.1255/jnirs.245>.
#'
#' @docType data
#'
#' @usage data(Meat)
#'
#' @format A list with two components:
#' \describe{
#' \item{x}{Homogenized raw meat spectra. A matrix with 231 rows and 1050 columns.}
#' \item{y}{A vector containing the true class memberships.}}
#'
#' @keywords datasets
#'
#' @references Murphy, Dean and Raftery (2010) <doi:10.1214/09-AOAS279>
#'
#' @source McElhinney, Downey and Fearn (1999) <doi:10.1255/jnirs.245>
#'
#' @examples
#' data(Meat)
#' Meat$x[1:5,1:5]
#' Meat$y
"Meat"
