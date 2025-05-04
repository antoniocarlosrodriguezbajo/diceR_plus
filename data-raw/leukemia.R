# Load the package
library(golubEsets)

# Load the data
data(Golub_Train)

# Extract the phenotypic data
pheno <- pData(Golub_Train)

# Create the class column
class <- ifelse(
  pheno$ALL.AML == "AML",
  "AML",
  ifelse(pheno$T.B.cell == "T-cell", "T-ALL", "B-ALL")
)

# Extract the gene expression
exp_data <- t(exprs(Golub_Train))

# Create the dataframe
leukemia <- data.frame(exp_data)
leukemia$class = factor(class)

# Save the processed data in the 'data/' folder for package use
usethis::use_data(leukemia, overwrite = TRUE, compress = "xz")
