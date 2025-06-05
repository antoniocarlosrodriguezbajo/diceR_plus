# Read gene expression
lung_cancer_var <- read.table("data-raw/lung_cancer-raw.txt", header = FALSE)

# Read and assign class
class <- scan("data-raw/lung_cancer_labels-raw.txt", what = character(), sep = "\t")

lung_cancer <- list(
  x = as.matrix(lung_cancer_var),  # Convierte el data frame en una matriz numérica
  y = as.factor(class)         # Convierte la variable de clase en un factor
)

# Save the processed data in the 'data/' folder for package use
usethis::use_data(lung_cancer, overwrite = TRUE, compress = "xz")
