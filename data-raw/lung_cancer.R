# Read gene expression
lung_cancer <- read.table("data-raw/lung_cancer-raw.txt", header = FALSE)
colnames(lung_cancer) <- paste0("X", 1:ncol(lung_cancer))

# Read and assign class
class <- scan("data-raw/lung_cancer_labels-raw.txt", what = character(), sep = "\t")
lung_cancer$class <- factor(class)

# Save the processed data in the 'data/' folder for package use
usethis::use_data(lung_cancer, overwrite = TRUE, compress = "xz")
