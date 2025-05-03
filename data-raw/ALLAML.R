library(R.matlab)

# Load the MATLAB data file
ALLAML_data <- readMat("data-raw/ALLAML-raw.mat")

# Convert X (gene expression) into a data frame
ALLAML <- as.data.frame(ALLAML_data$X)

# Rename columns to X1, X2, ..., X7130
colnames(ALLAML) <- paste0("X", seq_len(ncol(ALLAML)))

# Extract class labels (Y) separately
class <- factor(ALLAML_data$Y, labels = c("AML", "ALL"))

# Add Diagnosis column to the dataframe
ALLAML$class <- class

# Save the processed data in the 'data/' folder for package use
usethis::use_data(ALLAML, overwrite = TRUE, compress = "xz")
