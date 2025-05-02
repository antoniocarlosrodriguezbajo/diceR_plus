library(R.matlab)

# Load the MATLAB data file
prostate_ge_data <- readMat("data-raw/prostate_ge-raw.mat")

# Convert X (gene expression) into a data frame
prostate_ge <- as.data.frame(prostate_ge_data$X)

# Rename columns to X1, X2, ..., X5966
colnames(prostate_ge) <- paste0("X", seq_len(ncol(prostate_ge)))

# Extract class labels (Y) separately
Diagnosis <- factor(prostate_ge_data$Y, labels = c("normal", "tumor"))

# Add Diagnosis column to the dataframe
prostate_ge$Diagnosis <- Diagnosis

# Save the processed data in the 'data/' folder for package use
usethis::use_data(prostate_ge, overwrite = TRUE, compress = "xz")
