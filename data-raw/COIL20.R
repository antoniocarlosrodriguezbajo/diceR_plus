library(R.matlab)

# Load the MATLAB data file
COIL20_data <- readMat("data-raw/COIL20-raw.mat")

# Convert X (image features) into a data frame
COIL20 <- as.data.frame(COIL20_data$X)

# Rename columns to X1, X2, ..., X1024
colnames(COIL20) <- paste0("X", seq_len(ncol(COIL20)))

# Extract class labels (Y) separately
class <- factor(COIL20_data$Y)

# Add Diagnosis column to the dataframe
COIL20$class <- class

# Save the processed data in the 'data/' folder for package use
usethis::use_data(COIL20, overwrite = TRUE, compress = "xz")
