library(R.matlab)

# Load the MATLAB data file
warpAR10P_data <- readMat("data-raw/warpAR10P-raw.mat")

# Convert X (image features) into a data frame
warpAR10P <- as.data.frame(warpAR10P_data$X)

# Rename columns to X1, X2, ..., X2400
colnames(warpAR10P) <- paste0("X", seq_len(ncol(warpAR10P)))

# Extract class labels (Y) separately
Person <- factor(warpAR10P_data$Y)

# Add Person column to the dataframe
warpAR10P$Person <- Person

# Save the processed data in the 'data/' folder for package use
usethis::use_data(warpAR10P, overwrite = TRUE, compress = "xz")
