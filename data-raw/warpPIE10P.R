library(R.matlab)

# Load the MATLAB data file
warpPIE10P_data <- readMat("data-raw/warpPIE10P-raw.mat")

# Convert X (image features) into a data frame
warpPIE10P <- as.data.frame(warpPIE10P_data$X)

# Rename columns to X1, X2, ..., X2420
colnames(warpPIE10P) <- paste0("X", seq_len(ncol(warpPIE10P)))

# Extract class labels (Y) separately
Person <- factor(warpPIE10P_data$Y)

# Add Person column to the dataframe
warpPIE10P$Person <- Person

# Save the processed data in the 'data/' folder for package use
usethis::use_data(warpPIE10P, overwrite = TRUE, compress = "xz")
