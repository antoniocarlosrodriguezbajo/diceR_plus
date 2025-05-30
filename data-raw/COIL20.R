library(R.matlab)

# Load the MATLAB data file
COIL20 <- readMat("data-raw/COIL20-raw.mat")
names(COIL20) <- c("x", "y")

# Save the processed data in the 'data/' folder for package use
usethis::use_data(COIL20, overwrite = TRUE, compress = "xz")
