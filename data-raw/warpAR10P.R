library(R.matlab)

# Load the MATLAB data file
warpAR10P <- readMat("data-raw/warpAR10P-raw.mat")
names(warpAR10P) <- c("x", "y")

# Save the processed data in the 'data/' folder for package use
usethis::use_data(warpAR10P, overwrite = TRUE, compress = "xz")

