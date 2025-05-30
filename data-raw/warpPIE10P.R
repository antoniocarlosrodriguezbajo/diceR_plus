library(R.matlab)

# Load the MATLAB data file
warpPIE10P <- readMat("data-raw/warpPIE10P-raw.mat")
names(warpPIE10P) <- c("x", "y")

# Save the processed data in the 'data/' folder for package use
usethis::use_data(warpPIE10P, overwrite = TRUE, compress = "xz")
