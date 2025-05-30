library(R.matlab)

# Load the MATLAB data file
prostate_ge <- readMat("data-raw/prostate_ge-raw.mat")
names(prostate_ge) <- c("x", "y")

# Save the processed data in the 'data/' folder for package use
usethis::use_data(prostate_ge, overwrite = TRUE, compress = "xz")

