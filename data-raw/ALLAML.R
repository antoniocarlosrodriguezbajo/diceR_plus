library(R.matlab)

# Load the MATLAB data file
ALLAML <- readMat("data-raw/ALLAML-raw.mat")
names(ALLAML) <- c("x", "y")

# Save the processed data in the 'data/' folder for package use
usethis::use_data(ALLAML, overwrite = TRUE, compress = "xz")
