library(R.matlab)

# Load the MATLAB data file
GLIOMA <- readMat("data-raw/GLIOMA-raw.mat")
names(GLIOMA) <- c("x", "y")

# Save the processed data in the 'data/' folder for package use
usethis::use_data(GLIOMA, overwrite = TRUE, compress = "xz")
