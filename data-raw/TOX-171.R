library(R.matlab)

# Load the MATLAB data file
TOX_171 <- readMat("data-raw/TOX-171-raw.mat")
names(TOX_171) <- c("y", "x")

# Save the processed data in the 'data/' folder for package use
usethis::use_data(TOX_171, overwrite = TRUE, compress = "xz")
