# Load package
library(spls)

# Load the data
data(lymphoma)

# Save the processed data in the 'data/' folder for package use
usethis::use_data(lymphoma, overwrite = TRUE, compress = "xz")

