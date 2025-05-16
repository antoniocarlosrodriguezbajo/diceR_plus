# Load package
library(RPEClust)

# Load the data
data(Meat)

# Save the processed data in the 'data/' folder for package use
usethis::use_data(Meat, overwrite = TRUE, compress = "xz")
