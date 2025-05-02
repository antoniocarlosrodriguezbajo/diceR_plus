# Load required package
library(data.table)  # Efficient CSV reading for large datasets

# Read the CSV file, treating the first row as column names
meats <- fread("data-raw/meats-raw.txt", header = TRUE)

# Convert 'type' column to a factor
meats$type <- as.factor(meats$type)

# Save the processed data in the 'data/' folder for package use
usethis::use_data(meats, overwrite = TRUE, compress = "xz")
