# Load required package
library(data.table)  # Efficient CSV reading for large datasets

# Read the CSV file, treating the first row as column names
meats <- fread("data-raw/meats-raw.txt", header = TRUE)

# Rename 'type' to 'class'
setnames(meats, "type", "class")

# Convert 'class' column to a factor
meats$class <- as.factor(meats$class)

# Save the processed data in the 'data/' folder for package use
usethis::use_data(meats, overwrite = TRUE, compress = "xz")
