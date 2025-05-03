# Load the RData file
load("data-raw/lymphoma-raw.RData")

# Convert X (gene expression) into a data frame
lymphoma_data <- lymphoma
lymphoma <- as.data.frame(lymphoma_data$x)

# Rename columns to X1, X2, ..., X4026
colnames(lymphoma) <- paste0("X", seq_len(ncol(lymphoma)))

# Extract class labels (Y) separately
class <- factor(lymphoma_data$y, labels = c("DLBCL", "FL", "CLL"))

# Add Diagnosis column to the dataframe
lymphoma$class <- class

# Save the processed data in the 'data/' folder for package use
usethis::use_data(lymphoma, overwrite = TRUE, compress = "xz")
