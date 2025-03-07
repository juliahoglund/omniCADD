import os
from data_helper import load_dataset, save_npz_with_meta

# Define test data paths
test_files = ["path/to/test_data1", "path/to/test_data2"]
desired_columns = "path/to/desired_columns.csv"

# Load dataset
data_matrix, y_values, columns = load_dataset(test_files, desired_columns)

# Save dataset
output_file = "path/to/output_data.npz"
save_npz_with_meta(output_file, data_matrix, y_values, columns)

print("Test completed successfully.")
