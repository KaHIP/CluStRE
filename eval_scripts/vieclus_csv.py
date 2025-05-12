import os
import csv

# Define constants for fixed values
ALGORITHM_NAME = "vieclus"
CLUSTER_FRAC = "1"

# Input and output paths
input_folder = "./experiments/with_vieclus/vieclus"  # Change to your folder path
output_csv = "experiments/with_vieclus/csv_files/vieclus.csv"

# Initialize the CSV file
with open(output_csv, mode='w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write header row
    writer.writerow(["alg_name", "graph", "cluster_frac", "runtime", "memory", "solution_quality"])

    # Process each .txt file in the folder
    for txt_file in os.listdir(input_folder):
        if txt_file.endswith(".txt"):
            with open(os.path.join(input_folder, txt_file), 'r') as file:
                lines = file.readlines()

                # Extract the graph name from the file name
                graph_name = txt_file.split("_output.txt")[0]

                # Extract required values from specific lines
                try:
                    time_spent = lines[5].split()[-1]  # Line 6 (index 5)
                    modularity = lines[6].split()[-1]  # Line 7 (index 6)
                    memory_consumption = lines[7].split()[-1]  # Line 8 (index 7)

                    # Write the row
                    writer.writerow([ALGORITHM_NAME, graph_name, CLUSTER_FRAC, time_spent, memory_consumption, modularity])
                except IndexError:
                    print(f"Error processing file: {txt_file}. File may not have expected format.")

print(f"CSV file '{output_csv}' created successfully!")
