import os
import csv
import flatbuffers
from StreamCPIInfo import Partition, GraphMetadata, PartitionConfiguration, RunTime, PartitionMetrics, MemoryConsumption


def read_partition_from_file(file_path):
    """Reads a FlatBuffer .bin file and extracts relevant statistics."""
    with open(file_path, 'rb') as file:
        buf = file.read()
        partition = Partition.Partition.GetRootAsPartition(buf, 0)

        # Extract runtime information
        runtime = partition.Runtime()
        total_time = runtime.TotalTime()

        # Extract clustering metrics
        clustering_metrics = partition.ClusteringMetrics()
        score = clustering_metrics.Score()

        # Extract memory consumption
        memory_consumption = partition.MemoryConsumption()
        overall_max_rss = memory_consumption.OverallMaxRss()

        # Extract graph metadata
        graph_metadata = partition.GraphMetadata()
        graph_name = graph_metadata.Filename().decode('utf-8')

        # Extract cluster fraction
        partition_config = partition.PartitionConfiguration()
        cluster_frac = partition_config.ClusterFrac()

        return graph_name, total_time, overall_max_rss, score, cluster_frac


def process_flatbuffer_files(folder_path, output_csv_path):
    """Processes all FlatBuffer files in a folder and writes results to a CSV."""
    # CSV headers
    headers = ["alg_name", "graph", "cluster_frac", "runtime", "memory", "solution_quality"]

    # Extract algorithm name from the CSV file name
    alg_name = os.path.splitext(os.path.basename(output_csv_path))[0]

    # Prepare to write the output CSV
    with open(output_csv_path, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)  # Write the headers

        # Iterate through all files in the folder
        for filename in os.listdir(folder_path):
            if filename.endswith(".bin"):  # Adjust extension if needed
                file_path = os.path.join(folder_path, filename)
                try:
                    # Extract data from the FlatBuffer file
                    graph_name, runtime, memory, solution_quality, cluster_frac = read_partition_from_file(file_path)

                    # Write the row to the CSV
                    writer.writerow([alg_name, graph_name, cluster_frac, runtime, memory, solution_quality])
                except Exception as e:
                    print(f"Error processing file {filename}: {e}")


# Specify your folder and output CSV file paths
folder_path = "experiments/restream/no_vieclus/std/1/"  # Replace with your folder path
output_csv_path = "experiments/restream/csv_files/no_vieclus.csv"  # Replace with desired output CSV file path

# Process the files and generate the CSV
process_flatbuffer_files(folder_path, output_csv_path)

print(f"CSV file has been created at {output_csv_path}.")
