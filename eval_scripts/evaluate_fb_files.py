import sys
import os
import flatbuffers
from StreamCPIInfo import Partition, GraphMetadata, PartitionConfiguration, RunTime, PartitionMetrics, MemoryConsumption

def read_partition_from_file(file_path):
    """Reads a FlatBuffer .bin file and extracts relevant statistics."""
    with open(file_path, 'rb') as file:
        buf = file.read()
        partition = Partition.Partition.GetRootAsPartition(buf, 0)

        # Extract the required fields
        runtime = partition.Runtime()
        total_time = runtime.TotalTime()

        clustering_metrics = partition.ClusteringMetrics()
        score = clustering_metrics.Score()

        memory_consumption = partition.MemoryConsumption()
        overall_max_rss = memory_consumption.OverallMaxRss()

    return total_time, score, overall_max_rss

def process_all_flatbuffers(base_dir):
    """Processes all FlatBuffer .bin files in the given directory structure."""
    total_time_sum = 0
    score_sum = 0
    overall_max_rss_sum = 0
    file_count = 0

    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".bin"):  # Process only .bin files
                file_path = os.path.join(root, file)
                try:
                    total_time, score, overall_max_rss = read_partition_from_file(file_path)
                    total_time_sum += total_time
                    score_sum += score
                    overall_max_rss_sum += overall_max_rss
                    file_count += 1
                except Exception as e:
                    print(f"Error processing file {file_path}: {e}")

    if file_count == 0:
        print("No .bin files found.")
        return

    # Calculate averages
    avg_total_time = total_time_sum / file_count
    avg_score = score_sum / file_count
    avg_overall_max_rss = overall_max_rss_sum / file_count

    print(f"Processed {file_count} .bin files.")
    print(f"Average Total Time: {avg_total_time}")
    print(f"Average Score: {avg_score}")
    print(f"Average Overall Max RSS: {avg_overall_max_rss}")

if __name__ == "__main__":
    base_dir = "experiments/std/"  # Set this to the root directory containing the folders (e.g., cmp, std)
    process_all_flatbuffers(base_dir)
