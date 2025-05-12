#!/bin/bash

# Directory containing graph files
graph_test_dir="/home/speretz/experiments/testing_graphs/"
graph_tune_dir="/home/speretz/experiments/tuning_graphs/"

# Modes and their respective output directories
modes=("light_plus" "light" "evo" "strong")

#modes=("light")
# Output file to store commands
output_file="par_commands.txt"

# Clear the output file if it already exists
echo -n "" > "$output_file"

# List of graph files

graph_test=(
    "uk-2007-05.graph"                                          #web        SEA[9,10,11]        105,896,555         3,301,876,564
    "webbase-2001.graph"                                        #web        SEA[9,10,11]        118,142,155           854,809,761
    "sk-2005-sorted.graph"                                      #web        SEA[9,10,11]         50,636,154         1,810,063,330
    "eu-2005.graph"                                             #web        SEA[9,10,11]            862,664            1,613,8468
    "arabic-2005.graph"                                         #web        SEA[9,10,11]         22,744,080           553,903,073
    "it-2004.graph"                                             #web        SEA[9,10,11]         41,291,594         1,027,474,947
    "wiki-Talk.graph"                                           #web        MA Adil                 232,314             1,458,806

    "com-friendster.graph"                                      #social     SNAP                 65,608,366         1,806,067,135
    "com-amazon.graph"                                          #social     SNAP                    334,863               925,872
    "hollywood-2011.graph"                                      #social     BA wilwert            2,180,759           114,492,816
    "libimseti-sorted"                                          #social     BA wilwert              220,970            17,233,144    
    "Penn94.graph"                                              #social     BA wilwert               41,554             1,362,229
    "enwiki-2013.graph"                                         #social     BA wilwert            4,206,785            91,939,728

    "great-britain.osm.graph"                                   #road       DIMACS                7,733,822             8,156,517
    "italy.osm.graph"                                           #road       DIMACS                6,686,493             7,013,978
    "europe.osm.graph"                                          #road       DIMACS               50,912,018            54,054,660
    
    "coPapersDBLP.graph"                                        #Citations    MA adil               540,486            1,524,5729
    "citationCiteseer.graph"                                    #citations    MA Adil               268,495             1,156,647

    "rhg2b.graph"                                               #Rand Hyp.  KAGEN               100,000,000         1,999,544,833
    "rgg_n26.graph"                                             #Rand. Geo. KAGEN                67,108,864           574,553,645
)

graphs_tune=(
    "in-2004.graph"                                             #web       SEA[9,10,11]                     1,382,908       13,591,473
    "web-Google.graph"                                          #web       SEA[9,10,11]                       356,648        2,093,324
    "amazon-2008.graph"                                         #web       SEA[9,10,11]                       735,323        3,523,472
    
    "ljournal-2008.graph"                                       #social     MA Adil                         5,363,260       49,514,271 
    "Texas84.graph"                                             #social     BA wilwert                         36,371        1,590,655
    "Maryland58.graph"                                          #social     Network repository                 20,871          744,862

    "coAuthorsDBLP.graph"                                       #citations  MA Adil                           299,067          977,676

    "netherlands.osm.graph"                                     #road       DIMACS                          2,216,688        2,441,238
    "germany.osm.graph"                                         #road       DIMACS                         11,548,845       12,369,181

    "rhg1m10m.graph"                                            #Rand Geo.                                  1,000,000       10,047,330
)

# Loop through each mode and graph

for graph in "${graph_test[@]}"; do
    for mode in "${modes[@]}"; do
        input_graph="$graph_test_dir$graph"
        output_path="/home/speretz/experiments/clustre_test/$mode/"
        command="./build/clustre $input_graph --one_pass_algorithm=modularity --evaluate --ext_clustering_algorithm=VieClus --ext_algorithm_time=15 --mode=$mode --output_path=$output_path"
        echo "$command" >> "$output_file"
    done
done

# Inform user of completion
echo "Commands written to $output_file"

