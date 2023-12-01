# Abstract
Maximal biclique enumeration (MBE) in bipartite graphs is a fundamental problem in data mining with widespread applications.
Many recent works solve this problem based on the set-enumeration (SE) tree,
which sequentially traverses vertices to generate the enumeration tree nodes representing distinct bicliques,
then checks whether these bicliques are maximal or not.
However, existing MBE algorithms only expand bicliques with untraversed vertices to ensure distinction,
which often necessitate extensive node checks to eliminate non-maximal bicliques,
resulting in significant computational overhead during the enumeration process.
To address this issue, we propose an aggressive set-enumeration (ASE) tree that
aggressively expands all bicliques to their maximal form,
thus avoiding costly node checks on non-maximal bicliques.
This aggressive enumeration may produce multiple duplicate maximal bicliques,
but we efficiently eliminate these duplicates by
leveraging the connection between parent and child nodes 
and conducting low-cost node checking.
Additionally, we introduce an aggressive merge-based pruning (AMP) approach that
aggressively merges vertices sharing the same local neighbors.
This helps prune numerous duplicate node generations caused by subsets of merged vertices.
We integrate the AMP approach into the ASE tree,
and present the Aggressive Maximal Biclique Enumeration Algorithm (AMBEA).
Experimental results show that
AMBEA is 1.15X to 5.32X faster than its closest competitor
and exhibits better scalability and parallelization capabilities on larger bipartite graphs.
# Try out AMBEA
## Software Dependencies
- GNU Make 4.2.1
- CMake 3.22.0
- GCC/G++ 10.3.0
- Python 3.8.10

## Downloading
To download this repository, you can run the following commands.
```
wget https://anonymous.4open.science/r/MMBEA-public-D206/zipfile/AMBEA-public.zip
unzip AMBEA-public.zip
```

## Compiling
Using the following commands, one can easily compile the AMBEA and baselines. Then you will find the executable files under `bin/MBE_ALL`.
```
bash ./scripts/compile-ALL.sh
```

## Running

You can run AMBEA or baselines with the following command-line options.
```
./bin/MBE_ALL 
 -i: The path of input dataset file.
 -s: Select one GMBE version to run. 1: MBEA, 2: iMBEA, 3: FMBE, 4: PMBE, 5: AMBEA

#Execute this command to run ooMBEA on the specific dataset
./bin/mbbp [path of dataset]
```
## Experimental workflow
We provide a script to generate the experimental results of Figure 6 in the directory `scripts/`. You can execute the script as following.
```
bash ./scripts/overall-evaluation.sh 
```
Then you will find the experimental results in the file `scripts/results.txt`.
