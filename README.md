# In and Out: Optimizing Overall Interaction in Probabilistic Graphs under Clustering Constraints



## Folders
- datasets:  contains the original data as well as the preprocessed data (as described in the paper). Biggest networks (to be unzipped in datasets folder) can be downloaded at the following [link](https://drive.google.com/open?id=1r0krGQMm0QyUAbJlGZSxTmm2ULaajBKI)
- code: it contains this project code and links to the competing methods 
- output: it stores all results produced by the algorithms 

## Usage

From the folder 'optimize_interactions/', run the following command:
```bash
run_offline_algorithm.py [-h] -d DATASET -a {MIL,D-MIL} [-s SEED]
                                [-r RUNS] [--hill_climbing [HILL_CLIMBING]]
                                [-I HILL_CLIMBING_ITERS]
                                [--save_only_avg [SAVE_ONLY_AVG]]
```
### Dependencies
- numpy
- python-igraph


### Positional arguments
```
  -d DATASET, --dataset DATASET
                        Input dataset, whose name identifies a particular subfolder in 'datasets/'
  -a {MIL,D-MIL}, --alg {MIL,D-MIL}
                        Selected algorithm 
```
### Optional arguments
```
  -s SEED, --seed SEED  Random generation seed -- for reproducibility (default value 100)
  -r RUNS, --runs RUNS  Number of runs of the selected randomized algorithm (default value 50) 
  --hill_climbing [HILL_CLIMBING]
                        It enables the hill climbing postprocessing step
  -I HILL_CLIMBING_ITERS, --hill_climbing_iters HILL_CLIMBING_ITERS
                        Maximum number of iterations for the hill climbing
                        algorithm (default value 8)
  --save_only_avg [SAVE_ONLY_AVG]
                        if True (default value), only average results across several runs are
                        stored, otherwise the results for each run is saved in a specific sub-
                        folder in 'output/'
                    
```
Default values for -s, -r and -I correspond to the settings relating to the results presented in the paper.
