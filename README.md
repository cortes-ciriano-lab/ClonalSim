
# ClonalSim: A Python Wright-Fisher simulator for clonal expansions
Use for the MSI evolution paper and thesis chapters. 
## Description of the code
#### Master Simulator Function that includes running the following:
Step 1:  Python class that initiates a population class with features and a simulation method

Step 2:  Function that runs a population simulation and returns population, figures, probabilities used in each generation

Step 3:  Function that from a population array create a dictionary of clusters which represent clades of the tree

Step 4:  Function that takes in a population array, turns it in a network of nodes and then to a tree object 
         compatible with the treeswift python package for further manipulaton in an efficient way

Step 5:  Function that assigns edge lengths to the simulated tree based on a Poisson distribution with mean μ

Step 6:  Function that calculates the Lineage Through Time plots and statistics for a simulated tree

Step 7: function that reads in the observed data in a treeswift tree and returns Lineage Through Time plots and statistics

### Dependencies

* python3, miniconda
* all miniconda libraries can be installed with 

```
conda env create -n ColoSim --file ColoSim.yml
conda activate ColoSim

```

## Example Run for simulating data
```
python3 coloSim.py --N 2000 --generations 20 --mut_samples 60 --s 1.3 --mu 100
```
where N=population size, mu=mutation rate, mut_samples=number of mutant leaf nodes and s=fitness

## Example Run for reading observed data
```
python3 coloSim.py --N 2000 --generations 20 --mut_samples 60 --s 1.3 --mu 100
```


## Example Results of the Wright-Fisher model simulation

#### Example Results with high fitness
<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/output_s_1.4.png" alt="Results tree with s=1.8" width="80%">
<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/output_ltt_1.4.png" alt="Results tree with s=1.8" width="80%">
<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/output_tree_1.4.png" alt="Results tree with s=1.8" width="80%">

#### Example Results with low fitness
<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/output_s_0.5.png" alt="Results tree with s=0.5" width="80%">
<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/output_ltt_0.5.png" alt="Results tree with s=0.5" width="80%">
<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/output_tree_0.5.png" alt="Results tree with s=0.5" width="80%">



## Authors

Maria Kalyva
[@mkalyva] (https://twitter.com/mariakalyva1)


## License

This project is licensed - see the LICENSE.txt file for details

## Acknowledgments
All the members of the Cortes-Ciriano lab for the discussions and https://github.com/bacpop/PopPUNK for inspiration on how to use python treeswift classes.
