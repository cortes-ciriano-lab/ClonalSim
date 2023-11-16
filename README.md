
# ClonalSim: A Python simulator for clonal expansions

<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/ClonalSim.png" width="80%">

ClonalSim is an object-oriented program that effectively simulates a Wright-Fisher model with selection. Specifically, it is used to calculate a clonal expansion where the number of mutant cells in each generation are drawn with binomial sampling. It also allows for switch of mutational rates, to make it appropriate for MMRD(Mismatch Repair Deficient) hypermutated samples. This code is being developed as part of the MSI evolution project.

Moreover ClonalSim can turn an array or matrix straight a python Tree object with nodes. Then, the user can use the class which is compatible with the treeswift python package for downstream analysis such as visualisation, editing etc.

### NOTE: Currently ClonalSim is under-development for other user usage and its used for research only by the Cortes-Ciriano lab. An upgrade is coming soon with add ons.

## Simulations with different fitness

<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/figure_clonalsim.png" alt="" width="80%">


## Functionalities
Main functionalities include:
1. Population class with features and a method inside to simulate Wright-Fisher
2. Transfer of an array population to a tree network with nodes, ready to import to the python treeswift package and versatile for further analysis.
3. Read in function for observed data from a tsv file or string.
4. Function to assign branch lengths based on a Poisson distribution with mean of a given mutation rate 
5.  Function to do averaging across the branch lengths of the tree before summary statistics (ultrametric tree transformation)

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
conda env create -n ClonalSim --file ClonalSim.yml
conda activate ClonalSim

```

## Example Run for simulating data
```
python clonalSim.py --N 2000 --generations 20 --mut_samples 60 --s 1.3 --mu 100  --sim_number 4
```
where N=population size, mu=mutation rate, mut_samples=number of mutant leaf nodes, s=fitness and sim_number is the number of simulations

## Example Code for reading observed data, calculating summary LTT and plotting
```
obs_tree = read_observed_data(obs_tree_path) # function that will read in the data to a tree class
obs_tree.draw()
normalise_tree_lengths(obs_tree) # normalise tree before LTT statistics
obs_tree.ltt()
```
## Example Observed Data Tree from MPN patient Van Egeren, D. et al. (2021) ‘Reconstructing the Lineage Histories and Differentiation Trajectories of Individual Cancer Cells in Myeloproliferative Neoplasms’, Cell stem cell
<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/et1tree.png" alt="Patient with ET" width="40%">

<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/ET1LTT.png" alt="Patient with ET" width="40%">


### Best Simulated Tree

<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/ET1bestsim.png" alt="Patient with ET" width="40%">


<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/example_results/ET1bestsimltt.png" alt="Patient with ET" width="40%">


## Approximate Bayesian Computation Using ClonalSim

<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/clonalsim_framew.png" alt="ClonalSim ABC Framework" width="80%">


## Authors

Maria Kalyva
[@mkalyva] (https://twitter.com/mariakalyva1)


## License

This project is licensed - see the LICENSE.txt file for details

## Acknowledgments
All the members of the Cortes-Ciriano lab, David Helekal for the discussions and the https://github.com/bacpop/PopPUNK group for the inspiration on how to use python treeswift classes.
