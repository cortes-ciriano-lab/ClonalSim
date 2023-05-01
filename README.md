
# Wright-Fisher with selection python simulator
Use for the MSI evolution paper and thesis chapters. 

<img src="https://github.com/cortes-ciriano-lab/ColoSim/blob/main/MSI_dynamics.jpg" alt="MSI dynamics" width="20%">

## Description of the code
# Master Simulator Function that includes running the following:
Step 1:  initiate a population class with features and a simulation method

Step 2:  run population simulation and return population, figure, probabilities used in each generation

Step 3:  from a population array create a dictionary of clusters which represent clades of the tree

Step 4:  function that takes in a population array, turns it in a network of nodes and then to a tree object 
         compatible with the treeswift python package for further manipulaton in an efficient way

Step 5:  function that assigns edge lengths to the simulated tree based on a Poisson distribution with mean μ

Step 6:  calculation of Lineage Through Time plots and statistics for simulated tree

Step 7: function that reads in the observed data in a treeswift tree and returns Lineage Through Time plots and statistics

### Dependencies

* python3, miniconda
* all miniconda libraries can be installed with 

```
conda env create -n ColoSim --file ColoSim.yml
conda activate ColoSim

```

## Example Results of the Wright-Fisher model simulation


## Authors

Contributors names and contact info

ex. Maria Kalyva
ex. [@mkalyva] (https://twitter.com/mariakalyva1)


## License

This project is licensed under License - see the LICENSE.txt file for details

## Acknowledgments
All the members of the Cortes-Ciriano lab for the discussions and https://github.com/bacpop/PopPUNK for inspiration on how to use python treeswift classes.
