import numpy as np
import random
import abcpy
import scipy.stats as stats



####################################################
## Simulate a Wright-Fisher Model with selection ###
####################################################


# class Cell:
#     def __init__(self, cell_type, position):
#         self.type = cell_type
#         self.position = position

#     def __str__(self):
#         return f"Cell type: {self.type}"


class Population:
    def __init__(self, N, generations, mut_samples, s):
        self.N = N
        self.s = s
        self.generations = generations
        self.mut_samples = mut_samples
        self.generation_data = []

    def __str__(self):
        return f"Population size: {self.N}, Generations: {self.generations}, Mutant Samples: {self.mut_samples}"

    def simulate_population(self, params):
        """
        Simulate the population using the Wright-Fisher model with selection.
        """
        # Initialize the first population

        # Initialize the first population
        population = np.zeros(self.N, dtype=int)
        population[random.randint(0, self.N - 1)] = 1

        for gen in range(self.generations):
            cancer_p = (1 + s) * n / (N + n * s)
            offspring = np.random.binomial(n=1, p=cancer_p, size=N)
            self.generation_data.append(offspring)
            

class PhyTree:

    def __init__(self):
        self.graph = {}
        self.node_labels = {}

    def connect_random_cells(self, generation_data):
        """
        Connects random cells of the generations in generation_data.
        """

        # Initialize genealogy list where for each mutated cell I save the indexes of its randomly assigned genealogy
        genealogy = []

            if gen == self.generations[len(self.generations)]:
                last_gen_indices = np.where(gen == 1)[0]
                random_indices = np.random.choice(last_gen_indices, size=10)
                for index in random_indices:
                    curr_gen_index = index
                    indices_list = [curr_gen_index]
                    for prev_gen in range(gen - 1, -1, -1):
                        prev_gen_indices = np.where(self.generation_data[prev_gen] == 1)[0]
                        if len(prev_gen_indices) > 0:
                            prev_gen_index = np.random.choice(prev_gen_indices)
                            indices_list.append(prev_gen_index)
                            curr_gen_index = prev_gen_index
                        else:
                            break
                    genealogy.append(indices_list)

                return genealogy



    def assign_edge_lengths(self, mu):
        """
        Assigns a length in that edge between the nodes based on a Poisson distribution with mean rate μ.
        """
        for edge in self.graph.edges():
            length = np.random.poisson(mu)
            self.graph.edges[edge]["length"] = length


# ------------- Run Simulation ------------- #

# to add: for loop for multiple iterations & all simulations in list

def simulate_population_and_tree(N, generations, mut_samples, s, mu):
    pop = Population(N, generations, mut_samples, s)
    pop.simulate_population({})
    tree = PhyTree()
    tree.connect_random_cells(pop.generation_data)
    tree.assign_edge_lengths(mu)
    # write phylogenetic tree to file
    
    return tree


# ------------------------------------------- #

# Calculate LTT stats from simulated data
def LTT_statistics(synth):
    """
    Calculate the LTT statistics for simulated trees.
    """


def read_observed_data(observed):
    """
    Read observed tree and calculate LTT statistics and return or read in LTT statistics straight
    """
    # read tree from tsv file
    import treeswift

    # Define the path to the TSV file containing the tree
    tree_file = "path/to/tree.tsv"

    # Load the tree from the TSV file
    with open(tree_file) as f:
        tree_str = f.read()
    tree = treeswift.read_tree_newick(tree_str)

    # Calculate lineage through time plot statistics
    ltt = tree.lineage_through_time()

    # Save the results in a data structure
    results = {
        "tree": tree,
        "ltt": ltt,
    }
    print("Tree:", results["tree"])
    print("LTT statistics:", results["ltt"])
    return tree,ltt

# Distance function between observed and synthetic data
def distance_function():
    """
    Calculate the distance between the observed and synthetic data.
    """

    return np.abs(self.obs - synth)



def estimate_parameters(epsilon):
    """
    Estimate the parameters of the Wright-Fisher model with selection using the ABC algorithm.
    
    ### PSEUDO-CODE ####

    """
    s_prior = abcpy.Distribution(stats.uniform, 0, 2)
    model = abcpy.Model(generate_synthetic_data, [N, s_prior])
    data = abcpy.Data([obs], [distance_function])
    sampler = abcpy.Sampler(model, data, epsilon=epsilon)
    results = sampler.sample(n_samples)
    s_estimates = results.get_parameters()[1]
    return np.mean(s_estimates)
        