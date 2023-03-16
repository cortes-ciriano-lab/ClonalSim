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
        num_generations = len(p1.generation_data)
        for i in range(num_generations - 1):
            curr_gen = p1.generation_data[i]
            next_gen = p1.generation_data[i + 1]
            for j in range(p1.N):
                if curr_gen[j] == 1:
                    # connect to a random cell in the next generation
                    next_gen_indices = np.argwhere(next_gen == 1).flatten()

                    ### PSEUDO-CODE / Googled CODE ## 
                    if len(next_gen_indices) > 0:
                        next_gen_index = np.random.choice(next_gen_indices)
                        node1 = f"{i}_{j}"
                        node2 = f"{i+1}_{next_gen_index}"
                        self.add_edge(node1, node2)

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
    return tree


# ------------------------------------------- #



# Calculate LTT stats from simulated data
def LTT_statistics(synth):
    """
    Calculate the LTT statistics for simulated trees.
    """
    

# Distance function between observed and synthetic data
def read_observed_data(observed):
    """
    Read observed tree and calculate LTT statistics and return or read in LTT statistics straight
    """

    return LTT_stat_obs


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








        