import numpy as np
import random
import abcpy
import scipy.stats as stats



####################################################
## Simulate a Wright-Fisher Model with selection ###
####################################################


class Cell:
    def __init__(self, cell_type):
        self.type = cell_type

    def __str__(self):
        return f"Cell type: {self.type}"


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
            # to do

            #1 need to check the generation_data object
            #2 write a function to run the simulation
            # add cell class in the code
            # return 




class PhyTree:     # add code to create genealogy (basically lines between nodes and then randomly add the same mutation number)
    pass


def LTT_statistics(self, synth):
    """
    Calculate the LTT statistics for simulated trees.
    """
    


def distance_function(self, synth):
    """
    Calculate the distance between the observed and synthetic data.
    """

    return np.abs(self.obs - synth)



def estimate_parameters(self, epsilon):
    """
    Estimate the parameters of the Wright-Fisher model with selection using the ABC algorithm.
    """
    s_prior = abcpy.Distribution(stats.uniform, 0, 2)
    model = abcpy.Model(self.generate_synthetic_data, [self.N, s_prior])
    data = abcpy.Data([self.obs], [self.distance_function])
    sampler = abcpy.Sampler(model, data, epsilon=self.epsilon)
    results = sampler.sample(self.n_samples)
    s_estimates = results.get_parameters()[1]
    return np.mean(s_estimates)





        