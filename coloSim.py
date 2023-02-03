import numpy as np
import abcpy
import scipy.stats as stats

class WrightFisherModelWithSelection:
    def __init__(self, N, generations, obs, epsilon, n_samples):
        self.N = N
        self.generations = generations
        self.obs = obs
        self.epsilon = epsilon
        self.n_samples = n_samples

    def generate_summary_statistics(self, population):
        """
        Generate summary statistics from the population.
        """
        return np.mean(population) / 2

    def generate_synthetic_data(self, params):
        """
        Generate synthetic data using the Wright-Fisher model with selection.
        """
        N, s = params
        population = np.random.binomial(2, 0.5, N)
        for t in range(self.generations):
            avg_fitness = 1 + s * np.mean(population)
            offspring = np.random.binomial(2, 0.5 * avg_fitness / (1 + s), N)
            population = offspring
        return self.generate_summary_statistics(population)

    def distance_function(self, synth):
        """
        Calculate the distance between the observed and synthetic data.
        """
        return np.abs(self.obs - synth)

    def estimate_parameters(self):
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
        