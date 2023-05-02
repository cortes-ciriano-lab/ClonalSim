import matplotlib.pyplot as plt
import numpy as np
import random
#import abcpy
import treeswift
from treeswift import Tree, Node
import argparse

# create an argparse parser
parser = argparse.ArgumentParser(description="Simulate population and tree")

# add arguments for the simulation parameters
parser.add_argument("--N", type=int, help="population size")
parser.add_argument("--generations", type=int, help="number of generations to simulate")
parser.add_argument("--mut_samples", type=int, help="number of mutation samples")
parser.add_argument("--s", type=float, help="selection coefficient")
parser.add_argument("--mu", type=float, help="mutation rate")

# parse the command-line arguments
args = parser.parse_args()

class Population:
    def __init__(self, N, generations, s):
        self.N = N
        self.s = s
        self.generations = generations
        self.generation_data = []

    def __str__(self):
        return f"Population size: {self.N}, Generations: {self.generations}, Selection: {self.s}"

    def simulate_population(self):
        """
        Simulate the population using the Wright-Fisher model with selection.
        """

        # Initialize the first population
        population = np.zeros(self.N, dtype=int)
        population[random.randint(0, self.N - 1)] = 1
        self.generation_data.append(population)
        binom_prob_list  = []
        mut_n_list = []
        
        for gen in range(self.generations):
            mut_n = len(np.where(self.generation_data[gen] == 1)[0])
            mut_n_list.append(mut_n)
            
            cancer_p = (1 + self.s) * mut_n / (self.N + (mut_n * self.s))
            binom_prob_list.append(cancer_p)
                
            offspring = np.random.binomial(n=1, p=cancer_p, size=self.N)
            
            num_mutants = [np.count_nonzero(offspring == 1)]
            
            if num_mutants == 0:
                print("Stochastic Extinction")
                self.generation_data.append(offspring)
                num_mutants = [np.count_nonzero(generation == 1) for generation in self.generation_data]
                # Plot the number of mutants over time
                fig, ax = plt.subplots()
                ax.plot(range(len(num_mutants)), np.log(num_mutants))
                ax.set_xlabel("Time in Generations")
                ax.set_ylabel("Number of mutants ln(N)")
                ax.set_title(f"Mutant allele frequency over time (s={self.s})")
                plt.show()
                return(self.generation_data, binom_prob_list, mut_n_list, fig) 
            
            self.generation_data.append(offspring)
            
        # Plot the number of mutants over time
        # Count the number of individuals with a value of 1 in each generation
        num_mutants = [np.count_nonzero(generation == 1) for generation in self.generation_data]
        
        # Plot the number of mutants over time
        fig, ax = plt.subplots()
        ax.plot(range(len(num_mutants)), np.log(num_mutants))
        ax.set_xlabel("Time in Generations")
        ax.set_ylabel("Number of mutants ln(N)")
        ax.set_title(f"Mutant allele frequency over time (s={self.s})")
        plt.show()
        
        return(self.generation_data, binom_prob_list, mut_n_list, fig) 
    

def build_leaf_to_root_connections(tree_mask, mut_samples):
    """
    Simulate the population using the Wright-Fisher model with selection.
    """
    # List of each generation, where each generation is a dict from node_idx to set of leaves
    node_to_leaves = []
    desired_number = mut_samples

    # Initialize the list of generations
    for generation_idx, generation in enumerate(tree_mask):
        node_to_leaves.append({})
        for node_idx, node in enumerate(generation):
            if node == 1:
                node_to_leaves[-1][node_idx] = set()

    assert desired_number <= len(node_to_leaves[-1])
    node_to_leaves[-1] = {
        id: set() 
        for id in np.random.choice(
            list(node_to_leaves[-1].keys()), 
            desired_number, 
            replace=False
        )
    }

    # Go backward from leaf to root, randomly assigning each leaf to a parent
    for leaf_idx in node_to_leaves[-1].keys():
        leaf_to_follow = None
        for generation_idx, generation in reversed(
            list(enumerate(node_to_leaves[:-1]))
        ):
            # If we already have a leaf to follow, pick it's ancestor, otherwise pick a random one
            if leaf_to_follow is not None:
                parent_idx = next(
                    node
                    for node, leaves in generation.items()
                    if leaf_to_follow in leaves
                )
            else:
                parent_idx = random.choice(list(generation.keys()))

            # If the parent already has a leaf, pick it's ancestor in all previous generations to avoid cycles
            if len(node_to_leaves[generation_idx][parent_idx]) > 0:
                leaf_to_follow = list(node_to_leaves[generation_idx][parent_idx])[0]

            # Add the leaf to the parent
            node_to_leaves[generation_idx][parent_idx].add(
                (len(node_to_leaves) - 1, leaf_idx)
            )

    # Drop any non leaves that weren't connected
    # Create a dict from node coordinates to leaf coordinates
    result = {}
    for generation_idx, generation in enumerate(node_to_leaves):
        for node_idx, leaves in generation.items():
            if generation_idx == len(node_to_leaves) - 1 or len(leaves) > 0:
                result[str((generation_idx, node_idx))] = {str(leaf) for leaf in leaves}

    return result

def clusters_to_nodes(tree_clusters):
    """
    the following function goes from a nested dictionary of clades, creates the root node, adds node to that root 
    and then establishes parent and children relationships    """
    
    label_to_node = {cluster: Node(label=cluster) for cluster in tree_clusters}
    root_name = list(tree_clusters.keys())[0]
    
    # Connect each node with it's leaves 
    for parent_label, leaf_labels in tree_clusters.items():
        parent_node = label_to_node[parent_label]
        for leaf_label in leaf_labels:
            leaf_node = label_to_node[leaf_label]
            parent_node.add_child(leaf_node)

    # Connect each node with all descendants
    for node1 in label_to_node.values():
        for node2 in label_to_node.values():
            if node1 == node2:
                continue

            node_1_gen = int(node1.label.split(',')[0][1:]) # split the key on comma, take the first element, and remove the opening parenthesis
            node_2_gen = int(node2.label.split(',')[0][1:]) # split the key on comma, take the first element, and remove the opening parenthesis
            
            if abs(node_2_gen - node_1_gen) > 1 :
                continue
            
            possible_parent = node1 if len(node1.child_nodes()) > len(node2.child_nodes()) else node2
            possible_child = node2 if possible_parent is node1 else node1
            
            parent_leaf_labels = {node.label for node in possible_parent.child_nodes()}
            child_leaf_labels = {node.label for node in possible_child.child_nodes()}
            is_descendant = child_leaf_labels.issubset(parent_leaf_labels)

            if is_descendant and possible_child not in possible_parent.child_nodes() and possible_child.child_nodes():
                possible_parent.add_child(possible_child)

########  Remove non-direct descendants
    next_nodes = {label_to_node[root_name]}
    while next_nodes:
        node = next_nodes.pop()
        for child in node.child_nodes():
            next_nodes.add(child)
            child.set_parent(node)
            for grandchild in child.child_nodes():
                if grandchild in node.child_nodes():
                    node.remove_child(grandchild)

    # Create the tree using TreeSwift
    tree = Tree()
    tree.root = label_to_node[root_name]

    return tree

def assign_edge_lengths(mu, tree):
    """
    Iterate through the tree class and assign edge lengths based on a Poisson distribution with mean rate μ.
    """
    for node in tree.traverse_preorder():
        length = np.random.poisson(mu)
        node.set_edge_length(length)
    return tree

def normalise_tree_lengths(tree):
    current_gen = {tree.root}
    while current_gen:
        #gen_lengths = {node.get_edge_length() for node in current_gen}
        #average_length = sum(gen_lengths) / len(gen_lengths)
        next_gen = set()
        for node in current_gen:
            children = node.child_nodes()
            gen_lengths = {cnode.get_edge_length() for cnode in children}
            average_length = sum(gen_lengths) / len(gen_lengths)
            cnode.set_edge_length(average_length)
            next_gen.update(node.child_nodes())
        current_gen = next_gen


def read_observed_data(observed_data_path):
    """
    Read observed tree and calculate LTT statistics and return or read in LTT statistics straight
    """

    # Define the path to the file containing the tree
    tree_file = observed_data_path

    # Load the tree from the TSV file
    with open(tree_file) as f:
        tree_str = f.read()
    tree = treeswift.read_tree_newick(tree_str)
    
    # Calculate lineage through time plot statistics
    ltt = tree.lineages_through_time()

    # Save the results in a data structure
    results = {
        "tree": tree,
        "ltt": ltt,
    }
    print("Tree:", results["tree"])
    print("LTT statistics:", results["ltt"])
    return tree,ltt


##### ------------- Wright-Fisher Simulation ------------------------------ #

def simulate_population_and_tree(N, generations, mut_samples, s, mu):
    # initiate population
    popul = Population(N, generations, s) 
    # go from population array to tree_clusters dictionary
    gen, _prob, _mut, _fig = popul.simulate_population()
    # create genealogy and save in tree_clusters
    tree_clusters = build_leaf_to_root_connections(gen, mut_samples)
    # create phylo tree
    gen_tree = clusters_to_nodes(tree_clusters)
    from treeswift import read_tree_newick
    tree_string = gen_tree.newick()
    # read newick tree
    phy_tree = read_tree_newick(tree_string)
    # assign random edge (branch) lengths
    phy_tree_mut = assign_edge_lengths(mu, phy_tree)
    # visualise tree
    import matplotlib.patches as patches
    white_patch = patches.Patch(color='black', label=f"Phylogenetic tree (s={s})")
    plot = phy_tree_mut.draw(show_labels=False, handles=[white_patch])
    #normalise_tree_lengths(phy_tree_mut)
    # calculate ltt stats and plot using treeswift
    ltt_gen_tree = phy_tree_mut.lineages_through_time()
    # write tree to newick txt file
    return(plot,ltt_gen_tree)
    #gen_tree_expanded.write_tree_newick("output_gen_tree.tree.nwk", hide_rooted_prefix=True)

# initialize an empty list to store the results
results = []

def run_simulation_with_restart():
    for s in s_values:
        num_retries = 0
        while num_retries <= 2:
            try:
                result = simulate_population_and_tree(N=args.N, generations=args.generations, mut_samples=args.mut_samples, s=args.s, mu=args.mu)
                results.append(result)
            except AssertionError:
                num_retries += 1
                print("AssertionError occurred, restarting simulation...")


# call the function with the command-line arguments
#result = simulate_population_and_tree(N=args.N, generations=args.generations, mut_samples=args.mut_samples, s=args.s, mu=args.mu)


##### ------------- Approximate Bayesian Criterion ----------------- #
"""
Set Priors and run ABC analysis using abcpy
"""

##### Distance function between observed and synthetic data
def distance_function():
    """
    Calculate the distance between the observed and synthetic data.
    """
    pass



def estimate_parameters(epsilon):
    # set priors
    s_prior = abcpy.Distribution(stats.uniform, 0, 2) # selection
    
    pass
        

    