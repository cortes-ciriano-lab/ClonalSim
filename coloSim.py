import numpy as np
import random
import abcpy
import scipy.stats as stats
import treeswift
from treeswift import Tree, Node



####################################################
## Simulate a Wright-Fisher Model with selection ###
####################################################

class Population:
    def __init__(self, N, generations, mut_samples, s):
        self.N = N
        self.s = s
        self.generations = generations
        self.mut_samples = mut_samples
        self.generation_data = []

    def __str__(self):
        return f"Population size: {self.N}, Generations: {self.generations}, Mutant Samples: {self.mut_samples}"

    def simulate_population(self):
        """
        Simulate the population using the Wright-Fisher model with selection.
        """

        # Initialize the first population
        population = np.zeros(self.N, dtype=int)
        population[random.randint(0, self.N - 1)] = 1
        self.generation_data.append(population)
        #cancer_p_list  = []
        #offspring_list = []
        
        for gen in range(self.generations):
            mut_n = len(np.where(self.generation_data[gen] == 1)[0])
            cancer_p = (1 + self.s) * mut_n / (self.N + mut_n * self.s)
            offspring = np.random.binomial(n=1, p=cancer_p, size=self.N)
            self.generation_data.append(offspring)
            #cancer_p_list.append(cancer_p)
            #offspring_list.append(offspring)
        #return(self.generation_data)
        # Plot the number of mutants over time
        
        # Count the number of individuals with a value of 1 in each generation
        num_mutants = [np.count_nonzero(generation == 1) for generation in self.generation_data]
        
        # Plot the number of mutants over time
        fig, ax = plt.subplots()
        ax.plot(range(len(num_mutants)), num_mutants)
        ax.set_xlabel("Generation")
        ax.set_ylabel("Number of mutants")
        ax.set_title("Mutant allele frequency over time")
        plt.show()
        
        return(self.generation_data, fig)
            



def build_leaf_to_root_connections(tree_mask):
    # List of each generation, where each generation is a dict from node_idx to set of leaves
    node_to_leaves = []

    # Initialize the list of generations
    for generation_idx, generation in enumerate(tree_mask):
        node_to_leaves.append({})
        for node_idx, node in enumerate(generation):
            if node == 1:
                node_to_leaves[-1][node_idx] = set()

    # Go backward from leaf to root, randomly assigning each leaf to a parent
    for leaf_idx in node_to_leaves[-1].keys():
        for generation_idx, generation in reversed(
            list(enumerate(node_to_leaves[:-1]))
        ):
            parent_idx = random.choice(list(generation.keys()))
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


if _name_ == "_main_":
    build_leaf_to_root_connections(tree_mask)


def clusters_to_nodes(tree_clusters):
    """
    the following function goes from a nested dictionary of clades, creates the root node, adds node to that root 
    and then establishes parent and children relationships    """
    
    label_to_node = {cluster: Node(label=cluster) for cluster in tree_clusters}

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

            possible_parent = node1 if len(node1.child_nodes()) > len(node2.child_nodes()) else node2
            possible_child = node2 if possible_parent is node1 else node1
            
            parent_leaf_labels = {node.label for node in possible_parent.child_nodes()}
            child_leaf_labels = {node.label for node in possible_child.child_nodes()}
            is_descendant = child_leaf_labels.issubset(parent_leaf_labels)

            if is_descendant and possible_child not in possible_parent.child_nodes() and possible_child.child_nodes():
                possible_parent.add_child(possible_child)

   # Remove non-direct descendants
    next_nodes = {label_to_node["1"]}
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
    tree.root = label_to_node["1"]

    return tree


def assign_edge_lengths(self, mu, tree):
    """
    Iterate through the tree class and assign edge lengths based on a Poisson distribution with mean rate μ.
    """
    for node in tree.nodes():
        length = np.random.poisson(mu)
        node.set_edge_length(length)
    return tree

def read_observed_data(observed):
    """
    Read observed tree and calculate LTT statistics and return or read in LTT statistics straight
    """

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


# ------------- Run Simulation ------------- #

# to add: for loop for multiple iterations & all simulations in list

def simulate_population_and_tree(N, generations, mut_samples, s, mu):
    # initiate population
    popul = Population(N, generations, mut_samples, s) 
    # initialise genealogy
    Gen = Genealogia() 
    # connect all random mutated cells
    Gen.connect_random_cells(popul) 
    # from matrix go to tree_clusters dictionary
    tree_clusters=genealogy_to_cluster(Gen.genealogy) 
    # create phylo tree
    gen_tree = clusters_to_nodes(tree_clusters)
    # assign random edge (branch) lengths
    gen_tree_expanded = assign_edge_lengths(gen_tree, mu)
    # calculate ltt stats and plot using treeswift
    ltt_gen_tree = gen_tree_expanded.lineages_through_time()
    # write tree to newick txt file
    gen_tree.write_tree_newick("output_gen_tree.tree.nwk", hide_rooted_prefix=True)



# run ABC 

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
        