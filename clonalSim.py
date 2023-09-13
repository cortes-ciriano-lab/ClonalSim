import matplotlib.pyplot as plt
import numpy as np
import random
import treeswift
from treeswift import Tree, Node
import argparse
import csv
from UltrametricConversion import traverse_and_run_average
from UltrametricConversion import transform_data
#from UltrametricConversion import normalise_data
from UltrametricConversion import handle_labelless_trees
from UltrametricConversion import handle_none_edge_trees
from typing import Dict, Iterable, List
import os

# create an argparse parser
parser = argparse.ArgumentParser(description="Simulate population and tree")

# add arguments for the simulation parameters
parser.add_argument("--N", type=int, help="population size")
parser.add_argument("--generations", type=int, required=False, help="number of generations to simulate")
parser.add_argument("--disease", type=int, help="number of generations the disease starts")
parser.add_argument("--mut_samples", type=int, required=False, help="number of mutation samples")
parser.add_argument("--s", type=float, help="selection coefficient")
parser.add_argument("--mu", type=float, help="mutation rate")
parser.add_argument("--epsilon", type=float, help="Epsilon Threshold")
parser.add_argument('--output_path', type=str, default='.', help='output path')
parser.add_argument('--observed_data_path', type=str, default='.', help='observed_data_path', required=False)

# parse the command-line arguments
args = parser.parse_args()


class Population:
    def __init__(self, N, generations, disease, s):
        self.N = N
        self.s = s
        self.generations = generations
        self.disease = disease
        self.generation_data = []

    def __str__(self):
        return f"Population size: {self.N}, Generations: {self.generations}, Disease_Onset: {self.disease}, Selection: {self.s}"

    def simulate_population(self):
        """
        Simulate the population using the Wright-Fisher model with selection.
        """

        # Initialize the first population
        population = np.zeros(self.N)
        self.generation_data.append(population)
        binom_prob_list = []
        mut_n_list = []

        for gen in range(1, self.generations + 1):
            if gen - 1 < len(self.generation_data) and len(self.generation_data) > 0:
                if gen < self.disease:
                    self.generation_data.append(np.zeros(self.N))
                    mut_n_list.append(0)
                    binom_prob_list.append(0)
                elif gen == self.disease:
                    new_population = np.copy(population)
                    new_population[random.randint(0, self.N - 1)] = 1
                    self.generation_data.append(new_population)
                    num_mutants = 1
                    mut_n_list.append(num_mutants)
                elif gen > self.disease:
                    mut_n = mut_n_list[-1]
                    cancer_p = ((1 + self.s) * mut_n) / (self.N + (mut_n * self.s))
                    if cancer_p > 1:
                        break
                    binom_prob_list.append(cancer_p)

                    offspring_l = np.random.binomial(size=1, n=self.N, p=cancer_p)
                    offspring_n = offspring_l.item()
                    random_indices = random.sample(range(self.N), k=offspring_n)
                    offspring = np.zeros(self.N)
                    offspring[random_indices] = 1

                    num_mutants = np.sum(offspring)
                    mut_n_list.append(num_mutants)

                    self.generation_data.append(offspring)


        return (self.generation_data, binom_prob_list, mut_n_list)


def build_leaf_to_root_connections(
        tree_mask: Iterable[Iterable[int]],
        mut_samples: int,
) -> Dict[str, List[str]]:
    """
    Simulate the population using the Wright-Fisher model with selection. The simulated
    geneologies allows only two descendants per cell.

    tree_mask: which is a list of lists, representing a matrix where each row is
    a generation and a value of 1 represents a mutated sample and 0 is a normal sample.

    Returns:
        A dictionary of num_samples mutated cells to their geneology (path to root)
        Cells are named with the convention <generation>-<index>
    """

    # Assertion that the last generation has a higher number of mut_cells than the desired_number of tips for the sim_tree
    assert sum(tree_mask[-1]) >= mut_samples

    max_children = 20

    # Convert the tree mask into list of generations which are sets of node names
    generation_to_nodes: List[Set[str]] = [
        {f"{gen_idx}-{node_idx}" for node_idx, node in enumerate(gen) if node == 1}
        for gen_idx, gen in enumerate(tree_mask)
    ]

    # Construct mappings of nodes to their children and parents
    node_to_children: Dict[str, Set[str]] = {}
    node_to_parent: Dict[str, str] = {}
    for generation in generation_to_nodes:
        for node in generation:
            node_to_children[node] = set()
            node_to_parent[node] = None

    # Find the root node (not necessarily in the first generation if there are empty generations)
    last_gen = None
    for generation_idx, generation in reversed(list(enumerate(generation_to_nodes))):
        if len(generation) == 0:
            if last_gen is not None:
                break
            else:
                continue
        last_gen = generation
        last_gen_idx = generation_idx
    assert len(last_gen) == 1
    root = next(iter(last_gen))
    root_gen = last_gen_idx

    # Simulate num_samples geneologies from root to a leaf
    leaves_to_path_to_root: Dict[str, List[str]] = {}
    while len(leaves_to_path_to_root) < mut_samples:

        path = [root]
        for generation in generation_to_nodes[root_gen + 1:]:
            node = path[0]

            # Possible children are nodes without parents, or nodes which
            # are already children of the current node
            possible_children = [n for n in generation if node_to_parent.get(n) in [None, node]]

            # If the node already has the maximum number of children, we
            # can't select a new one
            if len(node_to_children[node]) == max_children:
                possible_children = list(node_to_children[node])

            # Filter possible children to ones which don't yet have max children
            possible_children = [
                n for n in possible_children if len(node_to_children[n]) < max_children
            ]

            # If there are no possible children, this geneology was not possible
            if len(possible_children) == 0:
                break

            # Choose a random child and add it to the path to the root
            path.insert(0, random.choice(possible_children))  # insert adds child at index 0

        # If the path is unique, connect all parents and children and add the path
        if len(path) == len(generation_to_nodes) - root_gen:
            leaf = path.pop(0)
            if leaf not in leaves_to_path_to_root:
                child = leaf
                for node in path:
                    assert node_to_parent.get(child) in [None, node]
                    node_to_children[node].add(child)
                    node_to_parent[child] = node
                    child = node
                leaves_to_path_to_root[leaf] = path

    return leaves_to_path_to_root


def clusters_to_nodes(tree_clusters: Dict[Node, Iterable[Node]]) -> Tree:
    """
    Input: Dictionary from Leaf nodes to a list of nodes from that leaf to the root
    Output: Tree Object, all nodes with immediate parent and children
    """
    label_to_node = {leaf: Node(label=leaf) for leaf in tree_clusters}

    for leaf_label, parent_labels in tree_clusters.items():
        leaf_node = label_to_node[leaf_label]
        prev_parent = leaf_node
        # for parent_label in reversed(parent_labels):
        for parent_label in parent_labels:
            parent_node = label_to_node.get(parent_label)
            if parent_node is None:
                parent_node = Node(label=parent_label)
                label_to_node[parent_label] = parent_node

            # Connect the new parent to its child (Previous parent)
            if prev_parent not in parent_node.child_nodes():
                parent_node.add_child(prev_parent)
            prev_parent = parent_node

    # Select the element from label_to_node based on label starting from 0
    selected_node_label = None
    first_leaf, path_for_first_leaf = next(iter(tree_clusters.items()))
    selected_node_label = path_for_first_leaf[-1]  # this will give you the last node in the path for the first leaf
    print(f"This the root: {selected_node_label}")
    selected_node = label_to_node[selected_node_label]

    # Create the tree using TreeSwift
    tree = Tree()
    tree.root = selected_node

    return tree


def assign_edge_lengths(mu, tree, disease_onset):
    """
    Iterate through the tree class and assign edge lengths based on a Poisson distribution with mean rate μ.
    """

    one_cell_gens = disease_onset

    first_node = True  # to identify the root node

    for node in tree.traverse_preorder():
        if first_node:
            length = np.random.poisson(mu) * one_cell_gens
            first_node = False
            print(f'The root length is:{length}')
        else:
            length = np.random.poisson(mu)
            print(node.label)
            print(f'The node length is:{length}')
        node.set_edge_length(length)
    return tree


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
    handle_labelless_trees(tree)
    handle_none_edge_trees(tree)
    normalised_obs_tree = traverse_and_run_average(tree)
    # Calculate lineage through time plot statistics
    ltt = normalised_obs_tree.lineages_through_time(
        show_plot=False)  # , export_filename=f"{output_path}/Plot_obs_ltt_ultrametric_(s={s}).png")
    # ltt = normalised_obs_tree.lineages_through_time(show_plot=False, export_filename=f"{output_path}/Plot_obs_ltt_ultrametric_(s={s}).png")
    list_of_tuples_obs = [(key, value) for key, value in ltt.items()]
    data_transformed_obs = transform_data(list_of_tuples_obs)
    obs_tree_length = data_transformed_obs[-1][0]
    #norm_ltt = normalise_data(data_transformed_obs)

    return tree, ltt, obs_tree_length


from scipy.integrate import trapz

def calculate_epsilon(norm_data1, norm_data2, mut_samples, obs_tree_length):

    curve1 = norm_data1
    curve2 = norm_data2

    # Extract x and y values for each curve
    x_curve1, y_curve1 = zip(*curve1)
    x_curve2, y_curve2 = zip(*curve2)

    # Create a common range of x-values
    x_common = np.linspace(min(min(x_curve1), min(x_curve2)), max(max(x_curve1), max(x_curve2)), 1000)

    # Interpolate y-values for each curve at the common x-values
    y_interp_curve1 = np.interp(x_common, x_curve1, y_curve1)
    y_interp_curve2 = np.interp(x_common, x_curve2, y_curve2)

    # Find intervals where one curve is above the other
    fill_x = []
    fill_y1 = []
    fill_y2 = []

    for i in range(len(x_common) - 1):
        if y_interp_curve1[i] >= y_interp_curve2[i]:
            fill_x.extend([x_common[i], x_common[i + 1]])
            fill_y1.extend([y_interp_curve1[i], y_interp_curve1[i + 1]])
            fill_y2.extend([y_interp_curve2[i], y_interp_curve2[i + 1]])
        else:
            fill_x.extend([x_common[i], x_common[i + 1]])
            fill_y1.extend([y_interp_curve2[i], y_interp_curve2[i + 1]])
            fill_y2.extend([y_interp_curve1[i], y_interp_curve1[i + 1]])

    # Calculate the area between the curves using the trapezoidal rule
    area_between_curves_full = trapz(np.abs(np.array(fill_y1) - np.array(fill_y2)), fill_x)
    area_between_curves = area_between_curves_full / (mut_samples * obs_tree_length)

    # Plotting
    fig = plt.figure(figsize=(10, 6))
    plt.plot(x_curve1, y_curve1, label='MPN tree 1')
    plt.plot(x_curve2, y_curve2, label='MPN tree 1')
    plt.fill_between(fill_x, fill_y1, fill_y2, where=np.array(fill_y1) >= np.array(fill_y2), color='blue', alpha=0.3)
    plt.fill_between(fill_x, fill_y1, fill_y2, where=np.array(fill_y1) < np.array(fill_y2), color='red', alpha=0.3)
    plt.xlabel('Scaled Time')
    plt.ylabel('Lineages')
    plt.title('LTT Curves and Area Between Them')
    plt.legend()
    plt.grid(True)
    plt.ylim(0)  # Set the lower limit of y-axis to 0
    # Add the area between curves value to the plot
    plt.text(0, 0, f'Area: {area_between_curves}', fontsize=12)

    # Print the result
    print("The area between the curves is:", area_between_curves)

    return fig, area_between_curves

def euclidean_distance_dicts(dict1, dict2):

    # Reverse the dictionaries
    reversed_dict1 = {v: k for k, v in dict1.items()}
    reversed_dict2 = {v: k for k, v in dict2.items()}

    # Ensure that the reversed dictionaries have the same keys
    if set(reversed_dict1.keys()) != set(reversed_dict2.keys()):
        raise ValueError("The dictionaries have different values.")
    
    # Combine keys from both dictionaries
    all_keys = set(dict1.keys()) | set(dict2.keys())

    distance_squared = 0
    for key in all_keys:
        value1 = dict1.get(key, 0)
        value2 = dict2.get(key, 0)
        distance_squared += (value2 - value1) ** 2

    return distance_squared ** 0.5


def euclidean_distance_dicts(dict1, dict2):
    # Reverse the dictionaries
    reversed_dict1 = {v: k for k, v in dict1.items()}
    reversed_dict2 = {v: k for k, v in dict2.items()}

    # Ensure that the reversed dictionaries have the same keys
    if set(reversed_dict1.keys()) != set(reversed_dict2.keys()):
        raise ValueError("The dictionaries have different values.")

    # Combine keys from both dictionaries
    all_keys = set(dict1.keys()) | set(dict2.keys())

    distance_squared = 0
    for key in all_keys:
        value1 = dict1.get(key, 0)
        value2 = dict2.get(key, 0)
        distance_squared += (value2 - value1) ** 2

    return distance_squared ** 0.5


##### ------------- Wright-Fisher Simulation ------------------------------ ##########

def simulate_population_and_tree(N, generations, disease, mut_samples, s, mu, output_path, observed_d_path, epsilon):
    print("Simulating population...")
    # initiate population
    popul = Population(N, generations, disease, s)
    # go from population array to tree_clusters dictionary
    gen, prob, mut = popul.simulate_population()
    # fig.savefig(f"{output_path}/Simulation_{num_retries}_with_mutants_in_time_(s={s}).png")
    print("Population Done...")
    # create genealogy and save in tree_clusters
    print("Simulating Genealogy...")
    tree_clusters = build_leaf_to_root_connections(gen, mut_samples)
    # create phylo tree
    gen_tree = clusters_to_nodes(tree_clusters)
    from treeswift import read_tree_newick
    tree_string = gen_tree.newick()
    phy_tree = read_tree_newick(tree_string)
    print("Genealogy Done")
    # assign random edge (branch) lengths
    phy_tree_mut = assign_edge_lengths(mu, phy_tree, disease)

    # make tree ultrametric
    handle_labelless_trees(phy_tree_mut)
    handle_none_edge_trees(phy_tree_mut)
    normalised_tree = traverse_and_run_average(phy_tree_mut)
    print("Ultrametric tree done")

    # calculate ltt stats and plot using treeswift
    ltt_gen_tree = phy_tree_mut.lineages_through_time(show_plot=False)# export_filename=f"{output_path}/Plot_ltt_ultrametric_(s={s}).png")
    # normalise ltt stats
    list_of_tuples_tree = [(key, value) for key, value in ltt_gen_tree.items()]
    data_transformed = transform_data(list_of_tuples_tree)
    #norm_data = normalise_data(data_transformed)
    print("LTT Statistics Done")

    print("Reading Observed Data and Calculating LTT...")
    obs_tree, obs_ltt, obs_tree_length = read_observed_data(observed_d_path)
    fig_abc, abc = calculate_epsilon(obs_ltt, data_transformed, mut_samples, obs_tree_length)
    
    print("Reading Observed Data and Calculating LTT...")
    obs_tree, obs_ltt, obs_tree_length = read_observed_data(observed_d_path)

    fig_abc, abc = calculate_epsilon(obs_ltt, data_transformed, mut_samples, obs_tree_length)

    if abc is not None and epsilon is not None and abc < epsilon:
        fig_abc.savefig(f"{output_path}/Simulation_{N}_{disease}_with_abc_fig_(s={s}).png")
    print("Area Under the Curve calculated")

    eud = euclidean_distance_dicts(obs_ltt, ltt_gen_tree)

    return phy_tree_mut, abc, eud



max_retries = 1000
retry_count = 0

while retry_count < max_retries:
    try:
        result_tree, abc_epsilon, eud = simulate_population_and_tree(N=args.N, generations=args.generations,
                                                                disease=args.disease, mut_samples=args.mut_samples,
                                                                s=args.s, mu=args.mu, output_path=args.output_path,
                                                                observed_d_path=args.observed_data_path,
                                                                epsilon=args.epsilon)
        print(f"abc_epsilon: {abc_epsilon}")  # Debugging line

        # If abc_epsilon is less than a value, then create the file, write the header and the results
        if abc_epsilon < args.epsilon:

            # save simulated tree
            result_tree.write_tree_newick(
                f"{args.output_path}/Simulation_{args.N}_{args.generations}_{args.disease}_{args.mut_samples}_{args.s}_output_gen_tree.nwk",
                hide_rooted_prefix=False)

            file_path = f"{args.output_path}/Simulation_results_{args.N}_{args.generations}_{args.disease}_{args.mut_samples}_{args.s}.tsv"
            # Check if file already exists (i.e., has been written to in a previous run)
            header_needed = not os.path.exists(file_path)

            with open(file_path, "a", newline='') as f:
                if header_needed:
                    # Write the header with variable names
                    f.write(
                        "ABC_Epsilon\tEUD\tN\tGenerations\tDisease\tMut_Samples\tS\tMu\tOutput_Path\tObserved_Data_Path\n")

                f.write(
                    f"{abc_epsilon}\t{eud}\t{args.N}\t{args.generations}\t{args.disease}\t{args.mut_samples}\t{args.s}\t{args.mu}\t{args.output_path}\t{args.observed_data_path}\n")
            break
        else:
            raise ValueError("ABC_Epsilon is greater than or equal to given epsilon")
    except (AssertionError, ValueError) as error:
        print(f"Error occurred: {error}, restarting simulation...")
        retry_count += 1