import matplotlib.pyplot as plt
import numpy as np
import random
import treeswift
from treeswift import Tree, Node
import argparse
import csv
from UltrametricConversion import traverse_and_run_average
from UltrametricConversion import transform_data
from UltrametricConversion import remove_root_node
from UltrametricConversion import normalize_ltt_data
from UltrametricConversion import handle_labelless_trees
from UltrametricConversion import handle_none_edge_trees
from typing import Dict, Iterable, List
import os
import time
from memory_profiler import profile, memory_usage

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
parser.add_argument('--model_type', type=str, default='WF', help='WF OR MSI', required=False)
parser.add_argument('--MSI_onset', type=str, default='WF', help='when MSI starts', required=False)
parser.add_argument('--normal_rate_CRC', type=str, default='mutation rate for normal', help='WF OR MSI', required=False)

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

@profile
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
    # if sum(tree_mask[-1]) >= mut_samples:
    #     print("Not enough tips")
    #     print(sum(tree_mask[-1]))

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

@profile
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
    #print(f"This the root: {selected_node_label}")
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
            #print(f'The root length is:{length}')
        else:
            length = np.random.poisson(mu)
            # print(node.label)
            # print(f'The node length is:{length}')
        node.set_edge_length(length)
    return tree


def assign_edge_lengths_MSI(mu, tree, disease_onset, MSI_onset, normal_CRC_rate):
    """
    Iterate through the tree class and assign edge lengths based on a Poisson distribution with mean rate μ.
    """

    one_cell_gens = disease_onset

    first_node = True  # to identify the root node

    for node in tree.traverse_preorder():
        if first_node:
            # For MSI_onset times, do np.random.poisson(normal_CRC_rate) and save the results
            normal_lengths = [np.random.poisson(float(normal_CRC_rate)) for _ in range(int(MSI_onset))]

            # For (one_cell_gens - MSI_onset) times, do np.random.poisson(mu) and save those results
            msi_lengths = [np.random.poisson(mu) for _ in range(int(one_cell_gens) - int(MSI_onset))]

            # Sum all results from these two vectors and save them as length
            length = sum(normal_lengths) + sum(msi_lengths)

            first_node = False
            # print(f'The root length is:{length}')
        else:
            length = np.random.poisson(mu)
            # print(node.label)
            # print(f'The node length is:{length}')
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
    phy_tree_mut = remove_root_node(tree)
    normalised_obs_tree = traverse_and_run_average(phy_tree_mut)
    # Calculate lineage through time plot statistics
    ltt = normalised_obs_tree.lineages_through_time(
        show_plot=False)  # , export_filename=f"{output_path}/Plot_obs_ltt_ultrametric_(s={s}).png")
    # ltt = normalised_obs_tree.lineages_through_time(show_plot=False, export_filename=f"{output_path}/Plot_obs_ltt_ultrametric_(s={s}).png")
    list_of_tuples_obs = [(key, value) for key, value in ltt.items()]
    data_transformed_obs = transform_data(list_of_tuples_obs)
    obs_tree_length = data_transformed_obs[-1][0]

    return normalised_obs_tree, ltt, obs_tree_length



    # Calculate the area between the curves using the trapezoidal rule
    area_between_curves_full = trapz(np.abs(np.array(fill_y1) - np.array(fill_y2)), fill_x)
    area_between_curves = area_between_curves_full / (mut_samples * obs_tree_length)

    # Scale the area between curves
    scaled_area_between_curves = area_between_curves * scale_factor

    # Plotting
    # fig = plt.figure(figsize=(10, 6))
    # plt.plot(x_curve1, y_curve1, label='Observed tree')
    # plt.plot(x_curve2, y_curve2, label='Simulated tree')
    # plt.fill_between(fill_x, fill_y1, fill_y2, where=np.array(fill_y1) >= np.array(fill_y2), color='blue', alpha=0.3)
    # plt.fill_between(fill_x, fill_y1, fill_y2, where=np.array(fill_y1) < np.array(fill_y2), color='red', alpha=0.3)
    # plt.xlabel('Scaled Time')
    # plt.ylabel('Lineages')
    # plt.title('LTT Curves and Area Between Them')
    # plt.legend()
    # plt.grid(True)
    # plt.ylim(0)  # Set the lower limit of y-axis to 0
    # # Add the area between curves value to the plot
    # plt.text(0, 0, f'Scaled Area: {scaled_area_between_curves}', fontsize=12)
    
    # Print the result
    print("The scaled area between the curves is:", scaled_area_between_curves)

    return scaled_area_between_curves


# NEW!!! 
def extract_trunk_path(tree: Tree) -> List[float]:
    node = tree.root
    trunk_lengths = []
    while len(node.child_nodes()) == 1:
        child = node.child_nodes()[0]
        trunk_lengths.append(child.edge_length if child.edge_length else 0)
        node = child
    return trunk_lengths

def trunk_distance(trunk1: List[float], trunk2: List[float]) -> float:
    min_len = min(len(trunk1), len(trunk2))
    trunk1 = trunk1[:min_len]
    trunk2 = trunk2[:min_len]
    return sum(abs(a - b) for a, b in zip(trunk1, trunk2))



##### ------------- Wright-Fisher Simulation ------------------------------ ##########

@profile
def simulate_population_and_tree_MSI(N, generations, disease, mut_samples, s, mu, output_path,
                                     MSI_onset, normal_CRC_rate):
    print("Simulating population...")
    popul = Population(N, generations, disease, s)
    gen, prob, mut = popul.simulate_population()
    print("Population simulated.")

    print("Simulating genealogy...")
    tree_clusters = build_leaf_to_root_connections(gen, mut_samples)
    gen_tree = clusters_to_nodes(tree_clusters)

    from treeswift import read_tree_newick
    tree_string = gen_tree.newick()
    phy_tree = read_tree_newick(tree_string)
    print("Genealogy tree created.")

    phy_tree_mut = assign_edge_lengths_MSI(mu, phy_tree, disease, MSI_onset, normal_CRC_rate)
    handle_labelless_trees(phy_tree_mut)
    handle_none_edge_trees(phy_tree_mut)
    phy_tree_mut = remove_root_node(phy_tree_mut)
    normalised_tree = traverse_and_run_average(phy_tree_mut)

    return normalised_tree

max_retries = 3
retry_count = 0
sim_counter = 0

# Load and process the observed tree trunk
obs_tree, _, _ = read_observed_data(args.observed_data_path)
obs_trunk = extract_trunk_path(obs_tree)

def memory_profile(func, *args, **kwargs):
    mem_usage, retval = memory_usage((func, args, kwargs), retval=True, interval=0.1, timeout=None, include_children=True)
    print('Memory usage (in increments of 0.1 seconds):', mem_usage)
    print('Maximum memory usage: {:.2f} MiB'.format(max(mem_usage)))
    return retval

while retry_count < max_retries:
    try:
        if args.model_type == "MSI":
            result_tree = memory_profile(
                simulate_population_and_tree_MSI,
                N=args.N,
                generations=args.generations,
                disease=args.disease,
                mut_samples=args.mut_samples,
                s=args.s,
                mu=args.mu,
                output_path=args.output_path,
                MSI_onset=args.MSI_onset,
                normal_CRC_rate=args.normal_rate_CRC
            )
        else:
            result_tree = memory_profile(
                simulate_population_and_tree,
                N=args.N,
                generations=args.generations,
                disease=args.disease,
                mut_samples=args.mut_samples,
                s=args.s,
                mu=args.mu,
                output_path=args.output_path
            )

        # Compare trunks
        sim_trunk = extract_trunk_path(result_tree)
        distance = trunk_distance(sim_trunk, obs_trunk)
        print(f"Trunk distance: {distance}")

        # Save tree and results
        timestamp = int(time.time())
        prefix = f"{args.output_path}/Simulation_{args.N}_{args.generations}_{args.disease}_{args.mut_samples}_{args.s}"
        if args.model_type == "MSI":
            prefix += f"_{args.MSI_onset}_{args.normal_rate_CRC}"

        tree_path = f"{prefix}_ID{timestamp}_tree.nwk"
        result_tree.write_tree_newick(tree_path, hide_rooted_prefix=False)

        file_path = f"{args.output_path}/AllSimulationResults.tsv"
        header_needed = not os.path.exists(file_path)

        with open(file_path, "a", newline='') as f:
            if header_needed:
                f.write("Trunk_Distance\tN\tGenerations\tDisease\tMut_Samples\tS\tMu\tModel_Type\tMSI_onset\tCRC_rate\tTree_Path\n")
            f.write(f"{distance}\t{args.N}\t{args.generations}\t{args.disease}\t{args.mut_samples}\t{args.s}\t{args.mu}\t{args.model_type}\t{args.MSI_onset}\t{args.normal_rate_CRC}\t{tree_path}\n")

        sim_counter += 1
        retry_count += 1

    except AssertionError as error:
        print(f"Error occurred: {error}, retrying...")
        retry_count += 1
