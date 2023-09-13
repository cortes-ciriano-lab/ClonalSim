# Ultrametric conversion code #

from treeswift import Tree

def handle_labelless_trees(tree_obj):
    for idx, node in enumerate(tree_obj.traverse_postorder()):
        if not node.label:
            node.label = f"node-{idx}"

def handle_none_edge_trees(tree_obj):
    for node in tree_obj.traverse_postorder():
        if not node.get_edge_length():
            node.set_edge_length(0)

def traverse_and_run_average(tree_obj):
    # Construct a mapping from labels to nodes to modify the edge lengths by reference
    label_to_node = {node.label: node for node in tree_obj.traverse_preorder()}

    # Make a list of all subtrees, in ascending order of the number of nodes in the subtree
    subtrees = sorted(tree_obj.traverse_postorder(), key=lambda node: len(list(node.traverse_postorder())))

    # Make a mapping between subtrees and the length to the leaves
    node_to_subtree_length = {node.label: None for node in tree_obj.traverse_postorder()}

    # Iterate through all subtrees in ascending order
    for subtree in subtrees:
        if not subtree.child_nodes():
            node_to_subtree_length[subtree.label] = 0
        else:

            # Calculate the mean length of all branches from the root
            mean_branch_length = sum(
                child.get_edge_length() + node_to_subtree_length[child.label] 
                for child in subtree.child_nodes()
            ) / len(subtree.child_nodes())
            node_to_subtree_length[subtree.label] = mean_branch_length

            # Adjust each branch so that the subtree lengths match
            for child_subtree in subtree.child_nodes():
                new_subtree_length = node_to_subtree_length[child_subtree.label] + child_subtree.get_edge_length()
                if new_subtree_length != 0:
                    for child_subtree_node in child_subtree.traverse_postorder():
                        label_to_node[child_subtree_node.label].set_edge_length(
                            child_subtree_node.get_edge_length() / new_subtree_length * mean_branch_length)

    return tree_obj


def transform_data(data):
    data_transformed = []
    
    for i in range(len(data) - 1):
        current_x, current_y = data[i]
        next_x, _ = data[i+1]
        
        for x in range(int(current_x), int(next_x)):
            data_transformed.append((float(x), current_y))
    
    last_x, last_y = data[-1]
    # Assuming that we need to add one more point after the last x based on the example provided
    data_transformed.append((float(last_x), last_y))
    data_transformed.append((float(last_x) + 1, last_y))
    
    return data_transformed


# def normalise_data(data):
#     data_norm = []
    
#     for i in range(len(data)):
#         current_x, current_y = data[i]
        
#         # Normalise the 0th element of the tuple with the last 0th element of the data
#         current_x = current_x / data[-1][0]
#         current_y = current_y / data[-1][1]
        
#         data_norm.append((current_x, current_y))
    
#     return data_norm