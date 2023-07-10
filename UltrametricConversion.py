# Ultrametric conversion code #

def get_path_to_root(node):
    path = []
    while node is not None:
        path.append(node.label)
        node = node.parent
    return path

def find_node_by_label(tree, label):
    for node in tree.traverse_postorder():
        if node.label == label:
            return node
    return None

def get_path_to_leaf(node):
    if node.is_leaf():
        return [[node.label]]
    paths = []
    for child in node.child_nodes():
        for path in get_path_to_leaf(child):
            paths.append([node.label] + path)
    return paths

def process_path(path, subtree):
    # dictionary to store edge lengths
    edge_lengths = {}
    # add weights dict
    edges_weights = {}

    for node_label in path:  
        node = find_node_by_label(subtree, str(node_label))
        edge_name = node_label
        # get the edge length
        edge_length = node.get_edge_length()
        # add the edge length to the dictionary
        edge_lengths[edge_name] = edge_length 
        print(f"Edge: {edge_name}, Length: {edge_length}")

    total_path_value = sum(edge_lengths.values())
    print(f"Total: {total_path_value}")

    for edge, length in edge_lengths.items():
        if total_path_value == 0 and length == 0:
            continue
        edge_weight = length / total_path_value
        edges_weights[edge] = edge_weight

    for edge, weight in edges_weights.items():
        print(f"Edge: {edge}, Weight: {weight}")
    print("\n")
    return edge_lengths, edges_weights, total_path_value


def upd_tree_lengths(obs_tree_obj, left_path, right_path, subtree_obj):
    # calculate distances on each side and weights
    re, rw, rt = process_path(path=right_path, subtree=subtree_obj)
    le, lw, lt = process_path(path=left_path, subtree=subtree_obj)
    
    # take the average branch length
    avg_branch_l = (rt + lt) / 2

    # dictionary to store edge lengths
    upd_edge_lengths = {}
    
    print("The right path updated edge lengths are:")
    for edge, weight in rw.items():
        upd_edge_lengths[edge] = weight * avg_branch_l
        print(f"Edge: {edge}, Updated length: {upd_edge_lengths[edge]}")

    print("The left path updated edge lengths are:")
    for edge, weight in lw.items():
        upd_edge_lengths[edge] = weight * avg_branch_l
        print(f"Edge: {edge}, Updated length: {upd_edge_lengths[edge]}")
        
    # now update the edge lengths on the tree
    for node in obs_tree_obj.traverse_postorder():
        if node.get_label() in upd_edge_lengths:
            val= upd_edge_lengths[node.get_label()]
            node.set_edge_length(length=val)

    print("Final averaged tree")
    obs_tree_obj.draw(show_labels=True)
    return obs_tree_obj

# def upd_tree_lengths(obs_tree_obj, paths, subtree_obj):
#     # dictionary to store edge lengths
#     upd_edge_lengths = {}

#     # list to store total path values
#     total_path_values = []

#     # calculate distances on each path and weights
#     for path in paths:
#         e, w, t = process_path(path=path, subtree=subtree_obj)
#         total_path_values.append(t)

#         print(f"The path {path} edge lengths are:")
#         for edge, weight in w.items():
#             if edge not in upd_edge_lengths:
#                 upd_edge_lengths[edge] = weight * t
#                 print(f"Edge: {edge}, Length: {upd_edge_lengths[edge]}")

#     # take the average branch length
#     avg_branch_l = sum(total_path_values) / len(total_path_values)

#     # normalize the updated edge lengths by the average branch length
#     for edge in upd_edge_lengths:
#         upd_edge_lengths[edge] = avg_branch_l

#     # now update the edge lengths on the tree
#     for node in obs_tree_obj.traverse_postorder():
#         if node.get_label() in upd_edge_lengths:
#             val = upd_edge_lengths[node.get_label()]
#             node.set_edge_length(length=val)

#     print("Final averaged tree")


def extract_subtree_paths(subtree):
    paths_to_node_list = []

    child_nodes = list(subtree.root.child_nodes())  

    if len(child_nodes) < 2:
        print("Warning: Less than two child nodes!")
    else:
        left_child = subtree.root.child_nodes()[0]
        right_child = subtree.root.child_nodes()[1]

        right_paths = get_path_to_leaf(right_child)
        paths_to_node_list.append(right_paths)

        left_paths = get_path_to_leaf(left_child)
        paths_to_node_list.append(left_paths)

        return paths_to_node_list, right_paths, left_paths

# def average_subtree(tree_obj, subtree_side):
#     for n in subtree_side.traverse_postorder(leaves=False):
#         print("Subtree Calculations...")
#         print("\n")
#         print(f"Subtree Node Now: {n}")

#         subtree_small = subtree_side.extract_subtree(n)
#         child_nodes_v2 = list(subtree_small.root.child_nodes())
        
#         if len(child_nodes_v2) < 2:
#             print("Warning: Less than two child nodes!")
#             continue

#         # take paths_to_node_list, right paths and left paths of subtree
#         paths_to_node_list, _, _ = extract_subtree_paths(subtree_small)

#         # Flatten the list of paths
#         #all_paths = [path for sublist in paths_to_node_list for path in sublist]

#         # Update tree lengths using all paths
#         upd_tree_lengths(obs_tree_obj=tree_obj, paths=paths_to_node_list, subtree_obj=subtree_small)

# def average_subtree(tree_obj, subtree_side):
#     for n in subtree_side.traverse_postorder(leaves=False):
#         print("Subtree Calculations...")
#         print("\n")
#         print(f"Subtree Node Now: {n}")

#         subtree_small = subtree_side.extract_subtree(n)
#         child_nodes_v2 = list(subtree_small.root.child_nodes())
        
#         if len(child_nodes_v2) < 2:
#             print("Warning: Less than two child nodes!")
#             continue

#         # take paths_to_node_list, right paths and left paths of subtree
#         paths_to_node_list, _, _ = extract_subtree_paths(subtree_small)

#         # Iterate over all pairs of paths and update tree lengths
#         for i in range(len(paths_to_node_list)):
#             for j in range(i+1, len(paths_to_node_list)):
#                 for path_i in paths_to_node_list[i]:
#                     for path_j in paths_to_node_list[j]:
#                         upd_tree_lengths(obs_tree_obj=tree_obj, left_path=path_i, right_path=path_j, subtree_obj=subtree_small) 

def average_subtree(tree_obj, subtree_side):
    for n in subtree_side.traverse_postorder(leaves=False):
        print("Subtree Calculations...")
        print("\n")
        print(f"Subtree Node Now: {n}")

        subtree_small = subtree_side.extract_subtree(n)
        child_nodes_v2 = list(subtree_small.root.child_nodes())
        
        if len(child_nodes_v2) < 2:
            print("Warning: Less than two child nodes!")
            continue
        #print("now printing standard tree")
        #subtree.draw(show_labels=True)

        # take paths_to_node_list, right paths and left paths of subtree
        paths_to_node_list, right_paths, left_paths = extract_subtree_paths(subtree_small)
        right_avg_path_sub = paths_to_node_list[0][0] # right path
        left_avg_path_sub = paths_to_node_list[1][0] # left path   
        # take an averaged path from the paths_list from each side at random
        upd_tree_lengths(obs_tree_obj=tree_obj, left_path=left_avg_path_sub, right_path=right_avg_path_sub, subtree_obj=subtree_small) 

def traverse_and_run_average(tree_obj):
    for n in tree_obj.traverse_postorder(leaves=False):
        children_n = n.child_nodes()
        
        if n.is_root() and all(child.is_leaf() for child in children_n):
            print("All children are leaves, performing operation...")
            subtree = tree_obj.extract_subtree(n)
            subtree.draw()
            new_branch_length = subtree.avg_branch_length()
            print(new_branch_length)
            # now update the edge lengths on the tree
            for node in tree_obj.traverse_postorder(internal=False):
                node.set_edge_length(length=new_branch_length)
        elif n.is_root():
            continue

        parent = n.get_parent()
        if parent == tree_obj.root and n.is_leaf() and n.get_edge_length() == 0:
            print("Node with zero length right after root found")
            break
        print("\n")
        print(f"Node Now: {n}")

        subtree = tree_obj.extract_subtree(n)
        children = list(subtree.root.child_nodes())
        if len(children) < 2:
            print("Less than two child nodes!")
            continue  
        left_child = children[0]
        right_child = children[1]
        #print("now printing standard tree")
        #subtree.draw(show_labels=True)

        # take paths_to_node_list, right paths and left paths of subtree
        paths_to_node_list, right_paths, left_paths = extract_subtree_paths(subtree)

        right_totals = []
        
        for path in right_paths:
            _right_lengths, _right_weights, right_total = process_path(path, subtree)
            right_totals.append(right_total)
            print(f"Right Total Path Value for path {path}: {right_total}")

        # Check if all right path totals are equal
        if len(set(right_totals)) == 1:
            print("All right paths have equal total path values.")
            right_avg_path = paths_to_node_list[0][0] # right path
        else:
            print("WARNING:! Not all right paths have equal total path values.")
            print("Averaging subtree...")
            right_subtree = tree_obj.extract_subtree(right_child)
            #right_subtree.draw()
            average_subtree(tree_obj = tree_obj, subtree_side =right_subtree)
            print("\n")
            print("Right Subtree averaged...")
            subtree = tree_obj.extract_subtree(n)
            paths_to_node_list, right_paths, left_paths = extract_subtree_paths(subtree)
            right_avg_path = paths_to_node_list[0][0]

        # Process each left path and check the total path value
        left_totals = []
        for path in left_paths:
            _left_lengths, _left_weights, left_total = process_path(path, subtree)
            left_totals.append(left_total)
            print(f"Left Total Path Value for path {path}: {left_total}")

        # Check if all left path totals are equal
        if len(set(left_totals)) == 1:
            print("All left paths have equal total path values.")
            left_avg_path = paths_to_node_list[1][0] # left path
        else:
            print("WARNING:! Not all left paths have equal total path values.")
            print("Averaging subtree...")
            left_subtree = tree_obj.extract_subtree(left_child)
            average_subtree(tree_obj, left_subtree)
            print("\n")
            print("Left Subtree averaged...")
            subtree = tree_obj.extract_subtree(n)
            paths_to_node_list, right_paths, left_paths = extract_subtree_paths(subtree)
            left_avg_path = paths_to_node_list[1][0]
            
        # take an averaged path from the paths_list from each side at random
        for rpath in right_paths:
            upd_tree_lengths_two(obs_tree_obj=tree_obj, left_path=left_avg_path, right_path=rpath, subtree_obj=subtree) 
        for lpath in left_paths:
            upd_tree_lengths_two(obs_tree_obj=tree_obj, left_path=lpath, right_path=rpath, subtree_obj=subtree) 
            

def transform_data(data):
    data_transformed = []
    
    for i in range(len(data) - 1):
        current_x, current_y = data[i]
        next_x, _ = data[i+1]
        
        for x in range(int(current_x), int(next_x)):
            data_transformed.append((x, current_y))

    last_x, last_y = data[-1]
    # Assuming that we need to add one more point after the last x based on the example provided
    data_transformed.append((int(last_x), last_y))
    data_transformed.append((int(last_x) + 1, last_y))
    
    return data_transformed

def normalise_data(data):
    data_norm = []
    
    for i in range(len(data)):
        current_x, current_y = data[i]
        
        # Normalise the 0th element of the tuple with the last 0th element of the data
        current_x = current_x / data[-2][0]
        
        data_norm.append((current_x, current_y))
    
    return data_norm
