import ast

def load_weights(filepath):
    with open(filepath, 'r') as f:
        bin_edges = ast.literal_eval(f.readline().strip())
        ratios = ast.literal_eval(f.readline().strip())

    if len(bin_edges) != len(ratios) + 1:
        raise ValueError("Mismatch between bin edges and weight count.")
    
    return bin_edges, ratios

def get_weight_from_file(energy, bin_edges, ratios):
    for i in range(len(ratios)):
        if bin_edges[i] <= energy < bin_edges[i + 1]:
            #weighted_energy = energy*ratios[i]
            #print(energy, ratios[i], weighted_energy)
            return ratios[i]
    return 1.0  # default weight if energy is out of bounds
