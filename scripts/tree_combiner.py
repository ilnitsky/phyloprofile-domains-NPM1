#!/usr/bin/env python3
"""
Fresh script to combine two Newick trees
"""

import sys
import re

def read_newick_file(filename):
    """Read a Newick tree from file"""
    with open(filename, 'r') as f:
        return f.read().strip()

def parse_simple_newick(newick_str):
    """Simple parser to extract leaf names from Newick format"""
    # Remove the trailing semicolon
    newick_str = newick_str.rstrip(';')
    
    # Find all leaf names (text before commas or closing parentheses)
    # This regex finds sequences that look like leaf names with IDs
    leaf_pattern = r'([A-Za-z_][A-Za-z0-9_.-]*-[A-Za-z0-9_.-]+)'
    leaves = re.findall(leaf_pattern, newick_str)
    
    return leaves

def combine_newick_trees(tree1_str, tree2_str):
    """Combine two Newick trees by creating a simple polytomy"""
    # Remove trailing semicolons
    tree1_clean = tree1_str.rstrip(';')
    tree2_clean = tree2_str.rstrip(';')
    
    # Create a combined tree with both as sister groups
    combined = f"({tree1_clean},{tree2_clean});"
    
    return combined

def main():
    if len(sys.argv) != 4:
        print("Usage: python tree_combiner.py <tree1.nwk> <tree2.nwk> <output.nwk>")
        sys.exit(1)
    
    tree1_file = sys.argv[1]
    tree2_file = sys.argv[2]
    output_file = sys.argv[3]
    
    print(f"Reading tree 1: {tree1_file}")
    tree1 = read_newick_file(tree1_file)
    
    print(f"Reading tree 2: {tree2_file}")
    tree2 = read_newick_file(tree2_file)
    
    print("Combining trees...")
    combined_tree = combine_newick_trees(tree1, tree2)
    
    print(f"Writing combined tree to: {output_file}")
    with open(output_file, 'w') as f:
        f.write(combined_tree)
    
    # Show some stats
    leaves1 = parse_simple_newick(tree1)
    leaves2 = parse_simple_newick(tree2)
    leaves_combined = parse_simple_newick(combined_tree)
    
    print(f"Tree 1 has {len(leaves1)} leaves")
    print(f"Tree 2 has {len(leaves2)} leaves") 
    print(f"Combined tree has {len(leaves_combined)} leaves")
    print("Done!")

if __name__ == "__main__":
    main() 