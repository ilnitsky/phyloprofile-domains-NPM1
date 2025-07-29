#!/usr/bin/env python3
"""
Script to combine two phylogenetic trees by merging them at common taxonomic nodes.
"""

import re
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict
import argparse

class TreeNode:
    def __init__(self, name: str = ""):
        self.name = name
        self.children = []
        self.parent = None
        
    def add_child(self, child):
        child.parent = self
        self.children.append(child)
        
    def is_leaf(self):
        return len(self.children) == 0
    
    def get_leaves(self):
        """Get all leaf nodes under this node"""
        if self.is_leaf():
            return [self]
        leaves = []
        for child in self.children:
            leaves.extend(child.get_leaves())
        return leaves
    
    def to_newick(self):
        """Convert tree to Newick format"""
        if self.is_leaf():
            return self.name
        else:
            children_str = ",".join([child.to_newick() for child in self.children])
            if self.name:
                return f"({children_str}){self.name}"
            else:
                return f"({children_str})"

def parse_newick(newick_str: str) -> TreeNode:
    """Parse a Newick format string into a tree structure"""
    # Remove the final semicolon and whitespace
    newick_str = newick_str.strip().rstrip(';')
    
    # Stack to keep track of current nodes
    stack = []
    current_node = TreeNode()
    root = current_node
    i = 0
    
    while i < len(newick_str):
        char = newick_str[i]
        
        if char == '(':
            # Start a new internal node
            new_node = TreeNode()
            current_node.add_child(new_node)
            stack.append(current_node)
            current_node = new_node
            
        elif char == ')':
            # End current internal node, read its name
            current_node = stack.pop()
            i += 1
            # Read the node name after the closing parenthesis
            name_start = i
            while i < len(newick_str) and newick_str[i] not in '(),;':
                i += 1
            if i > name_start:
                current_node.children[-1].name = newick_str[name_start:i]
            i -= 1  # Adjust for the increment at the end of the loop
            
        elif char == ',':
            # Start a new sibling node
            new_node = TreeNode()
            current_node.parent.add_child(new_node)
            current_node = new_node
            
        else:
            # Read leaf node name
            name_start = i
            while i < len(newick_str) and newick_str[i] not in '(),;':
                i += 1
            if i > name_start:
                current_node.name = newick_str[name_start:i]
            i -= 1  # Adjust for the increment at the end of the loop
            
        i += 1
    
    return root

def extract_taxonomic_paths(node: TreeNode, current_path: List[str] = None) -> Dict[str, List[str]]:
    """Extract all taxonomic paths from leaf nodes to root"""
    if current_path is None:
        current_path = []
    
    paths = {}
    
    if node.name:
        current_path = current_path + [node.name]
    
    if node.is_leaf():
        # This is a leaf node (species/protein)
        if node.name:
            paths[node.name] = current_path
    else:
        # Recurse into children
        for child in node.children:
            child_paths = extract_taxonomic_paths(child, current_path)
            paths.update(child_paths)
    
    return paths

def find_common_ancestors(paths1: Dict[str, List[str]], paths2: Dict[str, List[str]]) -> Set[str]:
    """Find common taxonomic groups between two sets of paths"""
    all_taxa1 = set()
    all_taxa2 = set()
    
    for path in paths1.values():
        all_taxa1.update(path)
    
    for path in paths2.values():
        all_taxa2.update(path)
    
    return all_taxa1.intersection(all_taxa2)

def build_combined_tree(paths1: Dict[str, List[str]], paths2: Dict[str, List[str]]) -> TreeNode:
    """Build a combined tree from two sets of taxonomic paths"""
    # Combine all paths
    all_paths = {}
    all_paths.update(paths1)
    all_paths.update(paths2)
    
    # Build the tree bottom-up
    root = TreeNode()
    node_map = {"": root}  # Map from taxonomic name to node
    
    # Sort paths by length (deepest first) to build from leaves up
    sorted_items = sorted(all_paths.items(), key=lambda x: len(x[1]), reverse=True)
    
    for leaf_name, path in sorted_items:
        current_node = root
        
        # Traverse/create the path from root to leaf
        for i, taxon in enumerate(path):
            if taxon not in node_map:
                new_node = TreeNode(taxon)
                node_map[taxon] = new_node
                current_node.add_child(new_node)
                current_node = new_node
            else:
                current_node = node_map[taxon]
                # If this node already exists but isn't a child of current_node,
                # we need to handle the tree structure
                if current_node.parent is None and current_node != root:
                    # This node needs to be attached to the current path
                    if current_node not in [child for child in current_node.parent.children if child.parent]:
                        current_node.parent.add_child(current_node)
    
    return root

def simplify_tree_structure(root: TreeNode) -> TreeNode:
    """Simplify tree by removing redundant single-child internal nodes"""
    def simplify_node(node):
        # First, recursively simplify all children
        for child in node.children[:]:  # Use slice to avoid modification during iteration
            simplify_node(child)
        
        # If this node has exactly one child and both have names, 
        # we might want to keep the structure for taxonomic clarity
        # Only collapse if the parent node has no name
        if len(node.children) == 1 and not node.name:
            child = node.children[0]
            if node.parent:
                node.parent.children.remove(node)
                node.parent.add_child(child)
            else:
                # This is the root, replace it
                return child
        
        return node
    
    return simplify_node(root)

def combine_trees(tree1_str: str, tree2_str: str) -> str:
    """Main function to combine two Newick format trees"""
    print("Parsing tree 1...")
    tree1 = parse_newick(tree1_str)
    
    print("Parsing tree 2...")
    tree2 = parse_newick(tree2_str)
    
    print("Extracting taxonomic paths...")
    paths1 = extract_taxonomic_paths(tree1)
    paths2 = extract_taxonomic_paths(tree2)
    
    print(f"Tree 1 has {len(paths1)} leaf nodes")
    print(f"Tree 2 has {len(paths2)} leaf nodes")
    
    common_ancestors = find_common_ancestors(paths1, paths2)
    print(f"Found {len(common_ancestors)} common taxonomic groups")
    
    print("Building combined tree...")
    combined_tree = build_combined_tree(paths1, paths2)
    
    print("Simplifying tree structure...")
    simplified_tree = simplify_tree_structure(combined_tree)
    
    return simplified_tree.to_newick() + ";"

def main():
    parser = argparse.ArgumentParser(description="Combine two phylogenetic trees")
    parser.add_argument("-1", "--tree1", required=True, help="First tree in Newick format (file or string)")
    parser.add_argument("-2", "--tree2", required=True, help="Second tree in Newick format (file or string)")
    parser.add_argument("-o", "--output", help="Output file for combined tree")
    
    args = parser.parse_args()
    
    # Read tree 1
    try:
        with open(args.tree1, 'r') as f:
            tree1_str = f.read().strip()
    except FileNotFoundError:
        # Assume it's a tree string directly
        tree1_str = args.tree1
    
    # Read tree 2
    try:
        with open(args.tree2, 'r') as f:
            tree2_str = f.read().strip()
    except FileNotFoundError:
        # Assume it's a tree string directly
        tree2_str = args.tree2
    
    # Combine trees
    combined_tree = combine_trees(tree1_str, tree2_str)
    
    # Output result
    if args.output:
        with open(args.output, 'w') as f:
            f.write(combined_tree)
        print(f"Combined tree written to {args.output}")
    else:
        print("Combined tree:")
        print(combined_tree)

if __name__ == "__main__":
    main() 