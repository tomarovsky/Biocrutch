from ete3 import Tree
import pandas as pd
import re
import argparse

def extract_newick_from_nexus(nexus_file):
    """Extracts tree string from NEXUS file"""
    with open(nexus_file) as f:
        for line in f:
            if line.strip().startswith("UTREE"):
                return line.split("=", 1)[1].strip().rstrip(";")
    raise ValueError("Failed to find tree string (UTREE) in file")

def process_tree(input_file, output_file):
    """Main tree processing function"""
    # Extract tree string
    newick_str = extract_newick_from_nexus(input_file)
    
    # Extract HPD data
    hpd_data = []
    matches = re.finditer(r'\)\s*\[\&95%HPD=\{([^}]+)\}\]', newick_str)
    for match in matches:
        hpd_values = list(map(float, match.group(1).split(',')))
        hpd_data.append((hpd_values[0], hpd_values[1]))
    
    # Create clean tree without annotations
    clean_str = re.sub(r'\[\&95%HPD=\{([^}]+)\}\]', '', newick_str)
    tree = Tree(clean_str + ";", format=1)
    
    # Collect node data
    rows = []
    for i, node in enumerate([n for n in tree.traverse("postorder") if not n.is_leaf()]):
        dist_to_root = node.get_farthest_leaf()[1]
        hpd_lower, hpd_upper = hpd_data[i] if i < len(hpd_data) else (dist_to_root, dist_to_root)
        
        # Create node description
        leaves = node.get_leaf_names()
        desc = f"MRCA of {' and '.join(sorted(leaves))}" if len(leaves) == 2 else \
               f"MRCA of {', '.join(sorted(leaves[:2]))} and others"
        
        rows.append({
            'Node_ID': i+1,
            'Description': desc,
            'Date': round(dist_to_root, 6),
            'HPD_Lower': round(hpd_lower, 6),
            'HPD_Upper': round(hpd_upper, 6)
        })

    # Save results
    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Process phylogenetic tree with HPD intervals')
    parser.add_argument('-i', '--input', required=True, help='Input NEXUS file (e.g., FigTree.tre)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    args = parser.parse_args()
    
    # Run processing
    process_tree(args.input, args.output)
