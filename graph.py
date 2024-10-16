import argparse
from tqdm import tqdm
import networkx as nx
import numpy as np
from Bio.PDB import PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa

def get_CA_coords(residue):
    if residue.has_id('CA'):
        return residue['CA'].get_coord()
    else:
        return None

def get_residue_id(residue):
    chain_id = residue.get_parent().id
    resseq = residue.get_id()[1]
    icode = residue.get_id()[2].strip()
    resname = residue.get_resname()
    return f"{resname}{resseq}{icode}_{chain_id}"

def main():
    parser = argparse.ArgumentParser(description='Identify surface amino acids and combinatory patches using graphs')
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('X', type=int, help='Number of amino acids in a patch (group size)')
    parser.add_argument('Y', type=float, help='Distance threshold (A)')
    parser.add_argument('-o', '--output', help='Output file to store valid combinations', default='output.txt')
    args = parser.parse_args()

    pdb_file = args.pdb_file
    X = args.X
    Y = args.Y
    output_file = args.output

    # Parse the PDB file
    parser_pdb = PDBParser(QUIET=True)
    structure = parser_pdb.get_structure('protein', pdb_file)

    # Compute solvent accessible surface area (SASA)
    sr = ShrakeRupley()
    sr.compute(structure, level='R')  # Compute at residue level

    # Collect surface residues and their coordinates
    surface_residues = []
    residue_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if not is_aa(residue):
                    continue  # Skip non-amino acid residues
                if residue.has_id('CA'):
                    sasa = residue.sasa
                    if sasa > 0.0:
                        surface_residues.append(residue)
                        coord = get_CA_coords(residue)
                        residue_coords.append(coord)

    # Check if there are enough residues to form combinations
    if len(surface_residues) < X:
        print(f"Not enough surface residues to form combinations of {X}.")
        return

    # Build graph
    print("Building residue contact graph...")
    G = nx.Graph()
    for idx, residue in enumerate(surface_residues):
        G.add_node(idx, residue=residue)

    # Add edges between residues within distance Y
    for i in tqdm(range(len(surface_residues))):
        for j in range(i+1, len(surface_residues)):
            dist = np.linalg.norm(residue_coords[i] - residue_coords[j])
            if dist <= Y:
                G.add_edge(i, j)

    print("Finding groups...")
    # Find cliques of size X
    cliques = [clique for clique in nx.find_cliques(G) if len(clique) == X]

    # Open the output file
    with open(output_file, 'w') as f_out:
        for clique in cliques:
            residues = [G.nodes[node]['residue'] for node in clique]
            res_ids = [get_residue_id(residue) for residue in residues]
            f_out.write(', '.join(res_ids) + '\n')

    print(f"Valid combinations have been written to {output_file}")

if __name__ == '__main__':
    main()
