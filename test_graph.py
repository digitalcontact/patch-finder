import unittest
import numpy as np
import networkx as nx
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from unittest.mock import MagicMock, patch
from graph import get_CA_coords, get_residue_id, main

class TestSurfaceResidueIdentification(unittest.TestCase):

    def setUp(self):
        # Create a mock residue with a CA atom
        self.residue_with_ca = Residue((' ', 1, ' '), 'GLY', ' ')
        ca_atom = Atom('CA', np.array([1.0, 2.0, 3.0]), 1.0, 1.0, ' ', 'CA', 1, element='C')
        self.residue_with_ca.add(ca_atom)

        # Create a mock residue without a CA atom
        self.residue_without_ca = Residue((' ', 2, ' '), 'ALA', ' ')

    def test_get_CA_coords_with_CA(self):
        # Test that get_CA_coords returns coordinates when CA atom is present
        ca_coords = get_CA_coords(self.residue_with_ca)
        self.assertTrue(np.array_equal(ca_coords, np.array([1.0, 2.0, 3.0])))

    def test_get_CA_coords_without_CA(self):
        # Test that get_CA_coords returns None when CA atom is not present
        ca_coords = get_CA_coords(self.residue_without_ca)
        self.assertIsNone(ca_coords)

    def test_get_residue_id(self):
        # Add the residue to a chain, and chain to a model/structure for parent association
        chain = Chain('A')
        chain.add(self.residue_with_ca)

        # Test that the residue ID is generated correctly
        residue_id = get_residue_id(self.residue_with_ca)
        self.assertEqual(residue_id, 'GLY1_A')

    def test_graph_construction(self):
        # Mock surface residues and coordinates
        surface_residues = [self.residue_with_ca, self.residue_with_ca]  # 2 residues
        coords = [np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 4.0])]  # within distance Y = 5

        # Build the graph
        G = nx.Graph()
        for idx, residue in enumerate(surface_residues):
            G.add_node(idx, residue=residue)

        # Add edges between residues within distance 5
        for i in range(len(surface_residues)):
            for j in range(i + 1, len(surface_residues)):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist <= 5:
                    G.add_edge(i, j, distance=dist)

        # Test that graph contains correct number of nodes and edges
        self.assertEqual(len(G.nodes), 2)
        self.assertEqual(len(G.edges), 1)

    def test_finding_cliques(self):
        # Create a graph with 3 nodes forming a triangle
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 0)])  # Fully connected 3-node graph

        # Find cliques of size 3
        cliques = [clique for clique in nx.find_cliques(G) if len(clique) == 3]

        # Test that we found exactly one clique of size 3
        self.assertEqual(len(cliques), 1)
        self.assertEqual(len(cliques[0]), 3)

    @patch('graph.PDBParser')
    @patch('graph.ShrakeRupley')
    @patch('graph.open', new_callable=unittest.mock.mock_open)
    @patch('graph.csv.writer')
    def test_full_pipeline(self, mock_csv_writer, mock_open, MockShrakeRupley, MockPDBParser):
        # Mock structure creation
        structure = Structure("test_structure")
        model = Model(0)
        chain = Chain('A')

        # Adding residues
        chain.add(self.residue_with_ca)
        model.add(chain)
        structure.add(model)

        # Mock PDBParser and ShrakeRupley
        mock_parser = MockPDBParser.return_value
        mock_parser.get_structure.return_value = structure
        mock_sr = MockShrakeRupley.return_value

        # Mock SASA computation
        self.residue_with_ca.sasa = 1.0  # Set a non-zero SASA value

        # Run the main function with mock arguments
        test_args = ['graph.py', 'dummy.pdb', '1', '3', '4.0', '8.0', '-o', 'output.csv']
        with patch('sys.argv', test_args):
            main()

        # Check if the structure was parsed
        mock_parser.get_structure.assert_called_once_with('protein', 'dummy.pdb')

        # Check if SASA was computed
        mock_sr.compute.assert_called_once_with(structure, level='R')

        # Check if the output file was written correctly
        mock_open.assert_called_once_with('output.csv', 'w', newline='')
        mock_csv_writer.return_value.writerow.assert_any_call(['Num_AA', 'Avg_Distance', 'Residues'])

if __name__ == '__main__':
    unittest.main()