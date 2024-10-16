# Surface Residue Patch Finder

This tool identifies surface amino acids from a protein structure (PDB file) and finds combinatory patches of amino acids using graph theory. The patches are defined based on a distance threshold between residues and a group size. It uses the BioPython library to parse protein structures, compute solvent-accessible surface areas (SASA), and NetworkX to identify residue contact networks.

## Features

- Parses a PDB file to extract amino acids.
- Calculates solvent-accessible surface areas (SASA) to determine surface residues.
- Builds a graph of surface residues where edges are added based on a distance threshold.
- Identifies cliques of a specified size, which represent groups of residues that are all within the distance threshold from each other.
- Outputs valid combinations of surface residues as residue IDs.

## Requirements

- Python 3.x
- `tqdm`
- `networkx`
- `numpy`
- `biopython`

Install the dependencies using:

```bash
pip install tqdm networkx numpy biopython
```

## Usage

Run the tool using the following command:

```bash
python surface_patch_finder.py <pdb_file> <X> <Y> [-o output_file]
```

### Arguments:

- `pdb_file`: Path to the input PDB file.
- `X`: Number of amino acids in a patch (group size).
- `Y`: Distance threshold in angstroms for residues to be considered part of the same patch.
- `-o, --output`: Optional. Output file to store valid combinations. Defaults to `output.txt`.

### Example:

```bash
python surface_patch_finder.py 1abc.pdb 3 5.0 -o patches.txt
```

This example identifies patches of 3 surface residues, where each residue is within 5.0 Å of another in the same patch. The valid patches will be written to `patches.txt`.

### Output:

The tool will output the valid combinations of surface residues to the specified file. Each line in the output will list the residue IDs of a patch, separated by commas. The residue ID format is:

```
<residue_name><residue_number><insertion_code>_<chain_id>
```

For example:

```
ARG12A_B, GLY15B_C, LEU20D_A
```

## How It Works

1. **Parsing the PDB**: The PDB file is parsed using BioPython’s `PDBParser`.
2. **Surface Residue Detection**: Solvent-accessible surface area (SASA) is calculated using the Shrake-Rupley algorithm. Only residues with non-zero SASA are considered as surface residues.
3. **Graph Construction**: A graph is constructed where each node represents a surface residue, and edges are added between nodes whose residues are within the specified distance threshold `Y`.
4. **Clique Finding**: Using NetworkX’s `find_cliques` function, all cliques of size `X` are identified, where each residue in the clique is connected to every other residue within the distance threshold.
5. **Output**: Valid combinations of residues (cliques) are written to the output file in a readable format.
