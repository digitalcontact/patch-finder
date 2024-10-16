# Surface Residue Patch Finder

This tool identifies surface amino acids from a protein structure (PDB file) and finds combinatory patches of amino acids using graph theory. The patches are defined based on a distance threshold between residues and a group size. It uses the BioPython library to parse protein structures, compute solvent-accessible surface areas (SASA), and NetworkX to identify residue contact networks.

## Features

- Parses a PDB file to extract amino acids.
- Calculates solvent-accessible surface areas (SASA) to determine surface residues.
- Builds a graph of surface residues where edges are added based on a distance threshold.
- Identifies cliques of a specified size, which represent groups of residues that are all within the distance threshold from each other.
- Outputs valid combinations of surface residues as residue IDs.
- Provides visualization of the output data using histograms, box plots, and bar plots.

## Requirements

- Python 3.11
- `tqdm`
- `networkx`
- `numpy`
- `biopython`
- `matplotlib`

Install the dependencies using:

```bash
pip install tqdm networkx numpy biopython matplotlib
```

## Usage

### Surface Patch Finder

Run the tool using the following command:

```bash
python graph.py <pdb_file> <min_X> <max_X> <min_Y> <max_Y> [-o output_file]
```

### Arguments:

- `pdb_file`: Path to the input PDB file.
- `min_X`: Minimum number of amino acids in a patch (group size).
- `max_X`: Maximum number of amino acids in a patch (group size).
- `min_Y`: Minimum distance threshold in angstroms for residues to be considered part of the same patch.
- `max_Y`: Maximum distance threshold in angstroms for residues to be considered part of the same patch.
- `-o, --output`: Optional. Output file to store valid combinations. Defaults to `output.csv`.

### Example:

```bash
python graph.py example/6os0_D.pdb 3 5 4.0 8.0 -o patches.csv
```

This example identifies patches of 3 to 5 surface residues, where each residue is within 4.0 Å to 8.0 Å of another in the same patch. The valid patches will be written to `patches.csv`.

### Output:

The tool will output the valid combinations of surface residues to the specified file. Each line in the output will list the number of amino acids, average distance, and residue IDs of a patch, separated by commas. The residue ID format is:

```
<residue_name><residue_number><insertion_code>_<chain_id>
```

For example:

```
3, 5.67, ARG12A_B, GLY15B_C, LEU20D_A
```

### Plotting the Output

To visualize the output data, use the `plot_output.py` script:

```bash
python plot_output.py <csv_file> <output_image_histogram> <output_image_boxplot> <output_image_bar>
```

### Arguments:

- `csv_file`: Path to the CSV file containing the output data.
- `output_image_histogram`: Output image file to save the histogram plot.
- `output_image_boxplot`: Output image file to save the box plot.
- `output_image_bar`: Output image file to save the bar plot.

### Example:

```bash
python plot_output.py patches.csv histogram.png boxplot.png barplot.png
```

This example generates a histogram, box plot, and bar plot from the data in `patches.csv` and saves them as `histogram.png`, `boxplot.png`, and `barplot.png` respectively.

## How It Works

1. **Parsing the PDB**: The PDB file is parsed using BioPython’s `PDBParser`.
2. **Surface Residue Detection**: Solvent-accessible surface area (SASA) is calculated using the Shrake-Rupley algorithm. Only residues with non-zero SASA are considered as surface residues.
3. **Graph Construction**: A graph is constructed where each node represents a surface residue, and edges are added between nodes whose residues are within the specified distance threshold `min_Y` to `max_Y`.
4. **Clique Finding**: Using NetworkX’s `find_cliques` function, all cliques of size between `min_X` and `max_X` are identified, where each residue in the clique is connected to every other residue within the distance threshold.
5. **Output**: Valid combinations of residues (cliques) are written to the output file in a readable format.
6. **Visualization**: The `plot_output.py` script reads the output CSV file and generates visualizations (histogram, box plot, and bar plot) of the average distances and number of amino acids in the groups.
