# Ligand-Based 3D screening with VSFlow
Ligand-based virtual screening based on the VSFlow pipeline ([Jung et al, 2023](https://doi.org/10.1186/s13321-023-00703-1))


## Installation

Create a conda environment and install the following dependencies

```bash
conda create -n vsflow python=3.11
cd vsflow-ligand-based-3d-screening
pip install -r requirements.txt
conda install -c conda-forge pymol-open-source
```

## Usage
We provide an example usage. We recommend using the [smiles-to-3d](https://github.com/ersilia-os/smiles-to-3d) repository to generate the sdf files. They must have a Name property.

```bash
python src/vsflow3d.py -q example/one_molecule.sdf -d molecules.sdf -o results
```

The results will be stored under the specified folder, providing the following metrics in a csv file:
* ComboScore
* Shape Similarity
* FP3D Similarity 

## License
Ersilia code is licensed under a GPLv3 License. The packages used are licensed according to their original licenses.
