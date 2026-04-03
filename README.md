# GetContacts

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Application for efficiently computing non-covalent contact networks in molecular structures and MD simulations. Following example computes all salt bridges, pi-cation, aromatic, and hydrogen bond interactions in a trajectory:
```bash
get-dynamic-contacts --topology my_top.psf \
                     --trajectory my_traj.dcd \
                     --itypes sb pc ps ts hb \
                     --output my_contacts.tsv
```
The output, `my_contacts.tsv`, is a tab-separated file where each line (except the first two) records frame, type, and atoms involved in an interaction:
```
# total_frames:20000 interaction_types:sb,pc,ps,ts,hb
# Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]]
0   sb     C:GLU:21:OE2    C:ARG:86:NH2
0   ps     C:TYR:36:CG     C:TRP:108:CG
0   ts     A:TYR:36:CG     A:TRP:108:CG
0   hbss   A:GLN:53:NE2    A:GLN:69:OE1
1   hbss   A:GLU:60:OE2    A:SER:56:OG
1   hbbb   C:LEU:87:N      C:LEU:83:O
1   hbbb   C:ARG:88:N      C:LEU:87:N
2   hbsb   A:LYS:28:N      A:HIS:27:ND1
2   hbsb   A:ASP:52:OD2    A:PHE:48:O
2   wb2    C:ASN:110:O     B:ARG:73:NH1    W:TIP3:1524:OH2    W:TIP3:506:OH2
2   wb2    C:ASN:110:O     C:SER:111:OG    W:TIP3:1524:OH2    W:TIP3:2626:OH2
3   wb     A:ASP:100:OD1   A:ILE:67:O      W:TIP3:6762:OH2
3   wb     A:ASP:100:OD1   B:ASN:105:ND2   W:TIP3:9239:OH2
4   sb     A:GLU:47:OE2    A:LYS:33:NZ
4   pc     A:LYS:9:NZ      A:TYR:21:CG
4   hbbb   A:ILE:12:N      A:THR:20:O
5   vdw    A:LYS:164:CG    A:THR:161:OG1
5   vdw    A:LYS:164:CD    A:THR:161:OG1
6   vdw    A:ALA:156:O     A:PHE:159:C
6   hbls   A:EJ4:1219:OAD  A:TYR:129:OH
6   hbls   A:ASP:128:OD2   A:EJ4:1219:OAQ
7   hblb   C:LEU:87:N      A:EJ4:1219:OAQ
7   lwb    A:EJ4:1219:NAN  A:LYS:214:NZ    A:HOH:1316:O
7   lwb2   A:EJ4:1219:NAN  A:LYS:214:NZ    A:HOH:1316:O       A:HOH:1315:O
7   lwb2   A:EJ4:1219:NAN  A:TYR:129:OH    A:HOH:1316:O       A:HOH:1311:O
...
```
Interactions that involve more than two atoms (i.e. water bridges and extended water bridges) have extra columns to denote the identities of the water molecules. For simplicity, all stacking and pi-cation interactions involving an aromatic ring will be denoted by the CG atom. 

Interaction types are denoted by the following abbreviations:
* **sb** - salt bridges 
* **pc** - pi-cation 
* **ps** - pi-stacking
* **ts** - T-stacking
* **vdw** - van der Waals
* **hb** - Hydrogen bond subtypes
  * **hbbb** - Backbone-backbone hydrogen bonds
  * **hbsb** - Backbone-sidechain hydrogen bonds
  * **hbss** - Sidechain-sidechain hydrogen bonds
  * **wb** - Water-mediated hydrogen bond
  * **wb2** - Extended water-mediated hydrogen bond
  * **hblb** - Ligand-backbone hydrogen bonds
  * **hbls** - Ligand-sidechain hydrogen bonds
  * **lwb** - Ligand water-mediated hydrogen bond
  * **lwb2** - Ligand extended water-mediated hydrogen bond

Generated contact-list files are useful as inputs to visualization and analysis tools that operate on interaction-networks:
 * [Flareplot](https://gpcrviz.github.io/flareplot) - Framework for analyzing interaction networks based on circular diagrams
 * [MDCompare](MDCompare) - Heatmap fingerprints revealing groups of similar interactions in multiple MD trajectories
 * [TICC](https://github.com/davidhallac/TICC) - Changepoint detection algorithm to identifying significant events in the dynamic contact network
 * [NetworkAnalysis](Applications) - Compute residue contact frequencies in a simulation, analyze contact network graphs, visualize in PyMol

## Dependencies

GetContacts requires:
* [vmd-python](https://github.com/Eigenstate/vmd-python) — **only available via conda-forge**, not PyPI
* numpy, scipy, matplotlib, scikit-learn, pandas, seaborn, cython

## Installation

### Recommended: conda + pip

`vmd-python` must be installed via conda-forge before the package itself. Create a dedicated environment:

```bash
conda create -n getcontacts python=3.12
conda activate getcontacts
conda install -c conda-forge vmd-python
```

Then install getcontacts with pip (editable for development, or from git):

```bash
# From a local clone (editable — source changes take effect immediately)
git clone https://github.com/getcontacts/getcontacts
pip install -e ./getcontacts

# Or directly from GitHub (no clone needed)
pip install git+https://github.com/getcontacts/getcontacts
```

### As a dependency of another package

Add to your `pyproject.toml` or `setup.cfg`:

```toml
dependencies = [
    "getcontacts @ git+https://github.com/getcontacts/getcontacts",
]
```

> **Note:** The consuming environment must still have `vmd-python` installed via
> `conda install -c conda-forge vmd-python` before running `pip install`.

### Conda build (local)

```bash
conda install conda-build
conda build conda-recipe/
conda install --use-local getcontacts
```

### Verifying the installation

After installation, the following commands should all be available on your `$PATH`:

```bash
get-dynamic-contacts --help
get-static-contacts --help
get-contact-frequencies --help
```

To run the full test suite:

```bash
cd getcontacts
pytest
```

## Simulation and structure file format

GetContacts is compatible with all topology and reimaged trajectory file formats readable by [VMD](https://www-s.ks.uiuc.edu/Research/vmd/).

## Running the Code

### Command line arguments

Required Arguments:

    --topology STRING 
        Path to topology
    --trajectory STRING
        Path to simulation trajectory fragment
    --output STRING
        Path to output file
    --itypes STRING [STRING ...]
        Specifies interaction types (see above)
          all, sb, pc, ps, ts, vdw, hb

Optional Arguments:

    --cores INT 
        Number of CPU cores for parallelization [default = 6]
    --ligand STRING [STRING ...]
        Resname of ligand molecule(s) [default = ""]
    --sele STRING 
        VMD selection query to compute contacts in specified region of protein 
        [default = "protein"]
    --sele2 STRING
        Second VMD selection query to compute contacts between two regions of the protein
        [default = None]
    --solv STRING 
        Solvent identifier in simulation [default = "TIP3"]

Arguments for adjusting geometric criteria:

    --sb_cutoff_dist FLOAT
        Cutoff for distance between anion and cation atoms [default = 4.0 Angstroms]
    --pc_cutoff_dist FLOAT
        Cutoff for distance between cation and centroid of aromatic ring [default = 6.0 Angstroms]
    --pc_cutoff_ang FLOAT
        Cutoff for angle between normal vector projecting from aromatic plane and vector from 
        aromatic center to cation atom [default = 60 degrees]
    --ps_cutoff_dist FLOAT
        Cutoff for distance between centroids of two aromatic rings [default = 7.0 Angstroms]
    --ps_cutoff_ang FLOAT
        Cutoff for angle between the normal vectors projecting from each aromatic plane 
        [default = 30 degrees]
    --ps_psi_ang FLOAT
        Cutoff for angle between normal vector projecting from aromatic plane 1 and vector between 
        the two aromatic centroids [default = 45 degrees]
    --ts_cutoff_dist FLOAT
        Cutoff for distance between centroids of two aromatic rings [default = 5.0 Angstroms]
    --ts_cutoff_ang FLOAT
        Cutoff for angle between the normal vectors projecting from each aromatic plane minus 90 
        degrees [default = 30 degrees]
    --ts_psi_ang FLOAT
        Cutoff for angle between normal vector projecting from aromatic plane 1 and vector between 
        the two aromatic centroids [default = 45 degrees]
    --hbond_cutoff_dist FLOAT
        Cutoff for distance between donor and acceptor atoms [default = 3.5 Angstroms]
    --hbond_cutoff_ang FLOAT
        Cutoff for angle between donor hydrogen acceptor [default = 70 degrees in simulation, 180 in structures]
    --hbond_res_diff INT
        Minimum residue distance for which to consider computing hbond interactions 
        [default = 0 in simulations 1 in structures]
    --vdw_epsilon FLOAT
        Amount of padding for calculating vanderwaals contacts [default = 0.5 Angstroms]
    --vdw_res_diff INT
        Minimum residue distance for which to consider computing vdw interactions [default = 2]

   
### Examples

Compute all pi-cation, pi-stacking, and van der Waals contacts throughout an entire simulation:

    get-dynamic-contacts --topology TOP.psf \
                         --trajectory TRAJ.dcd \
                         --output output_pc_ps_vdw.tsv \
                         --itypes pc ps vdw

Compute salt bridges and hydrogen bonds for residues 100 to 160, including those that involves a ligand:

    get-dynamic-contacts --topology TOP.pdb \
                         --trajectory TRAJ.nc \
                         --output loop-ligand_sb_hb.tsv \
                         --cores 12 \
                         --solv IP3 \
                         --sele "chain A and resid 100 to 160" \
                         --ligand EJ4 \
                         --itypes sb hb 

Compute all interactions formed at a protein-protein interface throughout simulation:

    get-dynamic-contacts --topology TOP.pdb \
                         --trajectory TRAJ.nc \
                         --output protein_interface_all.tsv \
                         --cores 12 \
                         --sele "chain A" \
                         --sele2 "chain B" \
                         --itypes all

Locate salt bridges and hydrogen bonds using modified distance cutoffs:

    get-dynamic-contacts --topology TOP.mae \
                         --trajectory TRAJ.dcd \
                         --output output_sb_hb.tsv \
                         --sb_cutoff_dist 5.0 \
                         --hbond_cutoff_dist 4.5 \
                         --itypes sb hb

For further examples, see the [GetContacts main page](https://getcontacts.github.io).

### Programmatic usage

After `pip install`, the tools are importable as a Python package:

```python
from getcontacts.get_dynamic_contacts import main as get_dynamic_contacts
from getcontacts.get_contact_frequencies import main as get_contact_frequencies

get_dynamic_contacts([
    "--topology",  "my_top.pdb",
    "--trajectory", "my_traj.dcd",
    "--output",    "contacts.tsv",
    "--itypes",    "hb", "sb",
    "--cores",     "8",
])

get_contact_frequencies([
    "--input_files",  "contacts.tsv",
    "--output_file",  "frequencies.tsv",
])
```

Or via subprocess (recommended when calling from within another package's
multiprocessing workers, to keep VMD state fully isolated):

```python
import subprocess

subprocess.run([
    "get-dynamic-contacts",
    "--topology",  "my_top.pdb",
    "--trajectory", "my_traj.dcd",
    "--output",    "contacts.tsv",
    "--itypes",    "all",
    "--cores",     "8",
], check=True)
```

## Python version compatibility

GetContacts requires **Python ≥ 3.7**. Internal parallelization uses the
`forkserver` multiprocessing context (falling back to `spawn`) so that VMD's
C library is never duplicated into worker processes via `fork`. This resolves
memory allocation errors that occurred on Python 3.12+ with the previous
`fork`-based approach.

## For developers

Unit tests are located in `unittests` subfolders and use the python `unittest` module. To run just the unit-tests, type

    python -m unittest

Integration tests verifying the functionality of executables are in the `tests`-folder. To run all tests make sure pytest is installed (e.g. `pip install pytest`) and run the following from the root:

    pytest
    
Pull-requests won't be accepted before all tests pass, but if there are any problems we are happy to help work them out.

The code aims to be [PEP 8](http://pep8.org/) compliant, but pull-requests wont be rejected if they're not. 

## Citation

When using GetContacts for publication, please cite: 

``https://getcontacts.github.io/``
