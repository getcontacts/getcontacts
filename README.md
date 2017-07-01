# MDContactNetworks
Library for computing dynamic non-covalent contact networks in proteins throughout MD Simulation


## Running the Code

The code is organized into two steps:

1. Computing non-covalent contacts throughout a single fragment of an MD Simulation
   
   __Input:__ Simulation fragment trajectory, Topology, Output, Contact Type. Optional Flags: Stride, Solvent, Chain, Ligand

   
   __Output:__ List of non-covalent contacts of specified type for each frame of simulation fragment 


2. Visualizing dynamic contact networks

   __Input:__ Directory containing non-covalent contacts computed for each fragment of simulation

   
   __Output:__ JSON file that has keys representing interactions between a pair of atoms. Values being the time points that the interaction was formed in simulation. Formatted to be visualized with flareplots (https://gpcrviz.github.io/flareplot/)


