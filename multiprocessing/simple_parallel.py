
import sys
import os
import mdtraj as md 
import time
from multiprocessing import Pool, Process, Value, Array




USAGE_STR = """

TOP="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/DOR-inactive-naltrindole-unpublished/condition-naltrindole-bound/step5_assembly.pdb"
TRAJ="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/DOR-inactive-naltrindole-unpublished/condition-naltrindole-bound/rep_1/Prod_0/Prod_0_reimaged.nc"
cd /scratch/PI/rondror/akma327/github/MDContactNetworks/multiprocessing
python simple_parallel.py $TOP $TRAJ

TOP="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/B2AR-inactive-carazalol-cell2013b/condition-12/step5_assembly_strip_waters.pdb"
TRAJ="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutput/B2AR-inactive-carazalol-cell2013b/condition-12/rep_1/Prod_0/Prod_0_rewrapped_strip_waters.dcd"
cd /scratch/PI/rondror/akma327/github/MDContactNetworks/multiprocessing
python simple_parallel.py $TOP $TRAJ


"""

def count_frames(input_args):
	i, tf, distance, water_index, residue_indices = input_args
	print(i, tf, distance, water_index, residue_indices)
	return (i, len(tf))

def parallel(input_args):
	print("Begin parallelizing...")
	tic = time.clock()
	
	pool = Pool(processes=4)
	result = pool.map(count_frames, input_args)

	toc = time.clock()
	print("Parallel Time: " + str(toc-tic))
	print(result)


def driver(TOP, TRAJ):
	
	t = md.load(TRAJ, top=TOP)
	input_args = []

	for i, x in enumerate(range(0, len(t), 5)):
		input_args.append((i, t[x:x+5], 1.0, '25', '272 277 303'))

	for i in input_args: print i 
	
	parallel(input_args)
	

if __name__ == "__main__":
	(TOP, TRAJ) = (sys.argv[1], sys.argv[2])
	driver(TOP, TRAJ)