#! /usr/bin/python3

import pandas as pd
pd.options.display.float_format = '{:.5f}'.format
import glob
import os
import shutil



def pull_motionpars(task, tsv, output_path):
	"""
	this function will pull the motionparameters the tsv passed to it and output a motion.1D 
	file to the corresponding task folder in /data/D2/1stlevelanalysis
	:param task: <string> name of task corresponding to tsv
	:param tsv: <string> confounds_timeseries.tsv containing motion parameters
	:param output_path: <string> output path for motion file
	"""
	motion1D_path = os.path.join(output_path, f'{task}_motion.tsv')
	confounds_tsv = pd.read_csv(tsv, sep='\t')
	motion_pars = confounds_tsv[['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']]
	motion_pars.to_csv(motion1D_path, header=None, index=None, sep='\t')
	print(f"\t created {motion1D_path} from {tsv}")


if __name__ == "__main__":
	import sys
	
	subject_directory = os.path.abspath(sys.argv[1])
	tsv_files_paths = os.path.join(subject_directory, 'sub-*', 'ses-1', 'func', '*desc-confounds_timeseries.tsv')
	tsv_files = glob.glob(tsv_files_paths)
	for tsv in tsv_files:
		for task in ['efnback1','efnback2','reward1','reward2']:
			if (task in tsv):
				subject_directory = os.path.dirname(os.path.dirname(os.path.dirname(tsv)))
				output_path = os.path.join(subject_directory, 'motion')
				os.path.mkdir(output_path)
				pull_motionpars(task, tsv, output_path) 
		
