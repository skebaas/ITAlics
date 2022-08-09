#! /usr/bin/python3

import pandas as pd
pd.options.display.float_format = '{:.5f}'.format
import glob
import os
import shutil



def pull_motionpars(tsv, output_path):
	"""
	this function will pull the motionparameters the tsv passed to it and output a motion.1D 
	file to the corresponding task folder in /data/D2/1stlevelanalysis
	:param tsv: <string> confounds_timeseries.tsv containing motion parameters
	:param output_path: <string> output path for motion.1D
	"""
	motion1D_path = os.path.join(output_path, 'motion.1D')
	confounds_tsv = pd.read_csv(tsv, sep='\t')
	motion_pars = confounds_tsv[['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']]
	#for col in motion_pars.iteritems():
	#	#col = pd.to_numeric(col, downcast="float")
	#	print(col)
	motion_pars.to_csv(motion1D_path, header=None, index=None, sep='\t')
	print(f"\t created {motion1D_path} from {tsv}")



#path to subjects
subjects_path = os.path.join('/data', 'D2', 'Preprocessed_OC', 'fmriprep', 'sub-*')
tsv_files_paths = os.path.join('/data', 'D2', 'Preprocessed_OC', 'fmriprep', 'sub-*', 'ses-1', 'func', '*desc-confounds_timeseries.tsv')
l1_subjects_path = os.path.join('/data', 'D2', '1stlevelanalysis', 'diamond_1*')
tsv_files = glob.glob(tsv_files_paths)
subjects = glob.glob(subjects_path)
l1_subjects = glob.glob(l1_subjects_path)

print(tsv_files_paths)
for subject in subjects:
	subjName = os.path.basename(subject)
	subjID = subjName.split('-')[-1]
	gold_format = f"diamond_1.{subjID}"
	output_path = os.path.join('/data', 'D2', '1stlevelanalysis', gold_format)
	print(subjID)
	for tsv in tsv_files:
		if ('efnback1' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'efnback_1')
			pull_motionpars(tsv, output_folder)
		if ('efnback2' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'efnback_2')
			pull_motionpars(tsv, output_folder)
		if ('reward1' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'reward_1')
			pull_motionpars(tsv, output_folder)
		if ('reward2' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'reward_2')
			pull_motionpars(tsv, output_folder)


	
