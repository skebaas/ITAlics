#!/usr/bin/python3

def f_exists(foldername):
	import os 
	import glob
	import shutil
	"""
	this function checks to see if a folder already exists.  If it does, it will output 'True'. Else it will be 'False'
	:param foldername: <string> path to folder
	"""
	if os.path.exists(foldername):
		print(f"WARNING: {foldername} already exists!")
		return True
	else:
		print(f"Creating {foldername}...")
		return False
def pull_motionpars(tsv, output_path):
	import os 
	import glob
	import shutil
	import pandas as pd
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

def make_folders(subject_path, output_path):
	import os 
	import glob
	import shutil

	subject = os.path.join(subject_path)
	output_path = output_path
	subject_ID = os.path.basename(subject_path)


	'''variable holding only the ID number.  Will be used later'''
	just_ID = subject_ID.replace('sub-', '')
	gold_format = f'diamond_1.{just_ID}'
	subj_dir_output = os.path.join(output_path, gold_format)
	if not f_exists(subj_dir_output):
		os.mkdir(subj_dir_output)

	'''hold files for each task in the following lists:'''
	funcs_path = os.path.join(subject, 'ses-1', 'func','*preproc_bold.nii.gz')
	masks_path = os.path.join(subject, 'ses-1', 'func', '*brain_mask.nii.gz')
	anats_path = os.path.join(subject, 'ses-1', 'anat', '*run-1_desc-preproc_T1w.nii.gz')
	'''!!! the folowing line is only searching for subjects in the BPD folder! This will eventually have to be changed!!!!'''
	txts1_path = os.path.join('/data', 'D2', 'Dicom', 'BPD', subject_ID, 'ses-1', 'BHV*', 'Scan', '*[1-2]*.txt')
	txts2_path = os.path.join('/data', 'D2', 'Dicom', 'NON', subject_ID, 'ses-1', 'BHV*', 'Scan', '*[1-2]*.txt')
	txts3_path = os.path.join('/data', 'D2', 'Dicom', 'BPD', subject_ID, 'ses-1', 'BHV*', 'Scan', 'subject*.txt')
	txts4_path = os.path.join('/data', 'D2', 'Dicom', 'BPD', subject_ID, 'ses-1', 'BHV*', 'Scan', 'subject*.txt')
	funcs = glob.glob(funcs_path)
	masks = glob.glob(masks_path)
	anats = glob.glob(anats_path)
	txts1 = glob.glob(txts1_path)
	txts2 = glob.glob(txts2_path)
	txts3 = glob.glob(txts3_path)
	txts4 = glob.glob(txts4_path)
	txts = txts1 + txts2 + txts3 + txts4

	'''Create anat folder for the subject'''
	anat_dest = os.path.join(subj_dir_output, 'anat')
	if not f_exists(anat_dest):
		os.mkdir(anat_dest)
	for anat in anats:
		output_anat = f"{anat_dest}/diamond_{just_ID}_anat.nii.gz"
		if not f_exists(output_anat):
			shutil.copy(anat, output_anat)


	'''following 2 lines are for testing purposes only'''
	#print(f" \t task files for {subject_ID} are {funcs}") 
	#print(f" \t mask files for {subject_ID} are {masks}")

	for task in tasks:
		task_output = os.path.join(subj_dir_output, task)
		#print(task)
		#print(task_output)
		try:
			os.mkdir(task_output)
		except FileExistsError:
			print(f'WARNING: {task_output} already exists!')
		#print(f" \t \t {task} folder created at {task_output}")

		for func in funcs:
			if task.replace('_', '') in func:
			#if ('reward1' in func) or ('reward2' in func) or ('efnback_1' in func) or ('efnback_2' in func) or ('dynface' in func):	
				#print(f" \t \t \t moving {func} to {task_output}")
				output_func = f"{task_output}/diamond_{just_ID}_{task}.nii.gz"
				if not f_exists(output_func):
					shutil.copy(func, output_func)
					'''consider adding function here to unzip files if errors occur with spm'''
				continue
			if 'dynface' in func and task == 'dynamic_faces':
				output_func = f"{task_output}/diamond_{just_ID}_dynamic_faces.nii.gz"
				if not f_exists(output_func):
					shutil.copy(func, output_func)

		for mask in masks:
			if task.replace('_', '') in mask:
				#print(f" \t \t \t moving {mask} to {task_output}")
				output_mask = f"{task_output}/diamond_{just_ID}_{task}_mask.nii.gz"
				if not f_exists(output_mask):
					shutil.copy(mask, output_mask)
					'''consider adding function here to unzip files if errors occur with spm'''
				continue
			if 'dynface' in mask and task == 'dynamic_faces':
				output_func = f"{task_output}/diamond_{just_ID}_dynamic_faces_mask.nii.gz"
				if not f_exists(output_func):
					shutil.copy(mask, output_func)

	'''We now have to pull the task.txt files from the DICOM folder'''
	for txt in txts:
		if ('Reward' in txt) and '1.txt' in txt:
			output_txt = f"{subj_dir_output}/reward_1/diamond_{just_ID}_reward_1_task.txt"
			#print(f" \t \t \t moving {txt} to {output_txt}")
			shutil.copy(txt, output_txt)
		if ('Reward' in txt) and '2.txt' in txt:
			output_txt = f"{subj_dir_output}/reward_2/diamond_{just_ID}_reward_2_task.txt"
			#print(f" \t \t \t moving {txt} to {output_txt}")
			shutil.copy(txt, output_txt)
		if ('EFNBACK' in txt) and '1.txt' in txt:
			output_txt = f"{subj_dir_output}/efnback_1/diamond_{just_ID}_efnback_1_task.txt"
			#print(f" \t \t \t moving {txt} to {output_txt}")
			shutil.copy(txt, output_txt)
		if ('EFNBACK' in txt) and '2.txt' in txt:
			output_txt = f"{subj_dir_output}/efnback_2/diamond_{just_ID}_efnback_2_task.txt"
			#print(f" \t \t \t moving {txt} to {output_txt}")
			shutil.copy(txt, output_txt)
		if ('subject' in txt) and '.txt' in txt:
			output_txt = f"{subj_dir_output}/dynamic_faces/diamond_{just_ID}_dynamic_faces_task.txt"
			#print(f" \t \t \t moving {txt} to {output_txt}")
			shutil.copy(txt, output_txt)

def get_motionpars(subject_path, output_dir):
	import pandas as pd
	pd.options.display.float_format = '{:.5f}'.format
	import glob
	import os
	import shutil

	subject = subject_path
	#tsv_files_paths = os.path.join('/data', 'D2', 'Preprocessed_OC', 'fmriprep', 'sub-*', 'ses-1', 'func', '*desc-confounds_timeseries.tsv')
	tsv_files_paths = os.path.join(subject, 'ses-1', 'func', '*desc-confounds_timeseries.tsv')
	#tsv_files_paths = os.path.join(subject, 'ses-1', 'func', '*desc-confounds_regressors.tsv')
	#l1_subjects_path = os.path.join('/data', 'D2', '1stlevelanalysis', 'diamond_1*')
	l1_subjects_path = os.path.join(output_dir, 'diamond_1*')
	tsv_files = glob.glob(tsv_files_paths)
	l1_subjects = glob.glob(l1_subjects_path)
	print(tsv_files_paths)
	subjName = os.path.basename(subject)
	subjID = subjName.split('-')[-1]
	gold_format = f"diamond_1.{subjID}"
	output_path = os.path.join(output_dir, gold_format)

	for tsv in tsv_files:
		if ('efnback1' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'efnback_1')
			if not f_exists(f"{output_folder}/motion.1D"):
				pull_motionpars(tsv, output_folder)
		if ('efnback2' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'efnback_2')
			if not f_exists(f"{output_folder}/motion.1D"):
				pull_motionpars(tsv, output_folder)
		if ('reward1' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'reward_1')
			if not f_exists(f"{output_folder}/motion.1D"):
				pull_motionpars(tsv, output_folder)
		if ('reward2' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'reward_2')
			if not f_exists(f"{output_folder}/motion.1D"):
				pull_motionpars(tsv, output_folder)
		if ('dynface' in tsv) and (subjName in tsv):
			output_folder = os.path.join(output_path, 'dynamic_faces')
			if not f_exists(f"{output_folder}/motion.1D"):
				pull_motionpars(tsv, output_folder)


if __name__ == "__main__":	
	import os 
	import glob
	import shutil
	import sys

	'''This script requires the input subject directory to be the output folder of fmriprep containing the Preprocessed subjects'''
	#ex: /data/D2/Preprocessed_OC/fmriprep/
	subject_path = sys.argv[1]
	#ex: /data/D2/1stlevelAnalysis/
	output_path = sys.argv[2]
	tasks = ['efnback_1', 'efnback_2', 'reward_1', 'reward_2', 'dynamic_faces']

	if not f_exists(output_path):
		os.mkdir(output_path)
	preproc_subject = glob.glob(subject_path)
	

	for subject in list(preproc_subject):
		make_folders(subject, output_path)
		get_motionpars(subject, output_path)

	'''
	NON_tasks_path = os.path.join('/data', 'D2', 'Dicom', 'NON', 'sub-*', 'ses-1', 'BHV*', 'Scan', '*[1-2]*.txt')
	NON_tasks = glob.glob(NON_tasks_path)
	'''
	#make_folders()
	#get_motionpars(inputPath, output_path)
