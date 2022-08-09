#! /usr/bin/python3
import nibabel as nib
import numpy as np
import pandas as pd

from nilearn import image, masking
from nilearn.input_data import NiftiMasker, NiftiLabelsMasker
import os.path
import glob 

'''values to change'''
atlas = "/home/skebaas/Downloads/bnatlasviewer/content/bnatlas.nii.gz"
atlas_txt = "/home/skebaas/Downloads/bnatlasviewer/content/bnatlas.nii.txt"
fmriprep_dir = '/home/skebaas/Preprocessed/'
atlas_file = glob.glob(atlas)
fmri_files = sorted(glob.glob( '/home/skebaas/Preprocessed/sub-3054/ses-1/func/sub-3054_ses-1_task-sop_run-001_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'))
#fmri_files = atlas_file + fmri_files
print(fmri_files)
csv_output = 'sop_run_1.csv'

atlas_img = nib.load(atlas)
roi_indices, roi_counts = np.unique(atlas_img.get_data(),
return_counts=True)
roi_indices = roi_indices[:247]
roi_counts = roi_counts[:247]
#ROI names
f = open(atlas_txt, 'r+')
rois = [line.strip('\n') for line in f.readlines()]
#remove the empty string at the end of the list
rois.pop()
f.close()

#Dictionary that will be converted to csv
'''key = subject_id
    value = dictionary where key is roi_name and value is number of voxels
'''
subject_dict = {}

#Search through subject in the fmri_files directory
for file_img in fmri_files:
    sub_id = os.path.basename(file_img)[:8]
    print(file_img)
    roi_dict = {}
    #roi_dict['subject'] = sub_id
    t=0
    for img in image.iter_img(file_img):
    	# img is now an in-memory 3D img
	    #sub_mean = image.mean_img(img)
	    sub_mean_mask = image.math_img("1. * i.astype(bool)", i=img)
	    #print(sub_mean_mask)

	    masker = NiftiLabelsMasker(atlas_img)
	    voxel_ratios = masker.fit_transform([sub_mean_mask])
	    
	    #Get number of voxels present for each roi
	    n_voxels = voxel_ratios[0,:] * roi_counts[1:]
	    print(voxel_ratios)
	    #n_voxels = voxel_ratios * roi_counts
	    i=1
	    index_name = f"ts-{t}"
	    #Create the dictionary with the roi name and value 
	    for roi in n_voxels:
	        roi_dict[rois[i-1]] = roi
	        #print(f'roi {i} has {roi} voxels')
	        i+=1
	    subject_dict[index_name] = roi_dict
	    print(f"subject_dict for {sub_id} at {i} is {subject_dict[index_name]}")
	    t += 1

#Convert the dataframe to a csv 
#data = pd.DataFrame(subject_dict).T.reset_index().rename(columns={'index':'sub-id'})[rois]
data = pd.DataFrame(subject_dict).T
print(f"writing output to {csv_output} in pwd")
data.to_csv(csv_output)
