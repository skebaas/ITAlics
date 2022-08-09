#! /usr/bin/python3
import nibabel as nib
import numpy as np
import pandas as pd

from nilearn import image, masking
from nilearn.input_data import NiftiMasker, NiftiLabelsMasker
import os.path
import glob 

'''values to change'''
atlas = "/home/schumerm/percent_coverage/bnatlas.nii.gz"
atlas_txt = "/home/schumerm/percent_coverage/bnatlas.nii.txt"
fmriprep_dir = '/data/D2/data/Preprocessed/NON/fmriprep/'
atlas_file = glob.glob(atlas)
fmri_files = sorted(glob.glob(fmriprep_dir + 'sub*/ses-1/func' + '/*task-dynface*preproc_bold.nii'))
fmri_files = atlas_file + fmri_files
csv_output = 'dynface_NON_run_1.csv'

print(fmri_files)
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
signal_dict = {}
#Search through subject in the fmri_files directory
for subject in fmri_files:
    sub_id = os.path.basename(subject)[:9]
    print(sub_id)
    roi_dict = {}
    roi_signal_dict = {}
    #roi_dict['subject'] = sub_id
    sub_mean = image.mean_img(subject)

    sub_mean_mask = image.math_img("1. * i.astype(bool)", i=sub_mean)
    masker = NiftiLabelsMasker(atlas_img)
    voxel_ratios = masker.fit_transform([sub_mean_mask])
    #Get total signal strength in each ROI and output to np.Array form
    signal_strength = masker.fit_transform([sub_mean])
    
    #Get number of voxels present for each roi
    n_voxels = voxel_ratios[0,:] * roi_counts[1:]
    avg_signal_strength = signal_strength[0]/n_voxels
    i=1
    #Create the dictionary with the roi name and value 
    for roi in n_voxels:
        roi_dict[rois[i-1]] = roi
        #print(f'roi {i} has {roi} voxels')
        i+=1
    subject_dict[sub_id] = roi_dict
    print(f"subject_dict for {sub_id} is {subject_dict[sub_id]}")
    j=1
    for roi in avg_signal_strength:
        roi_signal_dict[rois[j-1]] = roi
        j+=1
    signal_dict[sub_id] = roi_signal_dict
#Convert the dataframe to a csv 
#data = pd.DataFrame(subject_dict).T.reset_index().rename(columns={'index':'sub-id'})[rois]
data = pd.DataFrame(subject_dict).T
data2 = pd.DataFrame(signal_dict).T
print(f"writing output to {csv_output} in pwd")
print(f"writing output to sig_strength{csv_output} in pwd")
data.to_csv(csv_output)
data.to_csv(f"sig_strength{csv_output}")
