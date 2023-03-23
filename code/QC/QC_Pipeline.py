#!/usr/bin/env python3.8
import os
import shutil
import sys
import pandas as pd
import glob
class Task_DataFrame:
    """
    Initialize the class with the given task and fMRI preprocessing directory.

    This method initializes the class with a specified task and the directory
    containing fMRI preprocessing outputs. It computes various dataframes
    related to mean and max framewise displacement and translation values, as well
    as the number of excluded mean and max values. It also initializes variables
    for percent coverage tables.

    Parameters
    ----------
    task : str
        The task identifier (e.g., 'rest', 'emotion', etc.).
    fmriprep_dir : str
        The path to the directory containing the fMRI preprocessing outputs.

    Attributes
    ----------
    task : str
        The task identifier.
    fmriprep_dir : str
        The path to the directory containing the fMRI preprocessing outputs.
    task_df : pd.DataFrame
        The dataframe containing mean and max values for the given task.
    mean_highlighted : pd.DataFrame
        Mean values flagged for exclusion (greater than 0.7).
    max_highlighted : pd.DataFrame
        Max values flagged for exclusion (greater than 5).
    mean_table : pd.DataFrame
        Table of mean values filtered based on specified thresholds.
    max_table : pd.DataFrame
        Table of max values filtered based on specified thresholds.
    max_table_highlighted : pd.DataFrame
        Highlighted max values in the max_table flagged for exclusion.
    mean_table_highlighted : pd.DataFrame
        Highlighted mean values in the mean_table flagged for exclusion.
    mean_excluded : int
        The number of mean values excluded based on the specified thresholds.
    max_excluded : int
        The number of max values excluded based on the specified thresholds.
    percent_coverage_table : str
        An empty string to store the percent coverage table information later.
    percent_coverage_highlighted : str
        An empty string to store the highlighted percent coverage information later.

    Example
    -------
    >>> my_data_frame = Task_DataFrame("rest", "/path/to/fmriprep/outputs")
    """

    def __init__(self, task, fmriprep_dir):
        self.task = task
        self.fmriprep_dir = fmriprep_dir
        task_name = task
        print(f"{task_name} and {fmriprep_dir}")
        task_df = self.generate_mean_max(fmriprep_dir, task_name)
        self.task_df = task_df
        self.mean_highlighted = self.flag_exclusion_val(task_df, 0.7, 'mean') #flag mean values greater than 0.7
        self.max_highlighted = self.flag_exclusion_val(task_df, 5, 'max') #flag max values greater than 5
        self.mean_table = self.task_df[(task_df['framewise_displacement_max'] > 5) | (task_df['trans_z_max'] > 5) | (task_df['trans_y_max'] > 5) | (task_df['trans_x_max'] > 5)]
        self.max_table = self.task_df[(task_df['framewise_displacement_mean'] > 0.5) | (task_df['trans_z_mean'] > 0.5) | (task_df['trans_y_mean'] > 0.5) | (task_df['trans_x_mean'] > 5)]
        self.max_table_highlighted = self.flag_exclusion_val(task_df[(task_df['framewise_displacement_max'] > 5) | (task_df['trans_z_max'] > 5) | (task_df['trans_y_max'] > 5) | (task_df['trans_x_max'] > 5)], 5, 'max')
        self.mean_table_highlighted = self.flag_exclusion_val(task_df[(task_df['framewise_displacement_mean'] > 0.5) | (task_df['trans_z_mean'] > 0.5) | (task_df['trans_y_mean'] > 0.5) | (task_df['trans_x_mean'] > 0.5)], 0.5, 'mean')
        self.mean_excluded = self.mean_table['framewise_displacement_mean'].count()
        self.max_excluded = self.max_table['framewise_displacement_max'].count()
        self.percent_coverage_table = ''
        self.percent_coverage_highlighted = ''

    def make_directory(self, folder):
        import os
        if os.path.exists(folder):
            return True
        else:
            os.mkdir(folder)

    def get_file_info(self, file):
        """Takes in path to confound regressors file and returns subject name and task name"""
        import os
        sub_name = None
        ses_name = 'ses-1'
        task_name = None
        for line in os.path.basename(file).split('_'):
            if 'sub' in line:
                sub_name = line
            elif 'ses' in line:
                ses_name = line
            elif 'task' in line:
                task_name = line
            else:
                continue
        return sub_name, ses_name, task_name

    def pull_confounds(self, tsv, output_path, confounds=['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'framewise_displacement']):
            """
            this function will pull the regressors specified in 'confounds' from 'tsv' and output a motion.tsv to 'output_path'
            Return pandas dataframe containing confounds
            :param tsv: <string> confounds_timeseries.tsv containing motion parameters
            :param output_path: <string> output path for motion.tsv
            :param confounds: <list> regressors to be used in dataframe
            """
            import pandas as pd
            import os

            sub_name, ses_name, task_name = self.get_file_info(tsv)
            motion_file = os.path.join(output_path, f'{sub_name}_{ses_name}_{task_name}_motion.tsv')
            confounds_file = pd.read_csv(tsv, sep='\t')
            confounds = confounds_file[confounds]
            confounds.to_csv(motion_file, index=None, sep='\t')

            return confounds

    def generate_metrics(self, tsv, confounds_df, output_path):
        """
        This function will generate the following metrics for the confounds_df columns:
        count/mean/std/min/25%/50%/75%/max
        
        Uses 'tsv' to extract information about the name of the file
        
        """
        sub_name, ses_name, task_name = self.get_file_info(tsv)
        confound_metrics = confounds_df.describe()
        confound_metrics_file = os.path.join(output_path, f'{sub_name}_{ses_name}_{task_name}_metrics.tsv')
        confound_metrics.to_csv(confound_metrics_file, sep='\t')
        
        return confound_metrics

    def generate_mean_max(self, fmriprep_dir, task):
        fmriprep_path = os.path.abspath(fmriprep_dir)
        #print(fmriprep_path)
        confounds_files = glob.glob(os.path.join(fmriprep_path, 'sub-*','ses-*','func',f'*{task}*confounds_timeseries.tsv'))
        #print(confounds_files)
        subject_dict = {}
        for file in confounds_files:
            confound_dict = {}
            sub_name, ses_name, task_name = self.get_file_info(file)
            out_dir = f'{fmriprep_path}/QC_output'
            self.make_directory(out_dir)
            confounds = self.pull_confounds(file,out_dir)
            confounds_metrics = self.generate_metrics(file, confounds, out_dir)
            for col in confounds_metrics:
                confound_dict[f'{col}_mean'] = confounds_metrics[col]['mean']
                confound_dict[f'{col}_max'] = confounds_metrics[col]['max']
            subject_dict[sub_name] = confound_dict
        data = pd.DataFrame(subject_dict).T
        data.sort_index(inplace=True)
        return data

    def flag_exclusion_val(self, df, exclusion_val, metric='mean'):
        """
        This function takes in the output from 'percent_coverage()' and searches for values less than 'exclusion_val'
        and highlights the values.
        :param df: <DataFrame> Percent Coverage DataFrame
        :param exlusion_val: <float> Value to search 'df' for
        :param metric: <string> Either 'mean' or 'max'
        """
        def highlight_cols(s):
            color = ''
            if s > exclusion_val:
                color = 'red'
            return 'background-color: % s' % color

        metric_cols = []
        for col in df.columns:
            if metric in col:
                metric_cols.append(col)
        return df[metric_cols].style.applymap(highlight_cols)

    def percent_coverage(self, atlas, atlas_txt):
        fmriprep_dir = self.fmriprep_dir
        task_name = self.task
        """
        This function will return a pandas dataframe which contains ROIs specified by 'atlas'.
        The dataframe will show the percentage of voxels found compared to the maximum possible count found in each ROI.
        This will automatically be done on all subjects found in 'fmriprep_dir' with the specified 'task_name'
        
        :param atlas: <string> path to atlas file
        :param atlas_txt: <string> path to atlas file containing ROI names
        :param fmriprep_dir: <string> path to directory containing fmriprep pre-processed images
        :param task_name: <string> task that you wish to generate percent_coverage output for i.e. ('dynface, reward1, efnback2, ...')
        """
        
        import nibabel as nib
        import numpy as np
        import pandas as pd
        from nilearn import image, masking
        from nilearn.input_data import NiftiMasker, NiftiLabelsMasker
        import os.path
        import glob
        
        atlas_file = glob.glob(atlas)
        print(fmriprep_dir)
        fmri_files_path = os.path.join(fmriprep_dir, 'sub-*', 'ses-1', 'func',f'*{task_name}*preproc_bold.nii*')
        fmri_files = sorted(glob.glob(fmri_files_path))
        fmri_files = atlas_file + fmri_files
        print(fmri_files)
        csv_output = f'{fmriprep_dir}/QC_output/{task_name}_percentcoverage.tsv'
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
            #print(sub_id)
            roi_dict = {}
            roi_signal_dict = {}
            #roi_dict['subject'] = sub_id
            try:
                sub_mean = image.mean_img(subject)
            except:
                continue

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
            #print(f"subject_dict for {sub_id} is {subject_dict[sub_id]}")
            j=1
            for roi in avg_signal_strength:
                roi_signal_dict[rois[j-1]] = roi
                j+=1
            signal_dict[sub_id] = roi_signal_dict
        #Convert the dataframe to a csv 
        #data = pd.DataFrame(subject_dict).T.reset_index().rename(columns={'index':'sub-id'})[rois]
        data = pd.DataFrame(subject_dict).T
        data = data.div(data.iloc[0])
        data.to_csv(csv_output, sep='\t')
        
        self.percent_coverage_table = data

    def output_low_coverage(self, exclusion_val):
        """
        This function takes in the output from 'percent_coverage()' and searches for values less than 'exclusion_val'
        and highlights the values.
        :param pc_df: <DataFrame> Percent Coverage DataFrame
        :param exlusion_val: <float> Value to search 'pc_df' for
        """
        def highlight_cols(s):
            color = ''
            if s < exclusion_val:
                color = 'red'
            return 'background-color: % s' % color
        self.percent_coverage_highlighted = self.percent_coverage_table.style.applymap(highlight_cols)
        return self.percent_coverage_table.style.applymap(highlight_cols)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Create an instance of Task_DataFrame with the given task and fmriprep directory.')

    parser.add_argument('task', type=str, help='The task identifier (e.g., "rest", "efnback1", "reward2", "dynface", etc.)')
    parser.add_argument('fmriprep_dir', type=str, help='The path to the directory containing the fMRI preprocessing outputs.')

    args = parser.parse_args()
    print(args.fmriprep_dir)
    df_class = Task_DataFrame(args.task, args.fmriprep_dir)
    df_class.task_df.to_csv(f'{args.task}_motion.csv')
