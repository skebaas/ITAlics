# Neuroimaging  1st Levels Workflow using Nipype

This Python script is used to define a workflow for a neuroimaging task using the Nipype package in Python. The script takes three arguments: a configuration file specifying the details of the task, the name of the task you are using, and the directory where the data is stored.

The script sets up a Nipype workflow object with a base directory, and sets up the first level analysis node with the specified contrasts. If there are multiple runs, the script sets up merge points for the functional data, design matrices, movement data, and regressors.

For each run, the script sets up nodes for data processing and analysis, such as data smoothing and creating design matrices. If there is only one run, the script connects these components into a pipeline. If there are multiple runs, the script merges the data from each run before connecting to the first level analysis node.

The script also sets up a datasink node to store output files, and defines regions of interest for psychophysiological interaction (PPI) analysis. For each ROI, the script sets up a PPI node and an estimation of contrasts node, before outputting the contrast images to the datasink.

The output files to be saved in the datasink are defined in the configuration script.

## Dependencies

The script requires the following Python packages:

-   configparser
-   scipy.io
-   nipype
-   nipype.pipeline.engine
-   nipype.interfaces.utility
-   nipype.interfaces.fsl
-   nipype.interfaces.spm

The following software must be installed on the system:
- spm12
- fsl
- matlab/R2015

## Usage

The script can be run from the command line using the following command:


`python code/analysis/bin/task.py <path_to_config_file> <sequence_name> <path_to_subject_directory>`

Where:

-   `task.py`: Main python script to run
-   `<path_to_config_file>`: The path to the configuration file.
-   `<sequence_name>`: The name of the sequence for the neuroimaging task.
-   `<path_to_subject_directory>`: The path to the subject directory containing the data.

## Configuration File

The configuration file should be in the following format:


```python
[Default] 
data_dir = "/path/to/data/directory"
ROI_dir = "/path/to/ROI/directory" 
ROI_suffix = ".nii.gz" 
ROIS = {"ROI_1": "/path/to/ROI_1.nii.gz", "ROI_2": "/path/to/ROI_2.nii.gz"} output_name = "output_directory_name"  
[Task] 
sequence = "sequence_name"
output_name = "test_run"
runs = number_of_runs 
contrasts = {"contrast_1": "contrast_1_file", "contrast_2": "contrast_2_file", ...} 
ppi_contrasts = {"ppi_contrast_1": "ppi_contrast_1_file", "ppi_contrast_2": "ppi_contrast_2_file", ...} 
design_script = "/path/to/design_script.m" 
PPI_ROIS = [("ROI_name_1", "/path/to/ROI_1.nii.gz"), ("ROI_name_2", "/path/to/ROI_2.nii.gz"), ...]
```

Examples of the configuration file can be found in `code/analysis/bin/`
Where:

-   `data_dir`: The path to the directory containing the data.
-   `ROI_dir`: The path to the directory containing the ROI masks.
-   `ROI_suffix`: The suffix for the ROI mask file.
-   `ROIS`: A dictionary of ROI names and their corresponding mask file paths.
-   `output_name`: The name of the output directory.
-   `sequence`: The name of the sequence for the neuroimaging task.
-   `runs`: The number of runs for the neuroimaging task.
-   `contrasts`: A dictionary of contrasts and their corresponding contrast files.
-   `ppi_contrasts`: A dictionary of contrasts to be used on the PPIs
-   `design_script`: Path to the matlab design script you would like to use (task dependent)
-   `PPI_ROIS`: Rois that you want to use in the GPPI analysis
