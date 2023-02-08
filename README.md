README
Welcome to the repository for the neuroimaging pipeline. This repository contains several subdirectories, each with its own specific purpose.

1stslevels
The 1stslevels subfolder contains scripts for 1st-level processing of the Diamond1, Diamond2, and TBD datasets. These scripts are compatible with tasks such as reward, dynamic faces, and feedback. The scripts in this subfolder have been repurposed from the old "Gold" pipeline and the scripts that can be directly used are located in the bin subfolder of 1stslevels.

ICA_FMRIPREP
The ICA_FMRIPREP subfolder contains changes that have been made to the fmriprep20.2.6 version. These files need to replace the corresponding files found in fmriprep.

Physio
The Physio subfolder contains scripts for processing physiologic data in the first-level analyses.

QC
The QC subfolder contains all scripts related to quality control for raw data, pre-preprocessed data, and general motion quality control.

Scripts_Pre_PreProcessing
The Scripts_Pre_PreProcessing subfolder contains scripts used before running through fmriprep. This includes scripts for adding fieldmaps. This folder may be deprecated soon.

If you have any questions or concerns, please consult the documentation within each subfolder or reach out to the repository maintainer.
