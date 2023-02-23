import argparse
import pandas as pd
import numpy as np
import sys
import os
import subprocess
import glob
import shutil

def run_physio(physio_path, spm_path, working_dir, task_name):
    matlab_task = f'{working_dir}/PhysIO_template_{task_name}.m'
    subprocess.run(['matlab', '-nodesktop', '-nodisplay', '-nojvm', '-nosplash', '-r', 
    f"oldpath = path; path(\'{physio_path}\',oldpath); oldpath = path; path(\'{spm_path}\', oldpath);run(\'{matlab_task}\'); exit;"])
parser = argparse.ArgumentParser(description='Outputs a regression file to be input into the first level regressors for physiological noise correction employing Tapas.')
parser.add_argument('-s', metavar='Subject', help='Subject Dicom directory')
parser.add_argument('-w', metavar='Working Directory', help='Specify where you want physio derivatives to go')
args = parser.parse_args(sys.argv[1:])
subject_dir = str(args.s)
working_dir = str(args.w)
#subject_dir = '/data/D2/data/Dicom/BPD/sub-50003'
#working_dir = '/data/D2/data/jupyter-notebooks'
subject_name = os.path.basename(subject_dir)
working_dir = os.path.join(working_dir, subject_name)

#Environment Variables 
matlab_path = subprocess.check_output(args=['which', '-a', 'matlab'], encoding='utf-8')
matlab_path = matlab_path.rstrip()
matlab_path = os.path.realpath(matlab_path)

matlab_master_dir = os.path.dirname(matlab_path)
matlab_master_dir = os.path.dirname(matlab_master_dir)

physio_path = matlab_master_dir + '/toolbox/tapas-master/PhysIO/code'
spm_path = '/usr/local/software/spm12'
task_array = ["dynamicfaces", "efnback1", "efnback2", "reward1", "reward2"]

pulse_file = glob.glob(os.path.join(subject_dir,'ses-1','BHV*',"*.puls"))
resp_file = glob.glob(os.path.join(subject_dir,'ses-1','BHV*',"*.resp"))
dyn_dicoms = glob.glob(os.path.join(subject_dir,'ses-1','dynamic-faces_864*','*'))
efn1_dicoms = glob.glob(os.path.join(subject_dir,'ses-1','EFNback_1_864*', '*'))
efn2_dicoms = glob.glob(os.path.join(subject_dir,'ses-1','EFNback_2_864*', '*'))
rest_dicoms = glob.glob(os.path.join(subject_dir,'ses-1','resting-state_864*', '*'))
rew1_dicoms = glob.glob(os.path.join(subject_dir,'ses-1','reward_1_864*', '*'))
rew2_dicoms = glob.glob(os.path.join(subject_dir,'ses-1','reward_2_864*','*'))

if not os.path.exists(working_dir):
    os.makedirs(working_dir)
shutil.copyfile(pulse_file[0], f'{working_dir}/pulse.puls')
shutil.copyfile(resp_file[0], f'{working_dir}/resp.resp')
shutil.copyfile(dyn_dicoms[0], f'{working_dir}/dynamicfaces.dcm')
shutil.copyfile(efn1_dicoms[0], f'{working_dir}/efnback1.dcm')
shutil.copyfile(efn2_dicoms[0], f'{working_dir}/efnback2.dcm')
shutil.copyfile(rest_dicoms[0], f'{working_dir}/rest.dcm')
shutil.copyfile(rew1_dicoms[0], f'{working_dir}/reward1.dcm')
shutil.copyfile(rew2_dicoms[0], f'{working_dir}/reward2.dcm')

script_dir = os.path.dirname(os.path.abspath(__file__))
print(script_dir)
matlab_scripts = glob.glob(f'{script_dir}/*.m')
for script in matlab_scripts:
    shutil.copyfile(script, f'{working_dir}/{os.path.basename(script)}')

#Run the Physio Script and convert output to tsv
os.chdir(working_dir)
for x in task_array:
    print(f'Running {x}')
    run_physio(physio_path, spm_path, working_dir, x)
    df = pd.read_csv(f'multiple_regressors_{x}.txt', sep='\t', header=None)
    df.to_csv(f'multiple_regressors_{x}.tsv', sep='\t', header=False, index=False)

#Run basic physio qc for each subject
qc_files = glob.glob(f'{working_dir}/multiple_regressors_*.tsv')
for x in qc_files:
    subprocess.run([f'{script_dir}/physio_qc.py', x])
