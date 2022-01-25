#! /usr/bin/env python3

"""The heudiconv_unstable.simg script generates a hidden .heudiconv folder which contains valueable dicom_info 
that we use for QC purposes.  This script can output the paths of these hidden .heudiconv folders and copy dicom_info.tsv files to a specified output directory"""
###@author Alexander Skeba###

import argparse
import os
import glob
import shutil

def run_script(script, stdin=None):
    """Returns (stdout, stderr), raises error on non-zero return code"""
    import subprocess
    # Note: by using a list here (['bash', ...]) you avoid quoting issues, as the
    # arguments are passed in exactly this order (spaces, quotes, and newlines won't
    # cause problems):
    proc = subprocess.Popen(['bash', '-c', script],
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode:
        raise ScriptException(proc.returncode, stdout, stderr, script)
    return stdout, stderr

class ScriptException(Exception):
    def __init__(self, returncode, stdout, stderr, script):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        Exception().__init__('Error in script')

def find_heudiconv(dir_path):
    if not os.path.isdir(dir_path):
        print(f"ERROR: {dir_path} is not a valid directory")
        exit(1)
    bash_script = f'find {dir_path} -type d -name ".heudiconv" '
    print(f"Outputting paths containing '.heudiconv' folders in {dir_path} ...")
    run_script(bash_script)
    print("Exiting...")
    exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--search", help="find hidden .heudiconv folders in folder and output their paths",type=str)
    parser.add_argument("-i", "--in_dir", help="path to .heudiconv foler", type=str)
    parser.add_argument("-o","--out_dir", help="path to where you want files to be copied",type=str)
    args = parser.parse_args()
    
    if args.search:
        search_dir = args.search
        find_heudiconv(search_dir)
    if args.in_dir:
        input_dir = args.in_dir
        if not os.path.isdir(input_dir):
            print(f"ERROR: {input_dir} is not a valid directory")
            exit(1)
        input_dir = os.path.abspath(input_dir)
    else:
        print("ERROR:  Input Directory MUST be specified")
        parser.print_help()
        exit(1)
    if args.out_dir:
        output_dir = args.out_dir
        if not os.path.isdir(output_dir):
            print(f"{output_dir} does not exists. Creating it now...")
            os.mkdir(output_dir)
        output_dir = os.path.abspath(output_dir)
    else:
        print("ERROR: Output Directory MUST be specified")
        parser.print_help()

    subject_paths = glob.glob(input_dir + '/*')
    for subject in subject_paths:
        subject_name = os.path.basename(subject)
        session_paths = glob.glob(subject+'/ses-*')

        for session in session_paths:
            dicom_info_files = os.path.join(session,'info','dicominfo_*.tsv')
            dicom_info_path = glob.glob(dicom_info_files)[0]
            dicom_info_filename = os.path.basename(dicom_info_path)
            output_filename = f'{subject_name}_{dicom_info_filename}'
            output_file_path = os.path.join(output_dir, output_filename)
            if os.path.exists(output_file_path):
                print(f"WARNING: {output_file_path} already exists! Skipping...")
            else:
                print(f"copying file --> {dicom_info_path} to {output_file_path}...")
                shutil.copy(dicom_info_path, output_file_path)

