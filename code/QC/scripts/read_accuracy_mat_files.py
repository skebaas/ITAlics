#! /usr/bin/env python3
import sys
import glob
import scipy.io
import pandas as pd
import os

def read_accuracy_data(file_pattern):
    data = {'ID': [], 'Accuracy': []}  # dictionary to store data
    
    # loop over matching .mat files
    for i, file_name in enumerate(glob.glob(file_pattern)):
        # load .mat file
        mat_file = scipy.io.loadmat(file_name)
        
        # extract accuracy information
        accuracy = mat_file['accuracy'][0][0][0]
        
        # add data to dictionary
        data['ID'].append(i+1)
        data['Accuracy'].append(accuracy)
    
    # create pandas dataframe
    df = pd.DataFrame(data)
    
    return df

# example usage
if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] in ['-h', '--help']:
        print(f'Usage: {sys.argv[0]} <file_pattern>')
        print(f'Example {sys.argv[0]} acc_efnback1_*.mat')
        sys.exit()
    
    file_pattern = sys.argv[1]
    file_patter = os.path.abspath(file_pattern)
    df = read_accuracy_data(file_pattern)
    print(df)

