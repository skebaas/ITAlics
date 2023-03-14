#! /usr/bin/env python3
import sys
import glob
import scipy.io
import pandas as pd
import os

def read_accuracy_data(files):
    data = {'ID': [], 'Run': [],'Accuracy': [], 'Average': []}  # dictionary to store data
    print(files)
    # loop over matching .mat files
    for i, file_name in enumerate(files):
        # load .mat file
        mat_file = scipy.io.loadmat(file_name)
        id_name = os.path.basename(file_name).split('.')[1]
        run_name = os.path.basename(file_name).split('_')[1]
        # extract accuracy information
        accuracy = mat_file['accuracy'][0][0][0]
        print(type(accuracy))
        # add data to dictionary
        data['ID'].append(id_name)
        data['Run'].append(run_name)
        data['Accuracy'].append(accuracy)
        data['Average'].append(accuracy.mean())
    
    # create pandas dataframe
    df = pd.DataFrame(data)
    
    return df

# example usage
if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] in ['-h', '--help']:
        print(f'Usage: {sys.argv[0]} <file_pattern> <output_name>')
        print(f'Example {sys.argv[0]} acc_efnback1_*.mat accuracies_efnback.txt')
        sys.exit()
    
    file_pattern = str(sys.argv[1])
    output_name = str(sys.argv[2])
    #file_pattern = os.path.abspath(file_pattern)
    files = glob.glob(file_pattern)
    df = read_accuracy_data(files)
    df.set_index('ID', inplace=True)
    df.sort_index(inplace=True)
    df.to_csv(output_name)
    print(df)

