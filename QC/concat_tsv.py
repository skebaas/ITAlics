''' A simple program to convert a directory of TSVs Summarizing Dicom Files for a Dataset and combining them into a single Pandas DataFrame'''


from re import sub
import pandas as pd
import glob
import os
import sys

from natsort import natsorted


def main():
    arguments = sys.argv

    if(len(arguments) != 3): #Simple Argument Check, can make more comprehensive later
        print('Please pass a single directory containing TSVs for your first argument, and an output filepath as your second argument')

    print("Path to TSV Dir: " + arguments[1])
    print("Path to output file: " + arguments[2])


    tsv_dir = arguments[1] #Should be an absolute path to your directory of choice
    tsv_list = glob.glob(tsv_dir + '/*.tsv') #Can easily change to .csv

    sorted_list = natsorted(tsv_list)
    list_len = len(sorted_list) + 1

    print(f'Found {list_len} TSV Files')
    print(sorted_list)


    print("Iterating through list...")
    dataset_df = pd.DataFrame() # Creating empty df to store our compiled dataset
    for x in sorted_list:
        sub_basename = os.path.basename(x)
        print(sub_basename) #Tester Code

        sub_basename_split = sub_basename.split('_')

        sub_num = sub_basename_split[0]
        ses_num = sub_basename_split[2].split('.')[0].split('-')[1]

        print(sub_num) #Tester Code
        print(ses_num) #Tester Code

        sub_df = pd.read_csv(x, sep='\t') # Remove Sep parameter to work with normal comma separated files

        sub_df['subject_num'] = sub_num # All entries from this tsv should have the same subject number
        sub_df['session_num'] = ses_num # All entries from this tsv should have the same session number

        first_col = sub_df.pop('subject_num') # Moving the Subject and Session Columns to the front
        second_col = sub_df.pop('session_num')
        sub_df.insert(0, 'session_num', second_col)
        sub_df.insert(0, 'subject_num', first_col)
    
        dataset_df = dataset_df.append(sub_df) #Appending onto the whole dataframe
    
    dataset_df.to_csv(path_or_buf=arguments[2], sep='\t')



    

if __name__ == '__main__':
    main()