#! /usr/bin/python3

import json
import os
import glob
import sys, getopt
import pandas as pd 


def main(argv):
    # Define input and output file variables with default values
    inputfile = ''
    outputfile = 'Subjects_default'

    # Try to get input and output file options from command line arguments
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print("parse_json.py -i <subject directory> -o <outputfile>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("parse_json.py -i <subject directory> -o <outputfile>")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = os.path.abspath(arg)
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    # Create a DataFrame with specified columns
    df = pd.DataFrame(columns=["File_Name", "Coil_Name", "Phase_Direction", "Echo_Time"])

    # Get list of all subjects in the input directory
    subjects_path = os.path.join(inputfile, 'sub-*')
    subjects = glob.glob(subjects_path)

    # Iterate through each subject
    for subject in subjects:
        # Get list of all json files in the subject's session 1 functional directory
        files_path = os.path.join(subject, 'ses-1/func/*json')
        files = glob.glob(files_path)

        # Iterate through each json file
        for file in files:
            filename = os.path.basename(file)
            json = search_json(file)
            df = add_row(filename, json, df)

    # Output the DataFrame to a csv file
    print(f"Outputting data to {outputfile}.csv")
    df.to_csv(f'{outputfile}.csv')

def search_json(inputfile: str) -> dict:
    """
    Returns a json as a dictionary
    :param inputfile: Path to json file
    :return: json_file or -1 if there is a KeyError
    """
    try:
        with open(inputfile, 'r') as f:
            json_file = json.load(f)
            return json_file
    except KeyError:
        return -1

def add_row(filename: str, json_file: dict, df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a new dataframe after adding a row to the input dataframe 'df'
    :param filename: Path to json file
    :param json_file: Dictionary containing json information
    :param df: Pandas dataframe which the row will be appended to
    :return: df with a new row appended
    """
    names = filename.split("_")
    name = f"{names[0]}_{names[2]}_{names[-2]}"
    row = pd.Series(data={'File_Name': name, 'Coil_Name': json_file["ReceiveCoilName"], 'Phase_Direction': json_file["PhaseEncodingDirection"], 'Echo_Time': json_file["EchoTime"]})
    df = df.append(row, ignore_index=True)
    return df

if __name__ == "__main__":
	main(sys.argv[1:])
