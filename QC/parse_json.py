#! /usr/bin/python3

import json
import os
import glob
import sys, getopt
import pandas as pd 

def main(argv):
	inputfile = ''
	outputfile = 'Subjects_default'

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

	df = pd.DataFrame(columns=["File_Name", "Coil_Name", "Phase_Direction", "Echo_Time"])
	subjects_path = os.path.join(inputfile, 'sub-*')
	subjects = glob.glob(subjects_path)
	for subject in subjects:
		files_path = os.path.join(subject, 'ses-1/func/*json')
		files = glob.glob(files_path)
		for file in files:
			filename = os.path.basename(file)
			json = search_json(file)
			df = add_row(filename, json, df)
	print(f"Outputting data to {outputfile}.csv")
	df.to_csv(f'{outputfile}.csv')
def search_json(inputfile):
	'''
	Returns a json as a dictionary

	Parameters:
		inputfile (str): Path to json file

	Returns:
		json_file (obj): 
	'''
	try:
		with open(inputfile, 'r') as f:
			json_file = json.load(f)
			return json_file
	except KeyError:
		return -1

def add_row(filename, json_file, df):
	'''
	Returns a new dataframe after adding a row to the input dataframe 'df'

	Parameters:
		filename (str): Path to json file
		json_file (str): Dictionary containing json information
		df (DataFrame): Pandas dataframe which the row will be appended to 
	'''
	names = filename.split("_")
	name = f"{names[0]}_{names[2]}_{names[-2]}"
	row = pd.Series(data={'File_Name': name, 'Coil_Name': json_file["ReceiveCoilName"], 'Phase_Direction': json_file["PhaseEncodingDirection"], 'Echo_Time': json_file["EchoTime"]})
	df = df.append(row, ignore_index=True)
	return df

if __name__ == "__main__":
	main(sys.argv[1:])
