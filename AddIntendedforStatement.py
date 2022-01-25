#!/usr/bin/python3



def addIntendedforStatement(subject):
	subject = subject
	import os
	import sys
	import glob
	import json

	subject_name = os.path.basename(subject)
	print(f"adding 'IntendedFor' Statement to fieldmaps for {subject_name}")

	#Edit 'ses-1' to be the same as your session
	fmapsPath = os.path.join(subject, 'ses-1', 'fmap', '*.json')
	fmaps = glob.glob(fmapsPath)
	funcsPath = os.path.join(subject, 'ses-1', 'func', '*.nii*')
	funcs = glob.glob(funcsPath)


	#substring to be removed from absolute path of functional files
	pathToRemove = subject + '/'
	funcs = list(map(lambda x: x.replace(pathToRemove, ''), funcs))
	for fmap in fmaps:
		with open(fmap, 'r') as data_file:
			fmap_json = json.load(data_file)
		fmap_json['IntendedFor'] = funcs

		with open(fmap, 'w') as data_file:
			fmap_json = json.dump(fmap_json, data_file, indent=4)

if __name__ == "__main__":	
	import sys
	import os
	import glob
	import json

	subjectPath = sys.argv[1]
	subject = glob.glob(subjectPath)
	addIntendedforStatement(subject[0])
	
