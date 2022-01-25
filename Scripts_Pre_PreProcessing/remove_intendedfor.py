

import os 
import glob
import json
'''Remove IntendedFor statement from json files if a mistake was made'''


#grab the epi jsons and store them in a list

files = os.path.join('/data', 'D2', 'Nifti_OC', 'sub-*', 'ses-1', 'fmap', '*.json')
fmaps = glob.glob(files)



#next search through each fmap json and remove 'IntendedFor'

for file in fmaps:

	with open(file, 'r') as data_file:
		data = json.load(data_file)
	try:
		print(f"removing IntendedFor statement from {file}")
		del data['IntendedFor']
	except KeyError:
		print(f"ERROR: could not find IntendedFor statement for {file}")
	with open(file, 'w') as data_file:
		data = json.dump(data, data_file)
