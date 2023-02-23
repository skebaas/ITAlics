import os
import glob
import json

def addIntendedforStatement(subject):
    """
    Add 'IntendedFor' condition to fieldmaps for a given subject
    :param subject: The path of the subject
    """
    subject_name = os.path.basename(subject)
    print(f"Adding 'IntendedFor' condition to fieldmaps for {subject_name}")

    # Define the path for fieldmaps and functional files
    fmaps_path = os.path.join(subject, 'ses-1', 'fmap', '*.json')
    funcs_path = os.path.join(subject, 'ses-1', 'func', '*.nii*')

    # Get the list of fieldmaps and functional files
    fmaps = glob.glob(fmaps_path)
    funcs = glob.glob(funcs_path)

    # Remove the subject path from the functional files' path
    funcs = list(map(lambda x: x.replace(f'{subject}/', ''), funcs))

    # Add the 'IntendedFor' key to the fieldmaps' json file
    for fmap in fmaps:
        with open(fmap, 'r') as data_file:
            fmap_json = json.load(data_file)
        fmap_json['IntendedFor'] = funcs

        with open(fmap, 'w') as data_file:
            json.dump(fmap_json, data_file, indent=4)

if __name__ == "__main__":	
	import sys

	subjectPath = sys.argv[1]
	subject = glob.glob(subjectPath)
	addIntendedforStatement(subject[0])
	
