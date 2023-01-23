import os 
import glob
import json

def remove_intended_for(path):
    """
    Remove 'IntendedFor' statement from json files if a mistake was made
    :param path: The path where the BIDS formatted subjects are located
    """
    # Get the list of fieldmaps json files
    files = os.path.join(path, 'sub-*', 'ses-1', 'fmap', '*.json')
    fmaps = glob.glob(files)

    # Search through each fmap json and remove 'IntendedFor'
    for file in fmaps:
        with open(file, 'r') as data_file:
            data = json.load(data_file)
        try:
            print(f"Removing 'IntendedFor' statement from {file}")
            del data['IntendedFor']
        except KeyError:
            print(f"ERROR: Could not find 'IntendedFor' statement for {file}")
        with open(file, 'w') as data_file:
            json.dump(data, data_file)
	
def usage_example():
    """
    Example of how to use the 'remove_intended_for' function
    """
    print("path = '/data/D2/Nifti_OC' ")
    print("remove_intended_for(path)")

if __name__ == "__main__":	
	import sys

	path = sys.argv[1]
	remove_inteded_for(path)
