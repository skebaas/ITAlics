import scipy.io as sp
import re
import os
import glob as gl
import numpy
import nipype.interfaces.matlab as mlab
from nipype.interfaces.base import Bunch
import chardet
import subprocess
import glob
import shutil
import pandas as pd
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.afni as afni

def load_design_matrix(mat_files, trim=0):
    """
    Load design matrices from one or more Matlab .mat files.

    Args:
        mat_files (str): a path to a .mat file or a list of paths to .mat files.
        trim (int): the number of columns to trim from the end of the design matrix (default 0).

    Returns:
        A list of Bunch objects, where each Bunch corresponds to a loaded design matrix.
        Each Bunch has the following attributes:
        conditions(list): a list of condition names
        onsets(list): a list of lists of onset times
        durations(list): a list of lists of durations
        pmod(list): a list of Bunch objects, each corresponding to a parametric modulator
    """
    import scipy.io as sp
    import re
    import os
    import glob as gl
    import numpy
    import nipype.interfaces.matlab as mlab
    from nipype.interfaces.base import Bunch
    # convert numpy data array
    def convert_numpy(ar,toString=False):
        import numpy
        lst = []
        for a in ar.tolist():
            if toString:
                lst.append(str(a))
            elif type(a) == numpy.ndarray:
                lst.append(a.tolist())
            else:
                lst.append([a])
        return lst

    if not isinstance(mat_files, list):
        mat_files = [mat_files]

    bunches = []
	
    for mat_file in mat_files:
        # load design matrix 
        dm = sp.loadmat(mat_file,squeeze_me=True)
        names = convert_numpy(dm.get('names'), True)
        onsets = convert_numpy(dm.get('onsets'))
        durations = convert_numpy(dm.get('durations'))

        # trim last columns from design matrix
        if trim > 0:
            names = names[0:-trim]
            durations = durations[0:-trim]
            onsets = onsets[0:-trim]
        bunch = Bunch(conditions=names,onsets=onsets,durations=durations)
        # load pmod if it exists
        if 'pmod' in dm:
            pmod = []
            for i in range(len(dm['pmod'])):
                pmod_i = dm['pmod'][i]
                if isinstance(pmod_i['name'], numpy.ndarray):
                    # handle case where pmod has multiple entries
                    names = [str(name) for name in pmod_i['name']]
                    params = pmod_i['param'].tolist()
                    polys = pmod_i['poly']
                else:
                    names = [str(pmod_i['name'])]
                    params = [pmod_i['param'].tolist()]
                    polys = [pmod_i['poly']]
                pmod.append(Bunch(name=names, param=params, poly=polys))
            bunch.pmod = pmod
        # create bunch to return
        #bunch = Bunch(conditions=names, onsets=onsets, durations=durations, pmod=pmod)
        bunches.append(bunch)
    return bunches

def create_design_matrix(matlab_function, eprime_file, sequence):
    """
    Creates a design matrix using a MATLAB function, an eprime file and a sequence.

    Args:
        matlab_function (str): Name of the MATLAB function to be executed
        eprime_file (str): File path of the eprime file
        sequence (str): Sequence identifier

    Returns:
        mat (str): File path of the nDM file
    """
    import os
    import glob as gl
    import chardet
    import subprocess
    import nipype.interfaces.matlab as mlab
    # Check if nDM file is already available
    mat = os.path.join(os.path.dirname(eprime_file), f"nDM*{sequence}*.mat")
    mat = gl.glob(mat)
    if len(mat) > 0:
        return mat[0]
    
    # Check encoding of eprime file
    # If UTF-16 encoding, then change it to ASCII (MATLAB-2015 cannot read UTF-16 encoding)
    with open(eprime_file, 'rb') as f:
        contents = f.read()
        detect_f = chardet.detect(contents) # Determine encoding
    if detect_f['encoding'] != 'ascii':
        decoded = contents.decode('utf-16')
        encoded = decoded.encode('ascii')
        with open(eprime_file, 'wb') as f:
            f.write(encoded)
        del contents
    # execute matlab script to generate nDM file
    m = mlab.MatlabCommand()
    m.inputs.mfile = False
    m.inputs.script = matlab_function+"(\'"+eprime_file+"\');"
    m.run()
    # Execute MATLAB script to generate nDM file
    #matlab_script = matlab_function + "(\'" + eprime_file + "\');"
    #subprocess.run(["matlab", "-nodisplay", "-nosplash", "-nodesktop", "-r", matlab_script], check=True)

    # Get nDM file (again)
    mat = os.path.join(os.path.dirname(eprime_file), f"nDM*{sequence}*.mat")
    mat = gl.glob(mat)
    if len(mat) > 0:
        mat = mat[0]

    return mat

def datasource(directory, sequence):
    """
    This function creates a Nipype data source node that retrieves input data for a given subject and sequence.
    The input data includes functional and structural MRI scans, as well as task-specific behavior files and motion parameters.

    Args:
        directory (str): The base directory where the input data is stored.
        sequence (str): The sequence/task name for which input data needs to be retrieved.

    Returns:
        A Nipype data source node that can be used as an input node in a Nipype pipeline.
    """

    # Get subject ID from the base directory
    subject = os.path.basename(directory)

    # Define the output fields that need to be retrieved
    outfields = ['func', 'struct']
    
    # Define the field template for functional and structural scans
    field_template = dict(func=f"ses-1/func/*{sequence}*preproc_bold.nii*",
                          struct="ses-1/anat/*desc-preproc_T1w.nii*")
    
    # Define template arguments for functional and structural scans
    template_args = dict(func=[[]], struct=[[]])

    # If the sequence is a resting-state scan, add mask and confounds to the output fields
    if sequence.startswith('resting'):
        field_template.update(mask="ses-1/anat/*space-MNI152NLin6Asym_res-2_desc-brain_mask.nii*",
                              confounds="ses-1/func/*task-rest*desc-confounds*.tsv")
        outfields.extend(['mask', 'confounds'])
    
    # If the sequence is a task-specific scan, add behavior files, mask file, and motion parameters to the output fields
    elif sequence.startswith('reward') or sequence.startswith('efnback') or sequence.startswith('dynface'):
        behavior_files = {'reward1': 'Diamond_Reward*-1', 'reward2': 'Diamond_Reward*-2', 'efnback1': 'EFNBACK*-1', 'efnback2': 'EFNBACK*-2', 'dynface': 'subject*'}
        behavior_file = behavior_files.get(sequence, None)
        if behavior_file is not None:
            field_template['behav'] = f"BHV*/Scan/{behavior_file}.txt"
            template_args['behav'] = [[]]
            outfields.append('behav')
        else:
            print(f"No behavior file found for {sequence}... Exiting now...")
        field_template['mask_file'] = f"ses-1/func/*{sequence}*brain_mask.nii*"
        template_args['mask_file'] = [[]]
        outfields.append('mask_file')
        field_template['movement'] = f"motion/{sequence}_motion.1D"
        template_args['movement'] = [[]]
        outfields.append('movement')

    # Create a Nipype data source node with the defined output fields and field templates
    datasource = pe.Node(interface=nio.DataGrabber(
                         infields=['subject_id', 'sequence'], outfields=outfields),
                         name=f"datasource_{sequence}")
    datasource.inputs.base_directory = os.path.abspath(directory)
    datasource.inputs.template = '*'
    datasource.inputs.field_template = field_template
    datasource.inputs.template_args = template_args
    datasource.inputs.subject_id = subject
    datasource.inputs.sequence = sequence
    datasource.inputs.sort_filelist = True

    return datasource

def create_motion_file(directory: str, sequence: str) -> None:
    """
    Function to create a motion file by extracting specific columns from a fmriprep TSV file.

    Args:
        directory (str): The path to the directory containing the fmriprep TSV file.
        sequence (str): The sequence name used to identify the TSV file.

    Returns:
        None: The function does not return anything, but saves a motion file in the specified directory.

    """

    # Define a helper function to check if a folder already exists.
    def create_folder(foldername: str) -> None:
        """
        Function to check if a folder already exists. If it does not exist, it will create one.

        Args:
            foldername (str): The name of the folder to be created.

        Returns:
            None: The function does not return anything, but creates a new folder if it does not already exist.

        """
        if not os.path.exists(foldername):
            os.mkdir(foldername)
    
    # Create a new folder to store the motion file.
    motion_folder = os.path.join(directory, 'motion')
    create_folder(motion_folder)

    # Define the filename and path for the new motion file.
    motion_file = os.path.join(motion_folder, sequence + '_motion.1D')

    # Find the fmriprep TSV file with the specified sequence name and read it into a pandas DataFrame.
    confounds_tsv = glob.glob(os.path.join(directory, 'ses-1', 'func', f'*{sequence}*desc-confounds_timeseries.tsv'))
    confounds_df = pd.read_csv(confounds_tsv[0], sep='\t')

    # Extract the specific columns for motion parameters and save them to a new file.
    motion_pars = confounds_df[['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']]
    motion_pars.to_csv(motion_file, header=None, index=None, sep='\t')

    # Print a message to indicate that the motion file was successfully created.
    print(f"Created motion file for {sequence} at {motion_file}")
