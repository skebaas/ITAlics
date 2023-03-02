import os
import configparser
import re
import glob as gl
import numpy
import scipy.io as sp
import nipype.interfaces.matlab as mlab
from nipype.interfaces.base import Bunch
from nipype import Workflow, Node, Function
from nipype.interfaces import fsl, DataSink
from nipype.interfaces.utility import IdentityInterface
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.afni as afni
import chardet
import subprocess
import scipy.io as sp
import nipype.interfaces.spm as spm

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

def create_smooth_despike_workflow(directory, sequence, base_dir, datasource):
    """
    This function creates a Nipype workflow that performs despiking, masking, smoothing, and outputting of fMRI data.

    Args:
        directory (str): The base directory where the input data is stored.
        sequence (str): The sequence/task name for which input data needs to be retrieved.
        base_dir (str): The base directory of the workflow.
        datasource (str): The Nipype datasource node.

    Returns:
        A Nipype data workflow.
    """
    # Create a new workflow
    workflow = pe.Workflow(name='smooth_despike_node_' + sequence, base_dir=base_dir)

    # Despiking node
    despike = pe.Node(afni.Despike(outputtype='NIFTI'), name='despike')
    workflow.connect(datasource, 'func', despike, 'in_file')

    # Mask application node
    apply_mask = pe.Node(fsl.ApplyMask(output_type='NIFTI'), name='apply_mask')
    workflow.connect(despike, 'out_file', apply_mask, 'in_file')
    workflow.connect(datasource, 'mask_file', apply_mask, 'mask_file')

    # Smoothing node
    smoothing = pe.Node(fsl.Smooth(fwhm=6.0, output_type='NIFTI'), name='smooth')
    workflow.connect(apply_mask, 'out_file', smoothing, 'in_file')

    # Output node
    outputnode = pe.Node(util.IdentityInterface(fields=['func', 'mask', 'struct', 'movement', 'behav', 'ufunc']), name='output')
    workflow.connect(smoothing, 'smoothed_file', outputnode, 'func')
    workflow.connect(datasource, 'struct', outputnode, 'struct')
    workflow.connect(datasource, 'mask_file', outputnode, 'mask')
    workflow.connect(datasource, 'behav', outputnode, 'behav')
    workflow.connect(datasource, 'func', outputnode, 'ufunc')
    workflow.connect(datasource, 'movement', outputnode, 'movement')

    # Set the base directory of the data source
    #datasource.inputs.base_directory = os.path.abspath(directory)

    return workflow

def level1analysis(configuration_file, name='level1'):
    import nipype.interfaces.spm as spm
    import nipype.interfaces.matlab as matlab
    import nipype.interfaces.utility as util
    import nipype.pipeline.engine as pe
    import nipype.algorithms.modelgen as model

    l1analysis = pe.Workflow(name=name)
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['movement', 'func', 'design_matrix', 'contrasts', 'mask']), name='input')

    # load config
    config = configparser.ConfigParser()
    config.read(configuration_file)

    # specify design matrix model
    modelspec = pe.Node(interface=model.SpecifySPMModel(), name="modelspec")
    modelspec.inputs.concatenate_runs = eval(config['modelspec']['concatenate_runs'])
    modelspec.inputs.time_repetition = eval(config['modelspec']['time_repetition'])
    modelspec.inputs.high_pass_filter_cutoff = eval(config['modelspec']['high_pass_filter_cutoff'])
    modelspec.inputs.input_units = eval(config['modelspec']['input_units'])
    l1analysis.connect(inputnode, 'movement', modelspec, 'realignment_parameters')
    l1analysis.connect(inputnode, 'func', modelspec, 'functional_runs')
    l1analysis.connect(inputnode, ('design_matrix', load_design_matrix, 0), modelspec, 'subject_info')

    # create design matrix
    level1design = pe.Node(interface=spm.Level1Design(), name="level1design")
    level1design.inputs.bases = eval(config['level1design']['bases'])
    level1design.inputs.timing_units = eval(config['level1design']['timing_units'])
    level1design.inputs.interscan_interval = eval(config['modelspec']['time_repetition'])
    level1design.inputs.microtime_onset = eval(config['level1design']['microtime_onset'])
    level1design.inputs.microtime_resolution = eval(config['level1design']['microtime_resolution'])
    level1design.inputs.model_serial_correlations = eval(config['level1design']['model_serial_correlations'])

    # Incorporate fmriprep mask in spm model
    l1analysis.connect(inputnode, 'mask', level1design, 'mask_image')
    l1analysis.connect(modelspec, 'session_info', level1design, 'session_info')

    # level 1 estimate
    level1estimate = pe.Node(interface=spm.EstimateModel(), name="level1estimate")
    level1estimate.inputs.estimation_method = eval(config['level1estimate']['estimation_method'])
    level1estimate.inputs.write_residuals = False
    l1analysis.connect(level1design, 'spm_mat_file', level1estimate, 'spm_mat_file')

    # no need for contrast for pppi model
    contrastestimate = pe.Node(interface=spm.EstimateContrast(), name="contrastestimate")
    contrastestimate.inputs.use_derivs = eval(config['contrastestimate']['use_derivs'])
    l1analysis.connect(inputnode, 'contrasts', contrastestimate, 'contrasts')
    l1analysis.connect(level1estimate,'spm_mat_file',contrastestimate,'spm_mat_file')
    l1analysis.connect(level1estimate,'beta_images',contrastestimate,'beta_images'),
    l1analysis.connect(level1estimate,'residual_image',contrastestimate,'residual_image')
    # output
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['spm_mat_file', 'spm_mat_file_con','con_images', 'spmT_images', 'residual_image','beta_images']), name='output')
    l1analysis.connect(level1estimate,'spm_mat_file',outputnode,'spm_mat_file')
    l1analysis.connect(level1estimate,'beta_images',outputnode,'beta_images')
    l1analysis.connect(contrastestimate,'spm_mat_file',outputnode,'spm_mat_file_con')
    l1analysis.connect(contrastestimate,'con_images',outputnode,'con_images')
    l1analysis.connect(contrastestimate,'spmT_images',outputnode,'spmT_images')

    return l1analysis

def create_motion_file(directory: str, sequence: str) -> None:
    """
    Function to create a motion file by extracting specific columns from a fmriprep TSV file.

    Args:
        directory (str): The path to the directory containing the fmriprep TSV file.
        sequence (str): The sequence name used to identify the TSV file.

    Returns:
        None: The function does not return anything, but saves a motion file in the specified directory.

    """
    import os
    import glob
    import shutil
    import pandas as pd

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

def task(directory, configuration_file):
    import wrappers as wrap
    fsl.FSLCommand.set_default_output_type('NIFTI')
    subject = os.path.basename(directory)
    config = configparser.ConfigParser()
    config.read(configuration_file)

    #Extract following variables from configuration file
    sequence = eval(config.get('Task', 'sequence'))
    runs = eval(config.get('Task', 'runs'))
    contrasts = eval(config.get('Task','contrasts',raw=False))
    ppi_contrasts = eval(config.get('Task','ppi_contrasts',raw=False))
    design_script = eval(config['Task']['design_script'])
    data_dir = eval(config['Default']['data_dir'])
    ROI_dir = eval(config['Default']['ROI_dir'])
    ROI_suffix = eval(config['Default']['ROI_suffix'])
    ROIS = eval(config['Default']['ROIS'])
    PPI_ROIS = eval(config['Task']['PPI_ROIS'])
    output_name = eval(config['Task']['output_name'])

    # Define base directory
    base_dir = os.path.abspath(os.path.join(directory, 'analysis', output_name, sequence))
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    out_dir = os.path.abspath(os.path.join(directory, 'output', output_name, sequence))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    task = Workflow(name=sequence, base_dir=base_dir)
    l1 = level1analysis(configuration_file)
    l1.inputs.input.contrasts = contrasts
    if runs > 1:
        # setup merge points
        merge_func = pe.Node(name="merge_func",interface=util.Merge(runs))
        merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(runs))
        merge_move = pe.Node(name="merge_movement",interface=util.Merge(runs))
        merge_regressors = pe.Node(name="merge_regressors",interface=util.Merge(runs))

    for run in range(1, runs+1):
        # Setup nodes to be used in Pipeline
        sequence = eval(config.get('Task', 'sequence'))
        run_str = str(run)
        if runs > 1:
            sequence = sequence+run_str
        ds = datasource(directory, sequence)
        smoothed_task = create_smooth_despike_workflow(directory, sequence, base_dir, ds)

        cc = Node(interface=wrap.mCompCor(), name='mCompCor'+run_str)
        cc.inputs.white_mask = ROIS['ROI_white']

        dm = Node(interface=Function(input_names=['matlab_function', 'eprime_file', 'sequence'],
                                    output_names=['design_matrix'],
                                    function=create_design_matrix), name='create_DM'+run_str)
        if run == 2:
            dm.inputs.matlab_function = design_script.replace('1','2')
        else:
            dm.inputs.matlab_function = design_script
        dm.inputs.sequence = sequence
        if runs == 1:
            # Connect components into a pipeline
            task.connect([(smoothed_task, l1, [('output.func', 'input.func')]),
                        (ds, dm, [('behav', 'eprime_file')]),
                        (smoothed_task, cc, [('output.ufunc', 'source'), ('output.mask', 'brain_mask'), ('output.movement', 'movement')]),
                        (dm, l1, [('design_matrix', 'input.design_matrix')]),
                        (cc, l1, [('regressors', 'input.movement')]),
                        (smoothed_task, l1, [('output.mask', 'input.mask')])])
        else:
            task.connect([(ds,cc,[('func','source'),('mask_file','brain_mask'),('movement','movement')])])
            task.connect(cc, 'regressors',merge_regressors,f'in{run_str}')
            task.connect(ds,'behav',dm,"eprime_file")	
            task.connect(smoothed_task,'output.func',merge_func,f'in{run_str}')
            task.connect(smoothed_task,'output.movement',merge_move,f'in{run_str}')
            task.connect(dm,'design_matrix',merge_nDM,f'in{run_str}')

    if runs > 1:
        task.connect([(merge_func, l1, [('out', 'input.func')]),
                    (merge_nDM, l1, [('out', 'input.design_matrix')]),
                    (merge_regressors, l1, [('out', 'input.movement')]),
                    (smoothed_task, l1, [('output.mask', 'input.mask')])])
    
    # Define datasink
    datasink = Node(interface=DataSink(), name='datasink')
    datasink.inputs.base_directory = out_dir

    for roi in PPI_ROIS:
        # define PPPI for each ROI
        pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
        pppi.inputs.voi_name = roi[0]
        pppi.inputs.voi_file = roi[1]
        pppi.inputs.subject = subject
        task.connect(l1,'output.spm_mat_file',pppi,'spm_mat_file')
        task.connect(ds,'mask_file',pppi,'mask_file')

        # estimate contrast for each ROI
        contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
        contrast.inputs.contrasts = ppi_contrasts
        task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
        task.connect(pppi,'beta_images',contrast,'beta_images')
        task.connect(pppi,'residual_image',contrast,'residual_image')

        # output con images to datasink
        task.connect(contrast,'con_images',datasink,"gPPI.pppi_"+roi[0]+".@_con_images")
        task.connect(contrast,"spm_mat_file",datasink,"gPPI.pppi_"+roi[0]+".@_spm_file")
        task.connect(contrast,"spmT_images",datasink,"gPPI.pppi_"+roi[0]+".@_spmT_images")


    # Define the output files to be saved in datasink
    output_names = ['spm_mat_file_con', 'con_images', 'spmT_images', 'beta_images']
    for output_name in output_names:
        task.connect([(l1, datasink, [(('output.' + output_name), f'level1.@{output_name}')])])
    return task

if __name__ == '__main__':
    import sys
    #create_motion_file('/data/github/ITAlics_Developmental/data/sub-50225/', 'dynface')	
    #df = task('/data/github/ITAlics_Developmental/data/sub-50225/', '/data/github/ITAlics_Developmental/code/analysis/bin/dynfaces_file.cfg')
    #create_motion_file('/data/github/ITAlics_Developmental/data/sub-50225/', 'efnback1')
    #create_motion_file('/data/github/ITAlics_Developmental/data/sub-50225/', 'efnback2')	
    #df = task('/data/github/ITAlics_Developmental/data/sub-50225/', '/data/github/ITAlics_Developmental/code/analysis/bin/efnback_file.cfg')
    config_file = sys.argv[1]
    sequence = sys.argv[2]
    subject = sys.argv[3]
    
    create_motion_file(os.path.abspath(subject), sequence)
    df = task(os.path.abspath(subject), os.path.abspath(config_file), sequence)
    df.run(plugin='MultiProc', plugin_args={'n_procs' : 8})