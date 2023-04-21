import os
import configparser
import scipy.io as sp
from nipype import Workflow, Node, Function
from nipype.interfaces import DataSink
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.spm as spm

#Local Imports
from workflows import create_smooth_despike_workflow, level1analysis
from util_encore_temp import load_design_matrix, create_design_matrix, datasource 
import wrappers as wrap

def task(directory, configuration_file, session=1):
    """
    This function defines a workflow for a neuroimaging task using the Nipype package in Python. 
    The function takes two arguments: the directory where the data is stored and a configuration 
    file specifying the details of the task.
    The function first imports necessary modules and defines variables by parsing the configuration file. 
    It then creates a workflow object with a base directory, and sets up the first level analysis node 
    with the specified contrasts. If there are multiple runs, the function sets up merge points for the 
    functional data, design matrices, movement data, and regressors.
    For each run, the function sets up nodes for data processing and analysis, such as data smoothing and 
    creating design matrices. If there is only one run, the function connects these components into a pipeline. 
    If there are multiple runs, the function merges the data from each run before connecting to the first level 
    analysis node.
    The function also sets up a datasink node to store output files, and defines regions of interest for 
    psychophysiological interaction (PPI) analysis. For each ROI, the function sets up a PPI node and an 
    estimation of contrasts node, before outputting the contrast images to the datasink.
    Finally, the function defines the output files to be saved in the datasink and returns the workflow object.
    
    Parameters:
        directory (str): A string with the path to the directory containing the data.
        configuration_file (str): A string with the path to the configuration file.
        session (int): Session number to be processed. Default=1

    Returns:
        task (pe.Workflow): The Nipype workflow for the entire 1st level analysis.
        base_dir (str): Path to analysis folder to be removed later

    Raises:
        FileNotFoundError: If the `directory` or `configuration_file` does not exist.
        ValueError: If the configuration file has missing or invalid parameters.

    """
    from util_encore_temp import load_design_matrix, create_design_matrix, datasource 
    print("TESTING")
    session = str(session)
    fsl.FSLCommand.set_default_output_type('NIFTI')
    subject = os.path.basename(directory)
    config = configparser.ConfigParser()
    config.read(configuration_file)
    script_dir = os.path.abspath(os.path.dirname(__file__))

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
    try: 
        BHV_struct = eval(config['Task']['BHV_struct'])
    except:
        BHV_struct = 'BHV*/Scan'

    # Define base directory
    base_dir = os.path.abspath(os.path.join(directory, f"ses-{session}",'analysis', output_name))
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    out_dir = os.path.abspath(os.path.join(directory, f"ses-{session}", 'output', output_name, sequence))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    task = Workflow(name=sequence, base_dir=base_dir)
    l1 = level1analysis(configuration_file)
    l1.inputs.input.contrasts = contrasts
    #If there are multiple runs, set up merge points
    if runs > 1:
        # setup merge points
        merge_func = pe.Node(name="merge_func",interface=util.Merge(runs))
        merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(runs))
        merge_move = pe.Node(name="merge_movement",interface=util.Merge(runs))
        merge_regressors = pe.Node(name="merge_regressors",interface=util.Merge(runs))
    #Loop over runs
    for run in range(1, runs+1):
        # Setup nodes to be used in Pipeline
        sequence = eval(config.get('Task', 'sequence'))
        run_str = str(run)
        if runs > 1:
            sequence = sequence+run_str
        # Define datasource and smoothed task nodes
        ds = datasource(directory, sequence, session, BHV_struct)
        smoothed_task = create_smooth_despike_workflow(directory, sequence, base_dir, ds)
        # Define mCompCor, create_DM, and connect nodes
        cc = Node(interface=wrap.mCompCor(), name='mCompCor'+run_str)
        cc.inputs.white_mask = ROIS['ROI_white']

        dm = Node(interface=Function(input_names=['matlab_function', 'eprime_file', 'sequence'],
                                    output_names=['design_matrix'],
                                    function=create_design_matrix), name='create_DM'+run_str)
        if 'paat' in sequence:
            dm.inputs.matlab_function = design_script
        elif run >= 2:
            dm.inputs.matlab_function = design_script.replace('1',str(run))
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
        task.connect(pppi,'beta_images',datasink,"gPPI.pppi_"+roi[0]+".@_beta_images") 
        task.connect(pppi,'residual_image',datasink,"gPPI.pppi_"+roi[0]+".@_residual_image")
        task.connect(contrast,"spm_mat_file",datasink,"gPPI.pppi_"+roi[0]+".@_spm_file")
        task.connect(contrast,"spmT_images",datasink,"gPPI.pppi_"+roi[0]+".@_spmT_images")


    # Define the output files to be saved in datasink
    output_names = ['spm_mat_file_con', 'con_images', 'spmT_images', 'beta_images', 'residual_image']
    for output_name in output_names:
        task.connect([(l1, datasink, [(('output.' + output_name), f'level1.@{output_name}')])])
    return task, base_dir

if __name__ == '__main__':
    import sys
    import os
    import shutil
    from util_encore_temp import create_motion_file
    from task import task
    import nipype
    
        
    # Define argument names and order
    script_name, config_file, sequence_name, subject_path, session = sys.argv
    
    # Convert session to integer, set default value of 1
    session = int(session) if session else 1
    
    # Call create_motion_file and task functions with the arguments
    # Define a dictionary of sequence names and corresponding motion file creation calls
    sequence_calls = {
        'efnback': ['efnback1', 'efnback2'],
        'reward': ['reward1', 'reward2'],
        'paat': ['paat1', 'paat2', 'paat3']
    }

    # Call create_motion_file with the appropriate arguments based on the sequence name
    if sequence_name in sequence_calls:
        for seq in sequence_calls[sequence_name]:
            create_motion_file(os.path.abspath(subject_path), seq, session)
    else:
        create_motion_file(os.path.abspath(subject_path), sequence_name, session)

    script_dir = os.path.dirname(os.path.abspath(script_name))
    df, base_dir = task(os.path.abspath(subject_path), os.path.abspath(config_file), session)
    print(df)
    print(os.path.abspath(subject_path))
    print(base_dir)
    print(os.path.abspath(config_file))
    #df.run(plugin='MultiProc', plugin_args={'n_procs': 8})
    df.run()
    shutil.rmtree(base_dir)
