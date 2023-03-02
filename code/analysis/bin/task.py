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

def create_smooth_despike_workflow(directory, sequence, base_dir, datasource):
    """
    This function creates a Nipype workflow that performs despiking, masking, smoothing, and outputting of fMRI data.

    Parameters:
        directory (str): The base directory where the input data is stored.
        sequence (str): The sequence/task name for which input data needs to be retrieved.
        base_dir (str): The base directory of the workflow.
        datasource (str): The Nipype datasource node.

    Returns:
        workflow (pe.Workflow): The Nipype workflow for despiking, masking, and smoothing.
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
    """
    This function takes in a configuration file path and an optional name for level1 fMRI analysis using SPM. 
    The configuration file contains parameters for model specification, level 1 design, and contrast estimation.
    The workflow is set up using several Nipype interfaces, including `SpecifySPMModel`, `Level1Design`,
    and `EstimateContrast`. The workflow takes inputs for movement, functional runs, design matrix, contrasts, and a mask.  
    It also outputs several files including: `spm_mat_file`, `spm_mat_file_con`, `con_images`, `spmT_images`, `residual_image` 
    and `beta_images`.

    Parameters:
        configuration_file (str):   The path to the configuration file containing parameters for the model specification, 
                                    level 1 design, and contrasts estimation
        name (str): The name of the workflow (default is `level`)
    Returns:
        l1analysis (pe.Workflow): The Nipype workflow for level 1 fMRI analysis using SPM
    """
    import nipype.interfaces.spm as spm
    import nipype.interfaces.matlab as matlab
    import nipype.interfaces.utility as util
    import nipype.pipeline.engine as pe
    import nipype.algorithms.modelgen as model
    from utils import load_design_matrix

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

def task(directory, configuration_file):
    import wrappers as wrap
    from utils import load_design_matrix, create_design_matrix, datasource 
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
    base_dir = os.path.abspath(os.path.join(directory, 'analysis', output_name))
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
    from utils import create_motion_file
    #create_motion_file('/data/github/ITAlics_Developmental/data/sub-50225/', 'dynface')	
    #df = task('/data/github/ITAlics_Developmental/data/sub-50225/', '/data/github/ITAlics_Developmental/code/analysis/bin/dynfaces_file.cfg')
    #create_motion_file('/data/github/ITAlics_Developmental/data/sub-50225/', 'efnback1')
    #create_motion_file('/data/github/ITAlics_Developmental/data/sub-50225/', 'efnback2')	
    #df = task('/data/github/ITAlics_Developmental/data/sub-50225/', '/data/github/ITAlics_Developmental/code/analysis/bin/efnback_file.cfg')
    config_file = sys.argv[1]
    sequence = sys.argv[2]
    subject = sys.argv[3]
    
    create_motion_file(os.path.abspath(subject), sequence)
    df = task(os.path.abspath(subject), os.path.abspath(config_file))
    df.run(plugin='MultiProc', plugin_args={'n_procs' : 8})
