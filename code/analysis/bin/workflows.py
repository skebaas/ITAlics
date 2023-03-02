import os
import configparser
import scipy.io as sp
from nipype import Workflow, Node, Function
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.spm as spm
import nipype.interfaces.utility as util
from nipype.interfaces.utility import IdentityInterface
import nipype.pipeline.engine as pe

def create_smooth_despike_workflow(directory, sequence, base_dir, datasource):
    """
    This function creates a Nipype workflow that performs despiking, masking, smoothing, and outputting of fMRI data.

    Args:
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
    This function takes in a configuration file path and an optional name forlevel1 fMRI analysis using SPM. 
    The configuration file contains parameters for model specification, leve 1 design, and contrast estimation.
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