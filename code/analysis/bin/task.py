import os
from nipype import Workflow, Node, Function
from nipype.interfaces import fsl, DataSink
from nipype.interfaces.utility import IdentityInterface
import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.afni as afni
import configparser

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
    subject = get_subject(directory)

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
        behavior_files = {'reward1': '-1', 'reward2': '-2', 'efnback1': '-1', 'efnback2': '-2', 'dynface': ''}
        behavior_file = behavior_files.get(sequence, None)
        if behavior_file is not None:
            field_template['behav'] = f"BHV*/Scan/*{sequence}*{behavior_file}.txt"
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
    datasource.inputs.base_directory = os.path.abspath(directory)

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
    modelspec.inputs.concatenate_runs = config['modelspec']['concatenate_runs']
    modelspec.inputs.time_repetition = config['modelspec']['time_repetition']
    modelspec.inputs.high_pass_filter_cutoff = config['modelspec']['high_pass_filter_cutoff']
    modelspec.inputs.input_units = config['modelspec']['input_units']
    l1analysis.connect(inputnode, 'movement', modelspec, 'realignment_parameters')
    l1analysis.connect(inputnode, 'func', modelspec, 'functional_runs')
    l1analysis.connect(inputnode, ('design_matrix', load_design_matrix, 0), modelspec, 'subject_info')

    # create design matrix
    level1design = pe.Node(interface=spm.Level1Design(), name="level1design")
    level1design.inputs.bases = config['level1design']['bases']
    level1design.inputs.timing_units = config['level1design']['timing_units']
    level1design.inputs.interscan_interval = config['modelspec']['time_repetition']
    level1design.inputs.microtime_onset = config['level1design']['microtime_onset']
    level1design.inputs.microtime_resolution = config['level1design']['microtime_resolution']
    level1design.inputs.model_serial_correlations = config['level1design']['model_serial_correlations']

    # Incorporate fmriprep mask in spm model
    l1analysis.connect(inputnode, 'mask', level1design, 'mask_image')
    l1analysis.connect(modelspec, 'session_info', level1design, 'session_info')

    # level 1 estimate
    level1estimate = pe.Node(interface=spm.EstimateModel(), name="level1estimate")
    level1estimate.inputs.estimation_method = config['level1estimate']['estimation_method']
    l1analysis.connect(level1design, 'spm_mat_file', level1estimate, 'spm_mat_file')

    # no need for contrast for pppi model
    contrastestimate = pe.Node(interface=spm.EstimateContrast(), name="contrastestimate")
    contrastestimate.inputs.use_derivs = config['contrastestimate']['use_derivs']
    l1analysis.connect(inputnode, 'contrasts', contrastestimate, 'contrasts')
    l1analysis.connect(level1estimate, ['spm_mat_file', 'beta_images', 'residual_image'], contrastestimate, ['spm_mat_file', 'beta_images', 'residual_image'])

    # output
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['spm_mat_file', 'con_images', 'spmT_images', 'residual_image']), name='output')
    l1analysis.connect(contrastestimate, ['spm_mat_file', 'con_images', 'spmT_images'], outputnode, ['spm_mat_file', 'con_images', 'spmT_images'])

    return l1analysis


def task(directory, sequence, configuration_file):
    fsl.FSLCommand.set_default_output_type('NIFTI')
    config = configparser.ConfigParser()
    config.read(configuration_file)

    #Extract following variables from configuration file
    contrasts = eval(config['Task']['contrasts'])

    # Define base directory
    base_dir = os.path.abspath(os.path.join(directory, 'analysis', output_name))
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    out_dir = os.path.abspath(os.path.join(directory, 'output', output_name))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Setup nodes to be used in Pipeline
    ds = datasource(directory, sequence)
    smoothed_task = create_smooth_despike_workflow(directory, sequence, base_dir, ds)

    l1 = Node(interface=level1analysis(configuration_file), name='level1analysis')
    l1.inputs.input.contrasts = contrasts

    cc = Node(interface=mCompCor(), name='mCompCor')
    cc.inputs.white_mask = conf.ROI_white

    dm = Node(interface=Function(input_names=['matlab_function', 'eprime_file', 'sequence'],
                                  output_names=['design_matrix'],
                                  function=create_design_matrix), name='create_DM')
    dm.inputs.matlab_function = 'dynfaces_task2dm'
    dm.inputs.sequence = 'dynface'

    # Connect components into a pipeline
    task = Workflow(name=sequence, base_dir=base_dir)
    task.connect([(smoothed_task, l1, [('output_func', 'input.func')]),
                  (ds, dm, [('behav', 'eprime_file')]),
                  (smoothed_task, cc, [('func', 'source'), ('mask_file', 'brain_mask'), ('movement', 'movement')]),
                  (dm, l1, [('design_matrix', 'input.design_matrix')]),
                  (cc, l1, [('regressors', 'input.movement')]),
                  (smoothed_task, l1, [('mask_file', 'input.mask')])])

    # Define datasink
    datasink = Node(interface=DataSink(), name='datasink')
    datasink.inputs.base_directory = out_dir

    # Connect nodes to datasink
    task.connect([(l1, datasink, [('output_directory', 'container')]),
                  (l1, datasink, [('output_res4d', 'res4d.@res4d'),
                                   ('output_copes', 'copes.@copes'),
                                   ('output_varcopes', 'varcopes.@')])])
    # Define the output files to be saved in datasink
    output_names = ['parameter_estimates', 'zstat', 'copes', 'residuals', 'dof_file', 'varcopes']
    substitutions = [('_subject_id_', ''), ('_model_name_', ''), ('_task_id_', sequence + '/')]
    for output_name in output_names:
        task.connect([(l1, datasink, [(('outputspec.' + output_name), output_name)])])
        datasink.inputs.substitutions = substitutions
        datasink.inputs.container = os.path.join(sequence, subject)
    
    return task
