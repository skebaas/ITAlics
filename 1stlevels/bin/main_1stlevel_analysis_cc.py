#!/usr/bin/env python2
# DIAMOND 2.0 Analysis Pipeline 
# CHANGED TO INCORPORATE FMRIPREP
# Author: Alexander Skeba, Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import sys
import os                                  
import re
import gold

## Predefined constants
# defaults are defined in gold.py
# you need to redefine constants either here
# or in each sequence
conf = gold.Config()
conf.CPU_CORES = 16
conf.time_repetition  = 1.5
	
# func BET	
conf.bet_mask = True
conf.bet_frac = 0.6
conf.bet_robust = True
conf.bet_vertical_gradient = -0.05

# struct BET
conf.struct_bet_mask = True
conf.struct_bet_frac = 0.5
conf.struct_bet_robust = True
conf.struct_bet_vertical_gradient = 0


	
conf.flirt_cost = 'mutualinfo'
conf.flirt_bins = 256
conf.flirt_dof = 12
conf.flirt_interp = 'trilinear'
conf.flirt_searchr_x = [-180, 180]
conf.flirt_searchr_y = [-180, 180]
conf.flirt_searchr_z = [-180, 180]

conf.coregister_cost_function = "nmi"
conf.coregister_separation = [4, 2]
conf.coregister_tolerance = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]
conf.coregister_fwhm = [7, 7]
conf.prepare_fieldmap_scanner = "SIEMENS"
conf.prepare_fieldmap_delta_TE = 2.46  
	
conf.fugue_dwell_time = 0.00064 # prizma scanner is default
conf.fugue_poly_order = 3

conf.dartel_fwhm = 6 #TODO Check value
conf.dartel_voxel_size =  (2, 2, 2)
conf.susan_brightness_threshold = 100 #200.0
conf.susan_fwhm = 6

# resting low/high band pass filter parametersf
# sigma = ((1/hz)/2.35)/TR
# hz (low) =  0.008 
# hz (high) =  0.08 
# 2.35 - gaussian constant
conf.filter_image_bptf = ' -bptf 35.461 3.546'
		
conf.modelspec_concatenate_runs   = False
conf.modelspec_high_pass_filter_cutoff = 60 # reward only
conf.modelspec_input_units = 'secs'
		
conf.level1design_bases = {'hrf':{'derivs': [0,0]}}
conf.level1design_timing_units = 'secs'
conf.level1estimate_estimation_method = {'Classical' : 1}
conf.contrastestimate_use_derivs = True
conf.level1design_microtime_onset = 1
conf.level1design_microtime_resolution = 16
#conf.level1design_model_serial_correlations = 'AR(1)'
conf.level1design_model_serial_correlations = 'FAST'





# default value to use fieldmap in the pipeline
useFieldmap=False
noPrint = True
def run_script(script, stdin=None):
    """Returns (stdout, stderr), raises error on non-zero return code"""
    import subprocess
    # Note: by using a list here (['bash', ...]) you avoid quoting issues, as the
    # arguments are passed in exactly this order (spaces, quotes, and newlines won't
    # cause problems):
    proc = subprocess.Popen(['bash', '-c', script],
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode:
        raise ScriptException(proc.returncode, stdout, stderr, script)
    return stdout, stderr

class ScriptException(Exception):
    def __init__(self, returncode, stdout, stderr, script):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        Exception().__init__('Error in script')

def gunzip_files(directory):
    import os
    import glob
    directory = os.path.abspath(directory)
    niftis = os.path.join(directory, 'ses-1', '*', '*.gz')
    bash_script = "for nifti in " + niftis + " ; do gunzip $nifti; done"
    print("unzipping nifti_gz files in " + directory)
    if glob.glob(niftis):
        run_script(bash_script)

def gzip_files(directory):
    print("zipping nifti files")
    bash_script = "find " +directory+ " -name '*nii' -exec gzip {} \;"
    run_script(bash_script)
    #bash_script = "find " +directory+ " -name '*img' -exec gzip {} \;"
    #run_script(bash_script)

def createFolder(foldername):
    """
    Function to check if a folder already exists.  If it does, it will not create one
    """
    import os
    import glob
    import shutil
    if not os.path.exists(foldername):
        os.mkdir(foldername)

def createMotionFile(directory, sequence):
    """
    Function to create a motion file by pulling specific columns from the fmriprep tsv

    """
    import os
    import glob
    import shutil
    import pandas as pd

    motionFolder = os.path.join(directory, 'motion')
    createFolder(motionFolder)  #WARNING! Calling on other method in script!
    motionFile = os.path.join(motionFolder, sequence+'_motion.1D')
    confoundsTSV = glob.glob(os.path.join(directory, 'ses-1', 'func', '*'+sequence+'*desc-confounds_timeseries.tsv'))

    confoundsDF = pd.read_csv(confoundsTSV[0], sep='\t')
    motionPars = confoundsDF[['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']]
    motionPars.to_csv(motionFile, header=None, index=None, sep='\t')
    print("Created motion file for "+sequence+" at " + motionFile)

"""
DIAMOND 1.0 Input Data Source

datasource is the module that uses wildcards to pull
struct		- anatomical NIFTI file
func 		- functional NIFIT file
behav		- behaviour file that can be either an EPRIME text file or something else
fieldmap_mag	- fieldmap magnitute NIFTI file
fieldmap_phase	- fieldmap phase NIFTI file

"""
def datasource(directory, sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o
	import glob as gl
	import os

	
	# define some variables beforehand
	subject=get_subject(directory)

	# Case for resting state
	if sequence.startswith('resting'):
		field_template = dict(func="ses-1/func/*task-rest*preproc_bold.nii*",mask="ses-1/anat/*space-MNI152NLin6Asym_res-2_desc-brain_mask.nii*", confounds="ses-1/func/*task-rest*desc-confounds*.tsv")
		template_args  = dict(func=[[]],mask=[[]],confounds=[[]])
		
		outfields=['func', 'mask', 'confounds']
		datasource = pe.Node(interface=nio.DataGrabber(
						 infields=['subject_id','sequence'], 
						 outfields=outfields),
	                     name = "datasource_"+sequence)
		
		datasource.inputs.base_directory = os.path.abspath(directory)
		datasource.inputs.template = '*'
		datasource.inputs.field_template = field_template
		datasource.inputs.template_args  = template_args
		datasource.inputs.subject_id = subject
		datasource.inputs.sequence = sequence
		datasource.inputs.sort_filelist = True
		
		return datasource
	
	# figure out if this is new or old naming convention file
	orig = False	
	m = gl.glob(os.path.join(directory,"anat/T1MPRAGE*[0-9].nii"))
	if len(m) > 0:
		orig = True


	outfields=['func', 'struct']


	# define templates for datasource for functional and structural images
	field_template = dict(func="ses-1/func/*"+sequence+"*preproc_bold.nii*",struct="ses-1/anat/*desc-preproc_T1w.nii*")
	template_args  = dict(func=[[]],struct=[[]])                


	# add behavior file to task oriented design
	if sequence.startswith('reward') or sequence.startswith('efnback') or sequence.startswith('dynface'):
		if sequence.startswith('reward1'):
			field_template['behav'] = "BHV*/Scan/*Reward*-1.txt"
		elif sequence.startswith('reward2'):
			field_template['behav'] = "BHV*/Scan/*Reward*-2.txt"
		elif sequence.startswith('efnback1'):
			field_template['behav'] = "BHV*/Scan/EFNBACK*-1.txt"
		elif sequence.startswith('efnback2'):
			field_template['behav'] = "BHV*/Scan/EFNBACK*-2.txt"
		elif sequence.startswith('dynface'):
			field_template['behav'] = "BHV*/Scan/subject*.txt"
		else:
			print('No bhv found for ' +sequence+'... Exiting now...')
		template_args['behav']  = [[]]
		outfields.append('behav')
		field_template['mask_file'] = "ses-1/func/*"+sequence+"*brain_mask.nii*"
		outfields.append('mask_file')
		template_args['mask_file'] = [[]]
		field_template['movement'] ="motion/"+sequence+"_motion.1D"
		template_args['movement'] = [[]]
		outfields.append('movement')

			
	# specify input dataset just pass through parameters
	datasource = pe.Node(interface=nio.DataGrabber(
						 infields=['subject_id','sequence'], 
						 outfields=outfields),
	                     name = "datasource_"+sequence)
	datasource.inputs.base_directory = os.path.abspath(directory)
	datasource.inputs.template = '*'
	datasource.inputs.field_template = field_template
	datasource.inputs.template_args  = template_args
	datasource.inputs.subject_id = subject
	datasource.inputs.sequence = sequence
	datasource.inputs.sort_filelist = True

	return datasource


"""
get subject name from a directory
"""
def get_subject(directory):
	m = re.search('diamond_[0-9]+.([0-9]+)',directory)
	if m:
		return m.group(1)
	return "subject"
"""
get a subset of a list, internal function used within pipeline
"""
def subset(x,i):
	return x[i]	

def smooth_despike_node(directory,sequence,base_dir,datasource):
	import nipype.pipeline.engine as pe          # pypeline engine
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import gold
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.fsl.maths as math  
	import nipype.interfaces.afni as afni	

	smooth_despike_node = pe.Workflow(name='smooth_despike_node_' + sequence)
	smooth_despike_node.base_dir = base_dir


	#Despiking

	despike = pe.Node(interface=afni.Despike(), name='despike')
	despike.inputs.outputtype = 'NIFTI'

	smooth_despike_node.connect(datasource,'func',despike,'in_file')

	
	# apply mask to 4d image
	apply_mask = pe.Node(interface=fsl.ApplyMask(), name="apply_mask")
	apply_mask.inputs.output_type = 'NIFTI'
	smooth_despike_node.connect(despike,'out_file',apply_mask,'in_file')	
	smooth_despike_node.connect(datasource,'mask_file',apply_mask,'mask_file')
	
    #smoothing    
	smoothing = pe.Node(interface=fsl.Smooth(), name="smooth")
	smoothing.inputs.fwhm = 6.0
	smoothing.inputs.output_type = 'NIFTI'
	smooth_despike_node.connect(apply_mask,'out_file',smoothing,'in_file')
	outputnode = pe.Node(interface=util.IdentityInterface(fields=['func','mask','struct', 'movement','behav','ufunc']),name='output')

	#Use this node to access the files you need from ds1
	smooth_despike_node.connect(smoothing,'smoothed_file',outputnode,'func')
	smooth_despike_node.connect(datasource,'struct',outputnode,'struct')
	smooth_despike_node.connect(datasource,'mask_file',outputnode,'mask')
	smooth_despike_node.connect(datasource,'behav',outputnode,'behav')
	smooth_despike_node.connect(datasource,'func',outputnode,'ufunc')
	smooth_despike_node.connect(datasource,'movement', outputnode, 'movement')
	
	datasource.inputs.base_directory = os.path.abspath(directory)


	
	return smooth_despike_node

"""
DIAMOND 2.0 Reward sequence

"""
def reward(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	import gold	
	
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.fsl.maths as math  
	import nipype.interfaces.afni as afni	     # afni
	conf.glm_design = '/usr/local/software/gold/linears/linear_rew.txt'


	fsl.FSLCommand.set_default_output_type('NIFTI')
	
	# define subject name
	subject = get_subject(directory)
	
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	# this is where you need to reset study specific settings
	conf.modelspec_high_pass_filter_cutoff = 60	
	conf.level1design_bases = {'hrf':{'derivs': [0,0]}}

	# define contrasts for 1st level
	contrasts = []
	contrasts.append(('RewardExpectancy','T', ['anticipationxreward_expectancy^1*bf(1)'],[1]))
	contrasts.append(('UncertainExpectancy','T', ['anticipationxuncertainty^1*bf(1)'],[1]))
	contrasts.append(('PredictionError','T', ['outcomexprediction_error^1*bf(1)'],[1]))

	# define contrasts for PPPI
	ppi_contrasts = []
	ppi_contrasts.append(('RewardExpectancy','T', ['PPI_anticipationxreward_expectancy^1'],[1]))
	ppi_contrasts.append(('PredictionError','T', ['PPI_outcomexprediction_error^1'],[1]))

	# get data input nodes from two runs
	ds1 = datasource(directory,sequence+"1")
	ds2 = datasource(directory,sequence+"2")

	
	#glm scaling for reward1
	reward1=smooth_despike_node(directory,sequence+"1",base_dir,ds1)
	
	#glm scaling for ds2
	reward2=smooth_despike_node(directory,sequence+"2",base_dir,ds2)

	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	# define a first level analysis pipeline
	l1 = gold.level1analysis(conf);
	l1.inputs.input.contrasts = contrasts

	# define a first level for PPPI
	l2 = gold.level1analysis(conf,0,"level1_pppi");
	l2.inputs.input.contrasts = contrasts

	
	# create DesignMatrix for run1
	dm1 = pe.Node(name="create_DM1",interface=Function(input_names=["matlab_function","eprime_file",'sequence'],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm1.inputs.matlab_function = "reward1_eprime2dm"
        dm1.inputs.sequence = 'reward1'
	
	# create DesignMatrix for run2
	dm2 = pe.Node(name="create_DM2",interface=Function(input_names=["matlab_function","eprime_file",'sequence'],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm2.inputs.matlab_function = "reward2_eprime2dm"
	dm2.inputs.sequence = 'reward2'
	# setup merge points
	merge_func = pe.Node(name="merge_func",interface=util.Merge(2))
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(2))
	merge_move = pe.Node(name="merge_movement",interface=util.Merge(2))
	merge_regressors = pe.Node(name="merge_regressors",interface=util.Merge(2))

	# connect behavour file from datasource to design matrix creator	
	task.connect(ds1,'behav',dm1,"eprime_file")	
	task.connect(ds2,'behav',dm2,"eprime_file")	

	# merge regressors from mCompCore
	task.connect(reward1,"output.movement",merge_regressors,"in1")	
	task.connect(reward2,"output.movement",merge_regressors,"in2")

	# merge movement files from preprocessing
	task.connect(reward1,'output.movement',merge_move,"in1")	
	task.connect(reward2,'output.movement',merge_move,"in2")
	
	# merge functional files from preprocessing
	task.connect(reward1,'output.func',merge_func,'in1')
	task.connect(reward2,'output.func',merge_func,'in2')
	
	# merge design matrix files
	task.connect(dm1,'design_matrix',merge_nDM,'in1')
	task.connect(dm2,'design_matrix',merge_nDM,'in2')

	# setup first level model 
	task.connect(merge_regressors,"out",l1,"input.movement")	
	task.connect(merge_func,'out',l1,'input.func')
	task.connect(merge_nDM,'out',l1,"input.design_matrix")
	task.connect(ds1,'mask_file',l1,"input.mask")
	
	# setup second level model for PPPI
	task.connect(merge_regressors,"out",l2,"input.movement")	
	task.connect(merge_func,'out',l2,'input.func')
	task.connect(merge_nDM,'out',l2,"input.design_matrix")
	task.connect(ds1,'mask_file',l2,"input.mask")

	# define datasink (this is where output files go)
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir

	
	# now define PPI ROIS

	pppi_rois  = [("Reward_VS",conf.ROI_VS_LR),
			  ("Reward_VLPFC",conf.ROI_leftVLPFC),
			  ("Reward_BA32",conf.ROI_BA32)]	
        

	# now do gPPI analysis
	for roi in pppi_rois:
		# setup PPPI node for a given ROI
		pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
		pppi.inputs.voi_name = roi[0]
		pppi.inputs.voi_file = roi[1]
		pppi.inputs.subject = subject
		#print(str(l2.outputs))
		task.connect(l2,'output.spm_mat_file',pppi,'spm_mat_file')
		
		# run estimate contrast on PPI results
		contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		contrast.inputs.contrasts = ppi_contrasts
		
		task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
		task.connect(pppi,'beta_images',contrast,'beta_images')
		task.connect(pppi,'residual_image',contrast,'residual_image')
		
		# deposit con images to datasing
		task.connect(contrast,'con_images',datasink,"data.pppi_"+roi[0]+"_con_images")
		task.connect(pppi,"spm_mat_file",datasink,"data.ppi_"+roi[0]+"_spm_file")
		#print(str(pppi.inputs))
		#print(str(pppi.outputs))
		

	# print and save the output of the preprocess pipeline
	gold.save_files(task,reward1.get_node('output'),datasink,["struct","mask"], not noPrint)
	gold.save_files(task,merge_move,datasink,[("out","movement")], not noPrint)
	gold.save_files(task,l1.get_node('input'),datasink,["func"], not noPrint)	
	gold.save_files(task,l1.get_node('output'),datasink,["spm_mat_file","con_images"], not noPrint)

	# create a pretty graphics
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task

"""
DIAMOND 2.0 EFNBACK sequence
"""
def efnback(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	import gold
	
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.fsl.maths as math  
	import nipype.interfaces.afni as afni	     # afni
	conf.glm_design = '/usr/local/software/gold/linears/linear_efn.txt'

	# redefine some of the parameters
	#conf.modelspec_high_pass_filter_cutoff = 60
	conf.modelspec_high_pass_filter_cutoff = 128
	conf.level1design_bases = {'hrf':{'derivs': [0,0]}}
	fsl.FSLCommand.set_default_output_type('NIFTI')
	

	subject = get_subject(directory)
	
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# some hard-coded sequence specific components
	contrasts = []	
	contrasts.append(("0back fear-0back noface","T",["zerofear*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("0back fear-0back neutral","T",["zerofear*bf(1)","zeroneutral*bf(1)"],[1,-1]))
	contrasts.append(("0back happy-0back noface","T",["zerohappy*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("0back happy-0back neutral","T",["zerohappy*bf(1)","zeroneutral*bf(1)"],[1,-1]))
	contrasts.append(("0back neutral-0back noface","T",["zeroneutral*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("2back fear-2back noface","T",["twofear*bf(1)","twoblank*bf(1)"],[1,-1]))
	contrasts.append(("2back fear-2back neutral","T",["twofear*bf(1)","twoneutral*bf(1)"],[1,-1]))
	contrasts.append(("2back happy-2back noface","T",["twohappy*bf(1)","twoblank*bf(1)"],[1,-1]))
	contrasts.append(("2back happy-2back neutral","T",["twohappy*bf(1)","twoneutral*bf(1)"],[1,-1]))
	contrasts.append(("2back neutral-2back noface","T",["twoneutral*bf(1)","twoblank*bf(1)"],[1,-1]))
	contrasts.append(("2back noface-0back noface","T",["twoblank*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("2back neutral-0back neutral","T",["twoneutral*bf(1)","zeroneutral*bf(1)"],[1,-1]))
	contrasts.append(("2back fear-0back fear","T",["twofear*bf(1)","zerofear*bf(1)"],[1,-1]))
	contrasts.append(("2back happy-0back happy","T",["twohappy*bf(1)","zerohappy*bf(1)"],[1,-1]))
	contrasts.append(("2back emotion-2back noface","T",["twofear*bf(1)","twohappy*bf(1)","twoblank*bf(1)"],[.5,.5,-1]))
	contrasts.append(("0back emotion-0back noface","T",["zerofear*bf(1)","zerohappy*bf(1)","zeroblank*bf(1)"],[.5,.5,-1]))
	contrasts.append(("2back emotion-2back neutral","T",["twofear*bf(1)","twohappy*bf(1)","twoneutral*bf(1)"],[.5,.5,-1]))
	contrasts.append(("0back emotion-0back neutral","T",["zerofear*bf(1)","zerohappy*bf(1)","zeroneutral*bf(1)"],[.5,.5,-1]))

	contrasts.append(("2back fear-0back noface","T",["twofear*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("2back happy-0back noface","T",["twohappy*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("2back neutral-0back noface","T",["twoneutral*bf(1)","zeroblank*bf(1)"],[1,-1]))
	contrasts.append(("2back emotion-0back noface","T",["twofear*bf(1)","twohappy*bf(1)","zeroblank*bf(1)"],[.5,.5,-1]))
	contrasts.append(("0back noface","T",["zeroblank*bf(1)"],[1]))
	contrasts.append(("0back neutral","T",["zeroneutral*bf(1)"],[1]))
	contrasts.append(("0back fear","T",["zerofear*bf(1)"],[1]))
	contrasts.append(("0back happy","T",["zerohappy*bf(1)"],[1]))
	contrasts.append(("2back noface","T",["twoblank*bf(1)"],[1]))
	contrasts.append(("2back neutral","T",["twoneutral*bf(1)"],[1]))
	contrasts.append(("2back fear","T",["twofear*bf(1)"],[1]))
	contrasts.append(("2back happy","T",["twohappy*bf(1)"],[1]))
	
	# get input data
	ds1 = datasource(directory,sequence+"1")
	ds2 = datasource(directory,sequence+"2")
	
	# run preprocess pipeline for both runs
	#glm_scale.inputs.out_res_name = 'residual_4d.nii.gz'
	efnback1 = smooth_despike_node(directory,sequence+"1",base_dir, ds1)
	efnback2 = smooth_despike_node(directory,sequence+"2",base_dir, ds2)
	
	
	# define first level pipeline
	l1 = gold.level1analysis(conf)
	l1.inputs.input.contrasts = contrasts
	
	# create DesignMatrix
	dm1 = pe.Node(name="create_DM1",interface=Function(input_names=["matlab_function","eprime_file", "sequence"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm1.inputs.matlab_function = "efnback1_eprime2dm"
	dm1.inputs.sequence = 'efnback1'

	dm2 = pe.Node(name="create_DM2",interface=Function(input_names=["matlab_function","eprime_file", "sequence"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm2.inputs.matlab_function = "efnback2_eprime2dm"
	dm2.inputs.sequence = 'efnback2'

	# setup merge points
	merge_func = pe.Node(name="merge_func",interface=util.Merge(2))
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(2))
	merge_move = pe.Node(name="merge_movement",interface=util.Merge(2))
	merge_regressors = pe.Node(name="merge_regressors",interface=util.Merge(2))
	
	
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	
	# merge regressors
	task.connect(efnback1,"output.movement",merge_regressors,"in1")	
	task.connect(efnback2,"output.movement",merge_regressors,"in2")	

	# merge results of both runs into one
	task.connect(ds1,'behav',dm1,"eprime_file")	
	task.connect(ds2,'behav',dm2,"eprime_file")	
	task.connect(efnback1,'output.func',merge_func,'in1')
	task.connect(efnback2,'output.func',merge_func,'in2')
	task.connect(efnback1,'output.movement',merge_move,'in1')
	task.connect(efnback2,'output.movement',merge_move,'in2')
	task.connect(dm1,'design_matrix',merge_nDM,'in1')
	task.connect(dm2,'design_matrix',merge_nDM,'in2')

	# setup first level
	task.connect(merge_regressors,"out",l1,"input.movement")	
	task.connect(merge_func,'out',l1,'input.func')
	task.connect(merge_nDM,'out',l1,"input.design_matrix")
	task.connect(ds1,'mask_file',l1,"input.mask")
	
	
	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir

	# print and save the output of the preprocess pipeline
	gold.save_files(task,efnback1.get_node('output'),datasink,["struct","mask"], not noPrint)
	gold.save_files(task,merge_move,datasink,[("out","movement")], not noPrint)
	gold.save_files(task,l1.get_node('input'),datasink,["func"], not noPrint)	
	gold.save_files(task,l1.get_node('output'),datasink,["spm_mat_file","con_images"], not noPrint)	


	# now define PPI
	ppi_contrasts = []	
	ppi_contrasts.append(("2back emotion-2back noface","T",["PPI_twofear","PPI_twohappy","PPI_twoblank"],[.5,.5,-1]))
	pppi_rois  = 	[("bilateral_amygdala",conf.ROI_amygdala_LR)]	


	# now do gPPI analysis
	for roi in pppi_rois:
		# define PPI for each ROI
		pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
		pppi.inputs.voi_name = roi[0]
		pppi.inputs.voi_file = roi[1]
		pppi.inputs.subject = subject
		task.connect(l1,'output.spm_mat_file',pppi,'spm_mat_file')
		
		# run estimate contrast for each ROI
		contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		contrast.inputs.contrasts = ppi_contrasts
		task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
		task.connect(pppi,'beta_images',contrast,'beta_images')
		task.connect(pppi,'residual_image',contrast,'residual_image')
		
		# extract con images to datasing
		task.connect(contrast,'con_images',datasink,"data.pppi_"+roi[0]+"_con_images")
		task.connect(pppi,"spm_mat_file",datasink,"data.ppi_"+roi[0]+"_spm_file")	
	# generate a pretty pipeline graphics	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')

	return task

"""
GOLD 2.0 dynamic faces sequence
"""
def dynamic_faces(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	import gold


	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.fsl.maths as math  
	import nipype.interfaces.afni as afni	     # afni
	
	# predefine hard-coded parameters
	conf.modelspec_high_pass_filter_cutoff = 256
	conf.level1design_bases = {'hrf':{'derivs': [1,0]}} #TODO check
	conf.glm_design = '/usr/local/software/gold/linears/linear_dyn.txt'
	fsl.FSLCommand.set_default_output_type('NIFTI')

	subject = get_subject(directory)




	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# some hard-coded sequence specific components
	contrasts = []	
	contrasts.append(("Anger","T",["Anger*bf(1)"],[1]))
	contrasts.append(("Anger_td","T",["Anger*bf(2)"],[1]))
	contrasts.append(("Fear","T",["Fear*bf(1)"],[1]))
	contrasts.append(("Fear_td","T",["Fear*bf(2)"],[1]))
	contrasts.append(("Sad","T",["Sad*bf(1)"],[1]))
	contrasts.append(("Sad_td","T",["Sad*bf(2)"],[1]))
	contrasts.append(("Happy","T",["Happy*bf(1)"],[1]))
	contrasts.append(("Happy_td","T",["Happy*bf(2)"],[1]))
	contrasts.append(("IDmorph","T",["IDmorph*bf(1)"],[1]))
	contrasts.append(("IDmorph_td","T",["IDmorph*bf(2)"],[1]))
	contrasts.append(("Shape","T",["Shape*bf(1)"],[1]))
	contrasts.append(("Shape_td","T",["Shape*bf(2)"],[1]))
	contrasts.append(("Anger > Shape","T",["Anger*bf(1)","Shape*bf(1)"],[1,-1]))
	contrasts.append(("Fear > Shape","T",["Fear*bf(1)","Shape*bf(1)"],[1,-1]))
	contrasts.append(("Sad > Shape","T",["Sad*bf(1)","Shape*bf(1)"],[1,-1]))
	contrasts.append(("Happy > Shape","T",["Happy*bf(1)","Shape*bf(1)"],[1,-1]))
	contrasts.append(("Emotion > Shape","T",["Anger*bf(1)","Fear*bf(1)","Sad*bf(1)","Happy*bf(1)","Shape*bf(1)"],[.25,.25,.25,.25,-1]))
	contrasts.append(("Maya_Contrast","T",["Anger*bf(1)","Fear*bf(1)","Sad*bf(1)","Happy*bf(1)"],[1,0,0,1]))

	# get datasoruce
	ds = datasource(directory,sequence)

	dynfaces = smooth_despike_node(directory, sequence, base_dir, ds)	

	# get first level analysis
	l1 = gold.level1analysis(conf)
	l1.inputs.input.contrasts = contrasts

	# mCompCor
	#cc = pe.Node(interface=wrap.mCompCor(), name="mCompCor")
	#cc.inputs.white_mask = conf.ROI_white
	

	# create DesignMatrix
	dm = pe.Node(name="create_DM",interface=Function(input_names=["matlab_function","eprime_file",'sequence'],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm.inputs.matlab_function = "dynfaces_task2dm"
	dm.inputs.sequence = 'dynface'

	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	#task.connect([(ds,cc,[('func','source'),('mask_file','brain_mask'),('movement','movement')])])

	#task.connect(dynfaces,"output.movement",l1,"input.movement")
	task.connect(dynfaces,'output.func',l1,'input.func')
	task.connect(ds,"movement",l1,"input.movement")	
	task.connect(ds,'behav',dm,"eprime_file")	
	task.connect(dm,'design_matrix',l1,'input.design_matrix')
	task.connect(ds,'mask_file',l1,"input.mask")
	
	# define datasink
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir

	# print and save the output of the preprocess pipeline
	gold.save_files(task,dynfaces.get_node('output'),datasink,["func","movement","struct","mask"], not noPrint)
	#gold.save_files(task,l1.get_node('input'),datasink,["func"], not noPrint)	
	gold.save_files(task,l1.get_node('output'),datasink,["spm_mat_file","con_images"],not noPrint)	
	
	# now define PPI
	ppi_contrasts = []	
	ppi_contrasts.append(("Anger","T",["PPI_Anger"],[1]))
	ppi_contrasts.append(("Fear","T",["PPI_Fear"],[1]))
	ppi_contrasts.append(("Sad","T",["PPI_Sad"],[1]))
	ppi_contrasts.append(("Happy","T",["PPI_Happy"],[1]))
	ppi_contrasts.append(("Shape","T",["PPI_Shape"],[1]))
	ppi_contrasts.append(("Emotion > Shape","T",["PPI_Anger","PPI_Fear","PPI_Sad","PPI_Happy","PPI_Shape"],[.25,.25,.25,.25,-1]))
	ppi_contrasts.append(("EmotionNeg > Shape","T",["PPI_Anger","PPI_Fear","PPI_Sad","PPI_Shape"],[.33,.33,.33,-1]))
	ppi_contrasts.append(("EmotionNeg > Happy","T",["PPI_Anger","PPI_Fear","PPI_Sad","PPI_Happy"],[.33,.33,.33,-1]))
	#ppi_contrasts.append(("Emotion > Shape","T",["Anger(1)","Fear(1)","Sad(1)","Happy(1)","Shape(1)"],[.25,.25,.25,.25,-1]))
	ppi_contrasts.append(("Anger_td","T",["PPI_Anger"],[1]))
	ppi_contrasts.append(("Maya_Contrast","T",["PPI_Anger","PPI_Fear","PPI_Sad","PPI_Happy"],[1,0,0,1]))
	
        pppi_rois  =    [("left_amygdala",conf.ROI_L_amyg),
                         ("right_amygdala",conf.ROI_R_amyg),
                         ("left_VLPFC",conf.ROI_L_VLPFC),
                         ("right_VLPFC",conf.ROI_R_VLPFC),
                         ("beckmann_region_1",conf.ROI_BR1),
                         ("beckmann_region_2",conf.ROI_BR2),
                         ("beckmann_region_3",conf.ROI_BR3),
                         ("beckmann_region_4",conf.ROI_BR4),
                         ("bilateral_amygdala",conf.ROI_amygdala_LR),
                         ("ROI_BilateralFG1", conf.ROI_Bilateral_FG1),
                         ("ROI_BilateralFG2", conf.ROI_Bilateral_FG2),
                         ("ROI_BilateralFG3", conf.ROI_Bilateral_FG3),
                         ("ROI_BilateralFG4", conf.ROI_Bilateral_FG4),
                         ("BilateralFusiform", conf.ROI_Bilateral_FG),
                        ("putamen",conf.ROI_putamen),
                        ("insula",conf.ROI_insula),
                        ("dlPFC",conf.ROI_dlPFC),
                        ("caudate",conf.ROI_caudate)]


	# now do gPPI analysis
	for roi in pppi_rois:
		# define PPPI for each ROI
		pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
		pppi.inputs.voi_name = roi[0]
		pppi.inputs.voi_file = roi[1]
		pppi.inputs.subject = subject
		task.connect(l1,'output.spm_mat_file',pppi,'spm_mat_file')
		
		# estimate contrast for each ROI
		contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		contrast.inputs.contrasts = ppi_contrasts
		task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
		task.connect(pppi,'beta_images',contrast,'beta_images')
		task.connect(pppi,'residual_image',contrast,'residual_image')
		
		# output con images to datasink
		task.connect(contrast,'con_images',datasink,"data.pppi_"+roi[0]+"_con_images")
		task.connect(pppi,"spm_mat_file",datasink,"data.ppi_"+roi[0]+"_spm_file")
	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task

def resting_state(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import custom_interfaces as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio 

	subject = get_subject(directory)




	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	
	# get datasoruce
	ds = datasource(directory,"restingstate")
	
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir

	rest = pe.Node(interface=wrap.Restingstate(), name='1stlevels')
	rest.inputs.output_file=subject+'output.csv'
	rest.inputs.input_wm_mask = '/usr/local/software/ITAlics/data/white_matter_mask_mni.nii'

	task.connect(ds,'mask', rest, 'input_mask')
	task.connect(ds,'func', rest, 'input_image')
	task.connect(ds,'confounds', rest, 'input_confound_csv')

	return task
	
	

# check sequence (should we process it or not)
def check_sequence(opt_list,directory,seq):
	seq_dir = seq
	

	# if no sequence specified, then run everything
	if len(opt_list) == 0 :
		return True

	# else if sequence specified, do it and ignore failed condition	
	elif "-"+seq in opt_list:
		return True

	return False

########################################################################################
#
# run pipeline if used as standalone script
#
########################################################################################
if __name__ == "__main__":	
	opts = "[-dynamic_faces|-efnback|-reward|-resting_state]"
	opt_list = []
	
	# get arguments
	if len(sys.argv) < 2:
		print "Usage: diamond.py "+opts+" <diamond subject directory>"
		sys.exit(1)
	
	# logging verbosity
	import time
	import nipype
	import logging
	from nipype import config
	import nipype.interfaces.matlab as mlab 
	
	# pick dataset that we'll be wroking on
	for arg in sys.argv:
		 if arg in opts:
		 	opt_list.append(arg)
	directory = sys.argv[len(sys.argv)-1]+"/"
        gunzip_files(directory+"ses-1/func/")
	gunzip_files(directory+"ses-1/anat/")
	
	# setup logging, display and other config
	#disp = os.environ['DISPLAY']
	log_dir = os.path.abspath(directory)+"/logs"
	if not os.path.exists(log_dir):
		os.mkdir(log_dir)
	cfg = dict(logging={'interface_level':'INFO',
						'workflow_level':'INFO',
						'log_to_file': True,
						'log_directory': log_dir},
    		   execution={'stop_on_first_crash': True,
                      	  'hash_method': 'timestamp',
                      	  'keep_inputs': True,
                      	  'remove_unnecessary_outputs': False})
	config.update_config(cfg)
	nipype.logging.update_logging(config)
	
	# setup matlab env
	mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodisplay -nosplash")
	mlab.MatlabCommand.set_default_terminal_output('stream')
	#mlab.MatlabCommand.set_default_paths(bin_dir)


	if check_sequence(opt_list,directory,"reward"):
		#log.info("\n\nREWARD pipeline ...\n\n")
		t = time.time()
		createMotionFile(directory, 'reward1')
		createMotionFile(directory, 'reward2')		
		reward = reward(directory,"reward")
		#reward.run()		
		reward.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
                folder=os.path.join(directory, 'analysis','reward')
                #gzip_files(folder)
		#log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"efnback"):
		#log.info("\n\nEFNBACK pipeline ...\n\n")
		t = time.time()		
		createMotionFile(directory, 'efnback1')
		createMotionFile(directory, 'efnback2')
		efnback = efnback(directory,"efnback")
		#efnback.run()		
		efnback.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
                folder=os.path.join(directory, 'analysis','efnback')
                gzip_files(folder)
		#log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
	if check_sequence(opt_list,directory,"dynamic_faces"):
		#log.info("\n\nDynamic_Faces pipeline ...\n\n")
		t = time.time()	
		createMotionFile(directory, 'dynface')	
		df = dynamic_faces(directory,"dynface")
		#df.run()
                #folder=os.path.join(directory, 'analysis','dynface')
                #gzip_files(folder)
		df.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		#log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"resting_state"):
		t = time.time()	
		df = resting_state(directory,"restingstate")
		df.run()
                #folder=os.path.join(directory, 'analysis','restingstate')
                #gzip_files(folder)
		#df.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		#log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

		
	#log.info("\n\npipeline complete\n\n")
