#!/usr/bin/env python2
# DIAMOND 1.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
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
conf.level1design_model_serial_correlations = 'AR(1)'





# default value to use fieldmap in the pipeline
useFieldmap=False
noPrint = True


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
	
	# define some variables beforehand
	subject=get_subject(directory)
	
	# figure out if this is new or old naming convention file
	orig = False	
	m = gl.glob(os.path.join(directory,"anat/T1MPRAGE*[0-9].nii"))
	if len(m) > 0:
		orig = True


	outfields=['func', 'struct']

	# ignore this section as it deals with original non-structured naming convention
	# for early diamond subjects
	if orig:
		# define templates for datasource
		field_template = dict(func=sequence+"/*.img",struct="anat/T1MPRAGE*[0-9].nii")
		template_args  = dict(func=[[]],struct=[[]])                


		# add behavior file to task oriented design
		if sequence.startswith('reward'):
			field_template['behav'] = sequence+"/*reward_*_task*.txt"
			template_args['behav']  = [[]]
		elif sequence.startswith('efnback'):
			field_template['behav'] = sequence+"/EFNBACK_NewEye*.txt"
			template_args['behav']  = [[]]
		elif sequence.startswith('dynamic_faces'):
			field_template['behav'] = sequence+"/subject*.txt"
			template_args['behav']  = [[]]
	else:

		# define templates for datasource for functional and structural images
		field_template = dict(func=sequence+"/*"+sequence+".[ni][im][ig]",struct="anat/*anat.nii")
		template_args  = dict(func=[[]],struct=[[]])                


		# add behavior file to task oriented design
		if sequence.startswith('reward') or sequence.startswith('efnback') or sequence.startswith('dynamic_faces'):
			field_template['behav'] = sequence+"/*task.txt"
			template_args['behav']  = [[]]
			outfields.append('behav')
			field_template['mask_file'] = sequence+"/*_mask.nii"
			outfields.append('mask_file')
			template_args['mask_file'] = [[]]
			field_template['movement'] = sequence+"/motion.1D"
			template_args['movement'] = [[]]
			outfields.append('movement')

		# if ASL then fetch functional reference file
		if sequence.startswith('asl'):
			field_template['func_ref'] = sequence+"/*"+sequence+"_ref.img"
			template_args['func_ref']  = [[]]
			outfields.append('func_ref')

		# if fieldmaps are being used, then define them
		if useFieldmap:
			field_template['fieldmap_mag']   = "field_map/*_mag.nii"
			field_template['fieldmap_phase'] ="field_map/*_phase.nii"
			template_args['fieldmap_mag']  = [[]]
			template_args['fieldmap_phase']  = [[]]
			outfields.append('fieldmap_mag')
			outfields.append('fieldmap_phase')
			
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

def glm_node(directory,sequence,base_dir,datasource):
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
	import nipype.interfaces.afni as afni	

	glm_node = pe.Workflow(name='GLM_NODE_' + sequence)
	glm_node.base_dir = base_dir


	#Despiking

	despike = pe.Node(interface=afni.Despike(), name='despike')
	despike.inputs.outputtype = 'NIFTI'

	glm_node.connect(datasource,'func',despike,'in_file')

	#glm scaling
	glm_scale = pe.Node(interface=fsl.GLM(),name='glm_scale')
	glm_scale.inputs.design = conf.glm_design
	glm_scale.inputs.out_res_name = 'residual_4d.nii'
	glm_scale.inputs.output_type = 'NIFTI'
	glm_node.connect(despike, 'out_file', glm_scale, 'in_file')
	

	# scale image so that mean 1000/original
	scale_image = pe.Node(interface=math.MathsCommand(),name='scale_image')
	#updated to 10000 from 1000 to account for SPM scaling issues
	scale_image.inputs.args = "-add 10000"
	glm_node.connect(glm_scale,'out_res',scale_image,'in_file')
	#preproc.connect(image_mean,('out_file',create_scale_args, "-add 1000"),scale_image,'args')



	# apply mask to 4d image
	apply_mask = pe.Node(interface=fsl.ApplyMask(), name="apply_mask")
	apply_mask.inputs.output_type = 'NIFTI'
	glm_node.connect(scale_image,'out_file',apply_mask,'in_file')	
	glm_node.connect(datasource,'mask_file',apply_mask,'mask_file')
	


	# smooth image using SPM
	#smoothing = pe.Node(interface=spm.Smooth(), name="smooth") 
	#smoothing.inputs.fwhm = [6, 6, 6]
	#glm_node.connect(apply_mask,'out_file',smoothing,'in_files')
        
	smoothing = pe.Node(interface=fsl.Smooth(), name="smooth")
	smoothing.inputs.fwhm = 6.0
	smoothing.inputs.output_type = 'NIFTI'
	glm_node.connect(apply_mask,'out_file',smoothing,'in_file')
	outputnode = pe.Node(interface=util.IdentityInterface(fields=['func','mask','struct', 'movement','behav','ufunc']),name='output')

	#Use this node to access the files you need from ds1
	glm_node.connect(smoothing,'smoothed_file',outputnode,'func')
	glm_node.connect(datasource,'struct',outputnode,'struct')
	glm_node.connect(datasource,'mask_file',outputnode,'mask')
	glm_node.connect(datasource,'behav',outputnode,'behav')
	glm_node.connect(datasource,'func',outputnode,'ufunc')
	glm_node.connect(datasource,'movement', outputnode, 'movement')
	
	datasource.inputs.base_directory = os.path.abspath(directory)
	
	#added 5/18/2020 for motion parameters 
	
	#print("glm_node for " + sequence + " is: " + str(glm_node.outputs))

	
	return glm_node

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
	ds1 = datasource(directory,sequence+"_1")
	ds2 = datasource(directory,sequence+"_2")

	
	#glm scaling for reward1
	reward1=glm_node(directory,sequence+"_1",base_dir,ds1)
	
	#glm scaling for ds2
	reward2=glm_node(directory,sequence+"_2",base_dir,ds2)

	# get preprocess pipelines for each run
	#pp1 = gold.glm_node(conf,useFieldmap,"preprocess_1")
	#pp2 = gold.preprocess_mni(conf,useFieldmap,"preprocess_2")
		
	# define a first level analysis pipeline
	l1 = gold.level1analysis(conf);
	l1.inputs.input.contrasts = contrasts

	# define a first level for PPPI
	l2 = gold.level1analysis(conf,0,"level1_pppi");
	l2.inputs.input.contrasts = contrasts

	
	# create DesignMatrix for run1
	dm1 = pe.Node(name="create_DM1",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm1.inputs.matlab_function = "reward_eprime2dm"
	
	# create DesignMatrix for run2
	dm2 = pe.Node(name="create_DM2",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm2.inputs.matlab_function = "reward_eprime2dm"
	
	# setup merge points
	merge_func = pe.Node(name="merge_func",interface=util.Merge(2))
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(2))
	merge_move = pe.Node(name="merge_movement",interface=util.Merge(2))
	merge_regressors = pe.Node(name="merge_regressors",interface=util.Merge(2))
	
	# mCompCor for first run
	#cc1 = pe.Node(interface=wrap.mCompCor(), name="mCompCor1")
	#cc1.inputs.white_mask = conf.ROI_white

	# mCompCor for second run
	#cc2 = pe.Node(interface=wrap.mCompCor(), name="mCompCor2")
	#cc2.inputs.white_mask = conf.ROI_white
	
		
	# now we are building a pipeline of nodes
	# connect components into a pipeline
	
	#edited by AS 5/15/2020	
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	
	# connect data from datasource to preprocessing pipeline for each run
	# with fieldmap and without
	if useFieldmap:	
		task.connect([(ds1,pp1,[('output.func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
		task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
	#else:
		#task.connect([(ds1,pp1,[('ufunc','input.func'),('struct','input.struct')])])
		#task.connect([(ds2,pp2,[('ufunc','input.func'),('struct','input.struct')])])
	
	# connect behavour file from datasource to design matrix creator	
	task.connect(ds1,'behav',dm1,"eprime_file")	
	task.connect(ds2,'behav',dm2,"eprime_file")	

	# connect output of preprocessing for each run to mCompCor to define regressors
	#task.connect([(reward1,cc1,[('output.func','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	#task.connect([(reward2,cc2,[('output.func','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	
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
	
	# setup second level model for PPPI
	task.connect(merge_regressors,"out",l2,"input.movement")	
	task.connect(merge_func,'out',l2,'input.func')
	task.connect(merge_nDM,'out',l2,"input.design_matrix")

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
	ds1 = datasource(directory,sequence+"_1")
	ds2 = datasource(directory,sequence+"_2")
	
	# run preprocess pipeline for both runs
	#glm_scale.inputs.out_res_name = 'residual_4d.nii.gz'
	efnback1 = glm_node(directory,sequence+"_1",base_dir, ds1)
	efnback2 = glm_node(directory,sequence+"_2",base_dir, ds2)
	
	
	# define first level pipeline
	l1 = gold.level1analysis(conf);
	l1.inputs.input.contrasts = contrasts
	
	# create DesignMatrix
	dm1 = pe.Node(name="create_DM1",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm1.inputs.matlab_function = "efnback_eprime2dm"

	dm2 = pe.Node(name="create_DM2",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm2.inputs.matlab_function = "efnback_eprime2dm"

	# setup merge points
	merge_func = pe.Node(name="merge_func",interface=util.Merge(2))
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(2))
	merge_move = pe.Node(name="merge_movement",interface=util.Merge(2))
	merge_regressors = pe.Node(name="merge_regressors",interface=util.Merge(2))
	
	# mCompCor for run 1
	#cc1 = pe.Node(interface=wrap.mCompCor(), name="mCompCor1")
	#cc1.inputs.white_mask = conf.ROI_white

	# mCompCor for run 2
	#cc2 = pe.Node(interface=wrap.mCompCor(), name="mCompCor2")
	#cc2.inputs.white_mask = conf.ROI_white
	
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	
	# setup the preprocessing with fieldmap and without
	if useFieldmap:	
		task.connect([(ds1,pp1,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
		task.connect([(ds2,pp2,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
	#else:
	#task.connect([(ds1,efnback1,[('func','input.func'),('struct','input.struct')])])
	#task.connect([(ds2,efnback2,[('func','input.func'),('struct','input.struct')])])

	# pipe results of preprocessing into CompCore
	#task.connect([(efnback1,cc1,[('output.func','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	#task.connect([(efnback2,cc2,[('output.func','source'),('output.mask','brain_mask'),('output.movement','movement')])])
	
	# merge regressors
	task.connect(efnback1,"output.movement",merge_regressors,"in1")	
	task.connect(efnback2,"output.movement",merge_regressors,"in2")	
	print("merge_nDM outputs are: " + str(merge_nDM.outputs))
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

	# get datasoruce
	ds = datasource(directory,sequence)

	dynfaces = glm_node(directory, sequence, base_dir, ds)	

	# get preprocess_mni pipeline
	#pp = gold.preprocess_mni(conf,useFieldmap)

	# get first level analysis
	l1 = gold.level1analysis(conf);
	l1.inputs.input.contrasts = contrasts
	
	# mCompCor
	#cc = pe.Node(interface=wrap.mCompCor(), name="mCompCor")
	#cc.inputs.white_mask = conf.ROI_white

	# create DesignMatrix
	dm = pe.Node(name="create_DM",interface=Function(input_names=["matlab_function","eprime_file"],
					output_names=["design_matrix"],function=gold.create_design_matrix))
	dm.inputs.matlab_function = "dynfaces_task2dm"

	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir

	#task.connect([(dynfaces,cc,[('output.func','source'),('output.mask','brain_mask'),('output.movement','movement')])])

	# setup merge points
	'''
	merge_func = pe.Node(name="merge_func",interface=util.Merge(2))
	merge_nDM = pe.Node(name="merge_nDM",interface=util.Merge(2))
	merge_move = pe.Node(name="merge_movement",interface=util.Merge(2))
	merge_regressors = pe.Node(name="merge_regressors",interface=util.Merge(2))
	'''

	task.connect(dynfaces,"output.movement",l1,"input.movement")	
	task.connect(ds,'behav',dm,"eprime_file")	
	task.connect(dynfaces,'output.func',l1,'input.func')
	task.connect(dm,'design_matrix',l1,'input.design_matrix')

	
	# setup the preprocessing with fieldmap and without
	#if useFieldmap:	
	#	task.connect([(ds,dynfaces,[('func','input.func'),('struct','input.struct'),
	#		('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
	#else:
	#	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])	

	# run compCore with results of preprocessing
	
	
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
	
	pppi_rois  = 	[("left_amygdala",conf.ROI_L_amyg),
			 ("right_amygdala",conf.ROI_R_amyg),
			 ("left_VLPFC",conf.ROI_L_VLPFC),
			 ("right_VLPFC",conf.ROI_R_VLPFC),	
 			 ("beckmann_region_1",conf.ROI_BR1),
			 ("beckmann_region_2",conf.ROI_BR2),
			 ("beckmann_region_3",conf.ROI_BR3),
			 ("beckmann_region_4",conf.ROI_BR4),
			 ("bilateral_amygdala",conf.ROI_amygdala_LR)	]

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
	

	# extract ROIs for 1st level
	#gold.extract_save_rois(task,l1.get_node('output'),('con_images',subset,0),datasink,"L1_AmygdalaROIs","DynamicFaces","n/a",pppi_rois)

	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task


"""
DIAMOND 1.0 Resting Sequence Ex: resting1/resting2
directory - dataset directory
sequence  - name of the sequence
subject   - optional subject name if None, embarc subject will be derived
ds		  - DataSource node for this dataset, if None embarc will be used

"""
def resting(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import CPAC					# import CPAC nuisance
	import nipype.interfaces.fsl as fsl          # fsl
	import wrappers as wrap	
	import nipype.interfaces.afni as afni		 # afni
	from nipype.interfaces.utility import Function
	import nipype.interfaces.io as nio           # Data i/o
	import nipype.algorithms.misc as misc
	import gold	
	
	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	
	# get subjects
	subject = get_subject(directory)
	
	# hard-coded custom parameters
	conf.modelspec_high_pass_filter_cutoff = 256
	
	# define resting ROI names
	resting_roi_names = ['LeftInsula','RightInsula','LeftAmygdala',
			     'RightAmygdala','LeftVS','RightVS','LeftBA9','RightBA9',
			     'BR1','BR2','BR3','BR4','BR9', 'leftVLPFC', 'leftPut', 'rightPut', 'leftCaud', 'rightCaud', 'Thal', 'vPCC', 'dPCC']
	resting_roi_images = [conf.ROI_L_insula,conf.ROI_R_insula,conf.ROI_L_amyg,conf.ROI_R_amyg,
			conf.ROI_VS_L,conf.ROI_VS_R,conf.ROI_BA9_L,conf.ROI_BA9_R,
			conf.ROI_BR1,conf.ROI_BR2,conf.ROI_BR3,conf.ROI_BR4,conf.ROI_BR9, conf.ROI_leftVLPFC, conf.ROI_putamen_L, conf.ROI_putamen_R, conf.ROI_caudate_head_L, conf.ROI_caudate_head_R, conf.ROI_thalamus, conf.ROI_vPCC, conf.ROI_dPCC]
	
	# get dataource an preprocess workflows
	ds = datasource(directory,sequence)
	pp = gold.preprocess_mni(conf,useFieldmap)
	
	# define nuisance node
	nu = pe.Node(interface=wrap.Nuisance(), name="nuisance")
	nu.inputs.white_mask = conf.ROI_white
	nu.inputs.time_repetition = conf.time_repetition

	# defin a node to only selects certain columns 
	column_select = pe.Node(interface=wrap.ColumnSelect(),name="column_select")
	column_select.inputs.selection = "18,24"
	column_select.inputs.complement = True
	
	# define GLM node
	glm = pe.Node(interface=fsl.GLM(), name="glm")
	glm.inputs.out_res_name = "residual.4d.nii.gz"
	
	# define GLM for running without global signal
	glm_NGS = pe.Node(interface=fsl.GLM(), name="glm_NGS")
	glm_NGS.inputs.out_res_name = "residual.4d.nii.gz"

	# filter an image
	filt = pe.Node(interface=fsl.ImageMaths(), name="filter")
	filt.inputs.op_string = conf.filter_image_bptf
	filt.inputs.terminal_output = 'none'
	

	# define ALFF workflows for a set of ...
	alff = dict()
	alff_nm = []
	for hl in [[0.01,0.1],[0.001,0.009],[0.0,0.04],[0.04,0.08],[0.08,0.12],[0.12,0.16],[0.16,0.20],[0.20,0.24]]: 
		nm = "ALFF"+str(hl[1]).replace("0.","_")
		alff_nm.append(nm)		
		alff[nm] = CPAC.alff.create_alff(wf_name=nm.lower())
		alff[nm].inputs.hp_input.hp = hl[0]
		alff[nm].inputs.lp_input.lp = hl[1]
	
	# define REHO workflow
	reho = CPAC.reho.create_reho()
	reho.inputs.inputspec.cluster_size = 27
	
	# define network centrality workflow
	nc = CPAC.network_centrality.create_resting_state_graphs(wf_name='network_centrality')
	nc.inputs.inputspec.method_option=0
	nc.inputs.inputspec.weight_options=[True, True]	
	nc.inputs.inputspec.threshold_option = 1
	nc.inputs.inputspec.threshold = 0.0744 
	nc.inputs.inputspec.template = conf.OASIS_labels
	zscore =  CPAC.network_centrality.get_zscore(wf_name='z_score')

	sca = dict()
	maskave = dict()
	gunzip = dict()
	
	for mask in ["BR9","LeftVS","RightVS","BR2","BR3", "leftVLPFC", "leftPut", "rightPut", "leftCaud", "rightCaud", "Thal", "vPCC", "dPCC"]:
		sca[mask] = CPAC.sca.create_sca(name_sca="sca_"+mask);
		maskave[mask] = pe.Node(interface=afni.Maskave(),name="roi_ave_"+mask)
		maskave[mask].inputs.outputtype = "NIFTI"
		maskave[mask].inputs.quiet= True
		maskave[mask].inputs.mask = resting_roi_images[resting_roi_names.index(mask)]
		gunzip[mask] = pe.Node(interface=misc.Gunzip(),name="gunzip_"+mask)
	
	# define ROI averages
	roiave = pe.MapNode(interface=afni.Maskave(),name="roi_ave",iterfield="mask")
	roiave.inputs.outputtype = "NIFTI"
	roiave.inputs.mask = resting_roi_images
	roiave.inputs.quiet= True

	# define ROI correlation 
	corroi = pe.Node(interface=wrap.CorrelateROIs(),name="corr_roi")
	corroi.inputs.roi_names = resting_roi_names
	corroi.inputs.task_name = "Resting_State"
	corroi.inputs.out_file = subject+"_"+sequence+"_outcomes_CORR.csv"

	# define datasing
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# start the workflow
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	
	# feed input into preprocessing with fieldmap and without
	if useFieldmap:	
		task.connect([(ds,pp,[('func','input.func'),('struct','input.struct'),
			('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
	else:
		task.connect([(ds,pp,[('func','input.func'),('struct','input.struct')])])	

	# feed the output of preprocessing into nuisance
	task.connect([(pp,nu,[('output.ufunc','source'),
				 ('output.mask','brain_mask'),
				 ('output.movement','movement')])])
	
	# run GLM with all regressors from nuisance and filter its output
	task.connect(nu,"regressors",glm,"design")
	task.connect(pp,"output.func",glm,"in_file")
	task.connect(glm,"out_res",filt,"in_file")
	
	# run GLM with columsn 18 and 24 excluded (global noise)
	task.connect(nu,"regressors",column_select,"in_file")
	task.connect(column_select,"out_file",glm_NGS,"design")
	task.connect(pp,"output.func",glm_NGS,"in_file")

	# run corrolation of ROI avarages from filtered GLM image
	task.connect(filt,"out_file",roiave,"in_file")
	task.connect(roiave,"out_file",corroi,"in_files")
		
	# feed GLM without global noise into ALFF 
	for nm in alff_nm:
		#task.connect(glm,'out_res',alff[nm],'inputspec.rest_res')
		task.connect(glm_NGS,'out_res',alff[nm],'inputspec.rest_res')
		task.connect(pp,'output.mask',alff[nm],'inputspec.rest_mask')	

	# run REHO
	task.connect(filt,"out_file",reho,"inputspec.rest_res_filt")
	task.connect(pp,"output.mask",reho,"inputspec.rest_mask")
	
	# run network centrality
	task.connect(glm,'out_res',nc,'inputspec.subject')
	task.connect(nc,'outputspec.centrality_outputs',zscore,'inputspec.input_file')
	task.connect(pp,'output.mask',zscore,'inputspec.mask_file')

	for mask in ["BR9","LeftVS","RightVS","BR2","BR3", "leftVLPFC", "leftPut", "rightPut", "leftCaud", "rightCaud", "Thal", "vPCC", "dPCC"]:
		task.connect(filt,"out_file",maskave[mask],"in_file")
		task.connect(filt,"out_file",sca[mask],"inputspec.functional_file")
		task.connect(maskave[mask],"out_file",sca[mask],"inputspec.timeseries_one_d")
		task.connect(sca[mask],("outputspec.Z_score",subset,0),gunzip[mask],'in_file')
		task.connect(sca[mask],"outputspec.Z_score",datasink,"data.sca."+mask)
	
	# output results into datasink
	task.connect(reho,"outputspec.z_score",datasink,"data.reho")
	task.connect(corroi,"out_file",datasink,"csv.@par5")
	
	# alff_Z_img in 0.3.5 now in 0.3.6 falff_img	
	for nm in alff_nm:	
		task.connect(alff[nm],"outputspec.alff_Z_img",datasink,"data."+nm.lower())
	task.connect(zscore,"outputspec.z_score_img",datasink,"data.nc")
	

	# print and save the output of the preprocess pipeline
	gold.save_files(task,pp.get_node('output'),datasink,["func","movement","struct","mask"], not noPrint)	
	
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task



# this method is here as a legacy when all we want to do is preprocess the dataset and
# not run any sequence over it
def preprocess(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o
	import gold

	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# get components
	ds = datasource(directory,sequence)
	pp = gold.preprocess_mni(conf)
	#pp = gold.preprocess_mni_no_fieldmap(conf)	
	#pp.get_node("input").inputs.template = embarc.OASIS_template	
	
	# connect components into a pipeline
	task = pe.Workflow(name=sequence)
	task.base_dir = base_dir
	task.connect([(ds,pp,[('func','input.func'),('struct','input.struct'),
	('fieldmap_mag','input.fieldmap_mag'),('fieldmap_phase','input.fieldmap_phase')])])
				
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# print and save the output of the preprocess pipeline
	gold.print_save_files(task,pp.get_node('output'),datasink,("func","movement","struct","mask"))	
		
	task.write_graph(dotfilename=sequence+"-workflow")#,graph2use='flat')
	return task

"""
ASL sequence
"""
def asl(directory,sequence):
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.interfaces.io as nio           # Data i/o
	import gold
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.fsl as fsl          # fsl
	import nipype.interfaces.fsl.maths as math  
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.afni as afni	     # afni
	import nipype.workflows.fmri.spm.preprocess as dartel # preprocess
	import wrappers as wrap		
	from nipype.interfaces.nipy.preprocess import Trim

	# define base directory
	base_dir = os.path.abspath(directory+"/analysis/")
	if not os.path.exists(base_dir):
		os.makedirs(base_dir)
	out_dir = os.path.abspath(directory+"/output/"+sequence)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	# get components
	inputnode = datasource(directory,sequence)
	
	# connect components into a pipeline
	preproc = pe.Workflow(name=sequence)
	preproc.base_dir = base_dir
	
	# lets do a workflow	
	fsl.FSLCommand.set_default_output_type('NIFTI')
	
	
	# realign 4D functional
	realign = pe.Node(interface=spm.Realign(), name="realign")
	realign.inputs.register_to_mean = True
	preproc.connect(inputnode,"func",realign,"in_files")

	# skull strip mean functional image
	bet_mean = pe.Node(interface=fsl.BET(), name="bet_mean")
	bet_mean.inputs.mask = conf.bet_mask
	bet_mean.inputs.frac = conf.bet_frac
	bet_mean.inputs.robust = conf.bet_robust
	bet_mean.inputs.vertical_gradient = conf.bet_vertical_gradient
	preproc.connect(realign,'mean_image',bet_mean,'in_file') 

	# skull strip mean structural image
	bet_struct = pe.Node(interface=fsl.BET(), name="bet_struct")
	bet_struct.inputs.mask = conf.bet_mask
	bet_struct.inputs.frac = conf.bet_frac
	bet_struct.inputs.robust = conf.bet_robust
	bet_struct.inputs.vertical_gradient = conf.bet_vertical_gradient
	preproc.connect(inputnode,'struct',bet_struct,'in_file')	


	# coregister images
	coreg_func2struct = pe.Node(interface=spm.Coregister(),name="coreg_func2struct")
	coreg_func2struct.inputs.jobtype = "estimate"
	coreg_func2struct.inputs.cost_function = conf.coregister_cost_function
	coreg_func2struct.inputs.separation = conf.coregister_separation
	coreg_func2struct.inputs.tolerance = conf.coregister_tolerance
	coreg_func2struct.inputs.fwhm = conf.coregister_fwhm
	preproc.connect(bet_struct,'out_file',coreg_func2struct,'target')
	preproc.connect(realign,'realigned_files',coreg_func2struct,'apply_to_files')
	preproc.connect(bet_mean,'out_file',coreg_func2struct,'source')

	# coregister reference images
	coreg_ref2struct = pe.Node(interface=spm.Coregister(),name="coreg_ref2struct")
	coreg_ref2struct.inputs.jobtype = "estimate"
	coreg_ref2struct.inputs.cost_function = conf.coregister_cost_function
	coreg_ref2struct.inputs.separation = conf.coregister_separation
	coreg_ref2struct.inputs.tolerance = conf.coregister_tolerance
	coreg_ref2struct.inputs.fwhm = conf.coregister_fwhm
	preproc.connect(bet_struct,'out_file',coreg_ref2struct,'target')
	preproc.connect(inputnode,'func_ref',coreg_ref2struct,'apply_to_files')
	preproc.connect(inputnode,'func_ref',coreg_ref2struct,'source')


	# invoke ASL script to create CBF mean image
	pcasl = pe.Node(interface=wrap.pCASL(), name="pCASL")
	preproc.connect(coreg_func2struct,'coregistered_files',pcasl,'in_file')
	preproc.connect(coreg_ref2struct,'coregistered_files',pcasl,'ref_file')

	# convert image using fslmats
	convert_image = pe.Node(interface=math.MathsCommand(),name='convert_image')
	convert_image.inputs.args = "-mul 1"
	preproc.connect(pcasl,'cbf_image',convert_image,'in_file')
		
	
	# create dartel template
	dartel_template = dartel.create_DARTEL_template()
	dartel_template.inputs.inputspec.template_prefix = 'Template'
	preproc.connect(inputnode, 'struct',dartel_template,'inputspec.structural_files')


	# now lets do normalization with DARTEL
	norm_func =  pe.Node(interface=spm.DARTELNorm2MNI(modulate=True),name='norm_func')	
	norm_func.inputs.fwhm = conf.dartel_fwhm
	norm_func.inputs.voxel_size = conf.dartel_voxel_size
	preproc.connect(dartel_template,'outputspec.template_file',norm_func,'template_file')
	preproc.connect(dartel_template, 'outputspec.flow_fields', norm_func, 'flowfield_files')
	preproc.connect(convert_image,'out_file',norm_func,'apply_to_files')
	
	# now lets do normalization with DARTEL
	norm_struct =  pe.Node(interface=spm.DARTELNorm2MNI(modulate=True),name='norm_struct')
	norm_struct.inputs.fwhm = conf.dartel_fwhm  #TODO Check value
	preproc.connect(dartel_template,'outputspec.template_file',norm_struct,'template_file')
	preproc.connect(dartel_template, 'outputspec.flow_fields', norm_struct, 'flowfield_files')
	preproc.connect(bet_struct,'out_file',norm_struct,'apply_to_files')

	# skull strip mean structural image
	bet_cbf = pe.Node(interface=fsl.BET(), name="bet_cbf")
	bet_cbf.inputs.mask = conf.bet_mask
	bet_cbf.inputs.frac = conf.bet_frac
	bet_cbf.inputs.robust = conf.bet_robust
	bet_cbf.inputs.vertical_gradient = conf.bet_vertical_gradient
	preproc.connect(norm_func,'normalized_files',bet_cbf,'in_file')	

		
	# calculated brighness threshold for susan (mean image intensity * 0.75)
	#image_mean = pe.Node(interface=fsl.ImageStats(),name='image_mean')	
	#image_mean.inputs.op_string = "-M"
	#preproc.connect(bet_cbf,'out_file',image_mean,'in_file')


	# smooth image using SUSAN
	#susan = pe.Node(interface=fsl.SUSAN(), name="smooth")
	#susan.inputs.fwhm = conf.susan_fwhm
	#preproc.connect(norm_func,'normalized_files',susan,'in_file') 
	#preproc.connect(image_mean,('out_stat',gold.create_brightness_threshold),susan,'brightness_threshold') 
	

	# merge inputs from bet for scaling the image
	#op_merge = pe.Node(interface=util.Merge(2),name='op_merge')
	#preproc.connect(bet_cbf,'out_file',op_merge,"in1")	
	#preproc.connect(bet_cbf,"mask_file",op_merge,"in2")

	# scale image: 4D - 3D (mean image) + 1000 (within the mask)
	#scale_image = pe.Node(interface=math.MultiImageMaths(),name='scale_image')
	#scale_image.inputs.op_string = "-Tmean -mul -1 -add %s -add 1000 -mas %s"
	#preproc.connect(bet_cbf,'out_file',scale_image,'in_file')
	#preproc.connect(op_merge,"out",scale_image,'operand_files')
	
	# calculated brighness threshold for susan (mean image intensity * 0.75)
	image_sd = pe.Node(interface=fsl.ImageStats(),name='image_sd')	
	image_sd.inputs.op_string = "-S"
	#preproc.connect(scale_image,'out_file',image_sd,'in_file')
	preproc.connect(bet_cbf,'out_file',image_sd,'in_file')


	# smooth image using SUSAN
	susan = pe.Node(interface=fsl.SUSAN(), name="smooth")
	susan.inputs.brightness_threshold = conf.susan_brightness_threshold 
	susan.inputs.fwhm = conf.susan_fwhm
	#preproc.connect(scale_image,'out_file',susan,'in_file') 
	preproc.connect(bet_cbf,'out_file',susan,'in_file') 
	preproc.connect(image_sd,('out_stat',gold.create_brightness_threshold),susan,'brightness_threshold') 
	            


	outputnode = pe.Node(interface=util.IdentityInterface(fields=['func','mask','movement','struct']),name='output')
	preproc.connect(susan,'smoothed_file',outputnode,'func')
	preproc.connect(realign,'realignment_parameters',outputnode,'movement')
	preproc.connect(norm_struct,'normalized_files',outputnode,'struct')
	preproc.connect(bet_cbf,'mask_file',outputnode,'mask')

	# datasing			
	datasink = pe.Node(nio.DataSink(), name='datasink')
	datasink.inputs.base_directory = out_dir
	
	# print and save the output of the preprocess pipeline
	gold.save_files(preproc,outputnode,datasink,["func","movement","struct","mask"],not noPrint)	
		
	preproc.write_graph(dotfilename=sequence+"-workflow")
	return preproc




# check sequence (should we process it or not)
def check_sequence(opt_list,directory,seq):
	seq_dir = seq
	
	# check if directory exists
	if not (os.path.exists(directory+seq_dir) or os.path.exists(directory+seq_dir+"_1")):
		print "Error: data directory for "+seq+" does not exists, skipping .."
		print "Missing directory: "+directory+seq_dir
		return False

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
	opts = "[-dynamic_faces|-efnback|-reward|-resting_state|-asl|-fieldmap|-noprint|-trio]"
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
	
	# check directory
	if not os.path.exists(directory) and directory != "preprocess":
		print "Error: data directory "+directory+" does not exist"
		sys.exit(1)
	
	if "subject" == get_subject(directory):
		print "Error: "+directory+" is not a valid EMBARC data directory"
		sys.exit(1)
	
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
	log = nipype.logging.getLogger('workflow')
	l = nipype.logging.getLogger('interface').parent.handlers[0]
	l.setFormatter(logging.Formatter('%(name)-2s %(levelname)-2s:\t %(message)s'))
	###########
	
	# setup matlab env
	mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodisplay -nosplash")
	mlab.MatlabCommand.set_default_terminal_output('stream')
	#mlab.MatlabCommand.set_default_paths(bin_dir)
	
	if "-fieldmap" in opt_list:	
		useFieldmap = True
		opt_list.remove("-fieldmap")
	if "-noprint" in opt_list:	
		noPrint = True
		opt_list.remove("-noprint")
	
	# change dwell time based on scanner
	if "-trio" in opt_list:
		opt_list.remove("-trio")
		conf.fugue_dwell_time = 0.000779983     # For Diamond Trio (currently in gold.py)
	else:
		conf.fugue_dwell_time = 0.00064  	# For Diamond Prisma + Impress Prisma



	if check_sequence(opt_list,directory,"reward"):
		log.info("\n\nREWARD pipeline ...\n\n")
		t = time.time()		
		reward = reward(directory,"reward")
		#reward.run()		
		reward.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
	if check_sequence(opt_list,directory,"resting_state"):
		log.info("\n\nRESTING pipeline ...\n\n")
		t = time.time()		
		resting = resting(directory,"resting_state")
		resting.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		#resting.run()
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"efnback"):
		log.info("\n\nEFNBACK pipeline ...\n\n")
		t = time.time()		
		efnback = efnback(directory,"efnback")
		#efnback.run()		
		efnback.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	
	if check_sequence(opt_list,directory,"dynamic_faces"):
		log.info("\n\nDynamic_Faces pipeline ...\n\n")
		t = time.time()		
		df = dynamic_faces(directory,"dynamic_faces")
		df.run()		
		#df.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

	if check_sequence(opt_list,directory,"asl"):
		log.info("\n\nASL pipeline ...\n\n")
		t = time.time()		
		workflow = asl(directory,"asl")
		#workflow.run()		
		workflow.run(plugin='MultiProc', plugin_args={'n_procs' : conf.CPU_CORES})
		log.info("elapsed time %.03f minutes\n" % ((time.time()-t)/60))

		
	log.info("\n\npipeline complete\n\n")
