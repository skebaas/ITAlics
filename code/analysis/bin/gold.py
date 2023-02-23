#!/usr/bin/env python
# EMBARC 2.0 Analysis Pipeline
# Author: Eugene Tseytlin, Henry Chase (University of Pittsburgh)
#
import os   


"""
params is the structure with parameters that is study specific
"""
class Config:
	def __init__(self):
		import sys
		import os                                  
		import re
		
		self.CPU_CORES = 16
		self.time_repetition  = 1.5
		
		self.bet_mask = True
		self.bet_frac = 0.5 #0.6
		self.bet_robust = True
		self.bet_vertical_gradient = 0
		
		self.struct_bet_mask = True
		self.struct_bet_frac = 0.5 #0.6
		self.struct_bet_robust = True
		self.struct_bet_vertical_gradient = 0

		
		self.flirt_cost = 'mutualinfo'
		self.flirt_bins = 256
		self.flirt_dof = 12
		self.flirt_interp = 'trilinear'
		self.flirt_searchr_x = [-180, 180]
		self.flirt_searchr_y = [-180, 180]
		self.flirt_searchr_z = [-180, 180]

		self.coregister_cost_function = "nmi"
		self.coregister_separation = [4, 2]
		self.coregister_tolerance = [0.02, 0.02, 0.02, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001]
		self.coregister_fwhm = [7, 7]
		self.prepare_fieldmap_scanner = "SIEMENS"
		self.prepare_fieldmap_delta_TE = 2.46  #TODO: CHECK VALUE - similar for Encore??
	
		self.fugue_dwell_time = 0.000779983	
		self.fugue_poly_order = 3

		self.dartel_fwhm = 6 #TODO Check value
		self.dartel_voxel_size =  (2, 2, 2)
		self.susan_brightness_threshold = 750 #200.0
		self.susan_fwhm = 6
		
		self.modelspec_concatenate_runs   = False
		self.modelspec_high_pass_filter_cutoff = 60 # reward only
		self.modelspec_input_units = 'secs'
		
		self.level1design_bases = {'hrf':{'derivs': [0,0]}}
		self.level1design_timing_units = 'secs'
		self.level1estimate_estimation_method = {'Classical' : 1}
		self.contrastestimate_use_derivs = True
		self.level1design_microtime_onset = 1
		self.level1design_microtime_resolution = 16
		self.level1design_model_serial_correlations = 'FAST'

		bin_dir = os.path.dirname(os.path.realpath(__file__))
		data_dir = os.environ.get("SUPPORT_DATA_DIR",bin_dir+"/../data")

		#self.OASIS_template = data_dir+"/templates/OASIS-30_Atropos_template_in_MNI152_2mm.nii.gz"
		#self.OASIS_labels = data_dir+"/templates/OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_MNI152_2mm.nii.gz"
		#TODO: replace with normal values		
		self.OASIS_template = data_dir+"/templates/MNI_SPM_grey.nii.gz" 
		self.OASIS_labels = data_dir+"/templates/ROI_MNI_V4.nii.gz"



		ROI_dir = data_dir+"/MNI_rois/"
		ROI_suffix = "_mni.nii"

		#ROI_dir = data_dir+"/OASIS_rois/"
		#ROI_suffix = "_oasis.nii"

		self.FNIRT_config = os.getenv("FSLDIR")+"/etc/flirtsch/T1_2_MNI152_2mm.cnf"
		self.subject_site=""
		########################
		## predefine all ROIs ##
		########################
		self.ROI_white = ROI_dir+"white_matter_mask"+ROI_suffix
		self.ROI_L_insula = ROI_dir+"left_insula"+ROI_suffix
		self.ROI_R_insula = ROI_dir+"right_insula"+ROI_suffix
		self.ROI_L_amyg = ROI_dir+"left_amygdala"+ROI_suffix
		self.ROI_R_amyg = ROI_dir+"right_amygdala"+ROI_suffix
		self.ROI_VS_L = ROI_dir+"left_accumbens_area"+ROI_suffix
		self.ROI_VS_R = ROI_dir+"right_accumbens_area"+ROI_suffix
		self.ROI_VS_LR = ROI_dir+"bilateral_accumbens_area"+ROI_suffix
		self.ROI_BA9_L = ROI_dir+"left_rostral_middle_frontal"+ROI_suffix
		self.ROI_BA9_R = ROI_dir+"right_rostral_middle_frontal"+ROI_suffix
		self.ROI_BA10 = ROI_dir+"medial_brodmann_area_10"+ROI_suffix
		self.ROI_BR1 = ROI_dir+"beckmann_region_1"+ROI_suffix
		self.ROI_BR2 = ROI_dir+"beckmann_region_2"+ROI_suffix
		self.ROI_BR3 = ROI_dir+"beckmann_region_3"+ROI_suffix
		self.ROI_BR4 = ROI_dir+"beckmann_region_4"+ROI_suffix
		self.ROI_BR9 = ROI_dir+"beckmann_region_9"+ROI_suffix 
		self.ROI_PG_ACC = ROI_dir+"pregenual_ACC"+ROI_suffix
		self.ROI_D_ACC = ROI_dir+"dorsal_ACC"+ROI_suffix
		self.ROI_SG_ACC = ROI_dir+"subgenual_ACC"+ROI_suffix
		self.ROI_L_MFG_compensate = ROI_dir+"left_anterior_MFG_compensate"+ROI_suffix
		self.ROI_R_MFG_compensate = ROI_dir+"right_anterior_MFG_compensate"+ROI_suffix
		self.ROI_L_VLPFC = ROI_dir+"left_brodmann_area_47"+ROI_suffix
		self.ROI_R_VLPFC = ROI_dir+"right_brodmann_area_47"+ROI_suffix
		self.ROI_L_ant_insula =	ROI_dir+"left_anterior_insula"+ROI_suffix
		self.ROI_R_ant_insula = ROI_dir+"right_anterior_insula"+ROI_suffix
		self.ROI_putamen_L = ROI_dir+"left_putamen"+ROI_suffix
		self.ROI_putamen_R = ROI_dir+"right_putamen"+ROI_suffix
		#self.ROI_caudate_body_L = ROI_dir+"left_caudate"+ROI_suffix
		#self.ROI_caudate_body_R = ROI_dir+"right_caudate"+ROI_suffix
		self.ROI_caudate_head_L = ROI_dir+"left_caudate"+ROI_suffix
		self.ROI_caudate_head_R = ROI_dir+"right_caudate"+ROI_suffix
		self.ROI_amygdala_LR = ROI_dir+"bilateral_amygdala"+ROI_suffix
		self.ROI_leftVLPFC = ROI_dir+"leftVLPFC_2mm"+ROI_suffix
		self.ROI_BA32 = ROI_dir+"brodmann_area_32"+ROI_suffix
		self.ROI_BA24_32 = ROI_dir+"BA24_32"+ROI_suffix 
		self.ROI_BA24 = ROI_dir+"BA24"+ROI_suffix
		self.ROI_BA32Left = ROI_dir+"BA32Left"+ROI_suffix 
		self.ROI_BA25 = ROI_dir+"BA25"+ROI_suffix  
		self.ROI_thalamus = "/data/Thal_ROI.nii"

"""
load Matlab Design Matrix (nDM file)
mat_file -> nDM SPM design matrix file
trim     -> trim design matrix by number of columns

"""
def load_design_matrix(mat_file,trim=0):
	import os
	import re
	import scipy.io as sp
	import glob as gl
	import numpy
	import nipype.interfaces.matlab as mlab 
	from nipype.interfaces.base import Bunch

	# convert numpy data array
	def convert_numpy(ar,toString = False):
		lst = []
		for a in ar.tolist():
			if toString:
				lst.append(str(a))
			elif type(a) == numpy.ndarray:
				if a.size > 1:			
					lst.append(a.tolist())
				else:
					lst.append([a.tolist()])
			else:
				lst.append([a])
		return lst

	# if list of mat_files, then do a list
	mat_files = []	
	if isinstance(mat_file,list):
		mat_files = mat_file
	else:
		mat_files = [mat_file]
		
	bunches = []

	# go over mat files
	for mat_file in mat_files:
		# load design matrix 
		dm = sp.loadmat(mat_file,squeeze_me=True)

		names  = convert_numpy(dm.get('names'),True)
		onsets = convert_numpy(dm.get('onsets'))
		durations = convert_numpy(dm.get('durations'))
		
		# load up values and convert them
		# for PPPI remove last column for reward PPI; for ert last 3 columns
		# error, posterror, misc
		if trim > 0:
			names = names[0:(len(names)-trim)]
			durations = durations[0:(len(durations)-trim)]
			onsets = onsets[0:(len(onsets)-trim)]

		# create bunch to return
		bunch = Bunch(conditions=names,onsets=onsets,durations=durations)
		if 'pmod' in dm:
			pmod = []
			for i in range(0,len(dm.get('pmod'))):
				if isinstance(dm['pmod']['name'][i],unicode) or dm['pmod']['name'][i].size == 1:
					name = str(dm['pmod']['name'][i])
					param = dm['pmod']['param'][i].tolist()
					poly = dm['pmod']['poly'][i]
					pmod.append(Bunch(name=[name],param=[param],poly=[poly]))
				elif dm['pmod']['name'][i].size >  1:
					names = []
					params = []
					polys = []
					for j in range(0,dm['pmod']['name'][i].size):
						names.append(str(dm['pmod']['name'][i][j]))
						params.append(dm['pmod']['param'][i][j].tolist())
						polys.append(dm['pmod']['poly'][i][j])
					pmod.append(Bunch(name=names,param=params,poly=polys))
				else:
					pmod.append(None)
			bunch.pmod = pmod
		print(bunch)
		bunches.append(bunch)
	return bunches


"""
 Generic method to convert an E-Prime file or other tasks to nDM 
 design matrix
 Added in option for sequence for better naming
"""
def create_design_matrix(matlab_function, eprime_file, sequence):
	import os
	import re
	import glob as gl
	import nipype.interfaces.matlab as mlab 
	import chardet

	# get nDM file if available already
	mat = os.path.join(os.path.dirname(eprime_file),'nDM*'+sequence+'*.mat')
	mat = gl.glob(mat)
	#print('mat is ' + mat)
	if len(mat) > 0:
		return mat[0]
        
	#Check encoding of eprime_file
	#If UTF-16 encoding, then change it to ASCII (MATLAB-2015 cannot read UTF-16 encoding)
	with open(eprime_file, 'r') as f:
		contents = f.read()
		firstline_f = list(contents)[1]  #only need the first line to determine encoding
		detect_f = chardet.detect(firstline_f)
		print(detect_f)
	if detect_f['encoding'] != 'ascii':
		decoded = contents.decode('utf-16')
		encoded = decoded.encode('ascii')
		with open(eprime_file, 'w') as f:
			f.seek(0, os.SEEK_SET)
			f.truncate()
			f.write(encoded)
			del contents
	# execute matlab script to generate nDM file
	m = mlab.MatlabCommand()
	m.inputs.mfile = False
	m.inputs.script = matlab_function+"(\'"+eprime_file+"\');"
	m.run();

	# get nDM file (again)
	mat = os.path.join(os.path.dirname(eprime_file), 'nDM*'+sequence+'*.mat')
	mat = gl.glob(mat)
	if len(mat) > 0:
		mat = mat[0]

	return mat

"""
EMBARC 1.5 Level1 Analysis Pipeline
"""
def level1analysis(config,trim=0,name='level1'):
	import nipype.interfaces.spm as spm          # spm
	import nipype.interfaces.matlab as matlab      # how to run matlab
	import nipype.interfaces.utility as util     # utility
	import nipype.pipeline.engine as pe          # pypeline engine
	import nipype.algorithms.modelgen as model   # model specification
	
	#mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")
	l1analysis = pe.Workflow(name=name)
	inputnode = pe.Node(interface=util.IdentityInterface(fields=['movement','func','design_matrix','contrasts','mask']),name='input')

	
	# specify design matrix model
	modelspec = pe.Node(interface=model.SpecifySPMModel(), name= "modelspec")
	modelspec.inputs.concatenate_runs   = config.modelspec_concatenate_runs
	modelspec.inputs.time_repetition = config.time_repetition
	modelspec.inputs.high_pass_filter_cutoff = config.modelspec_high_pass_filter_cutoff
	modelspec.inputs.input_units = config.modelspec_input_units
	l1analysis.connect(inputnode,'movement',modelspec,'realignment_parameters')
	l1analysis.connect(inputnode,'func',modelspec,'functional_runs')
	l1analysis.connect(inputnode,('design_matrix',load_design_matrix,trim),modelspec,'subject_info')

	# create design matrix
	level1design = pe.Node(interface=spm.Level1Design(), name= "level1design")
	level1design.inputs.bases = config.level1design_bases
	level1design.inputs.timing_units = config.level1design_timing_units
	level1design.inputs.interscan_interval = config.time_repetition
	level1design.inputs.microtime_onset = config.level1design_microtime_onset
	level1design.inputs.microtime_resolution = config.level1design_microtime_resolution
	level1design.inputs.model_serial_correlations = config.level1design_model_serial_correlations 


	#Incorporate fmriprep mask in spm model
	l1analysis.connect(inputnode,'mask',level1design,'mask_image')
	l1analysis.connect(modelspec,'session_info',level1design,'session_info')

	# level 1 estimate
	level1estimate = pe.Node(interface=spm.EstimateModel(), name="level1estimate")
	level1estimate.inputs.estimation_method = config.level1estimate_estimation_method
	l1analysis.connect(level1design,'spm_mat_file',level1estimate,'spm_mat_file')	

	# no need for contrast for pppi model
	contrastestimate = pe.Node(interface = spm.EstimateContrast(), name="contrastestimate")
	contrastestimate.inputs.use_derivs = config.contrastestimate_use_derivs
	
	l1analysis.connect(inputnode,'contrasts',contrastestimate,'contrasts')
	l1analysis.connect(level1estimate,'spm_mat_file',contrastestimate,'spm_mat_file')
	l1analysis.connect(level1estimate,'beta_images',contrastestimate,'beta_images'),
	l1analysis.connect(level1estimate,'residual_image',contrastestimate,'residual_image')

	# output
	outputnode = pe.Node(interface=util.IdentityInterface(fields=['spm_mat_file','con_images','spmT_images','residual_image']),name='output')
	l1analysis.connect(contrastestimate,'spm_mat_file',outputnode,'spm_mat_file')
	l1analysis.connect(contrastestimate,'con_images',outputnode,'con_images')
	l1analysis.connect(contrastestimate,'spmT_images',outputnode,'spmT_images')

	              		              			

	return l1analysis		             



"""
EMBARC 1.0 Task Sequence Ex: ert/reward
create generic task analysis
pppi_trim_dm - trim the last N columns of design matrix for PPPI analysis
pppi_rois    - list of tuples (roi_name,roi_file) to create a workflow for 

"""
def task(pppi_trim_dm, pppi_rois):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	import nipype.interfaces.spm as spm          # spm
	from nipype.interfaces.utility import Function
	import nipype.algorithms.misc as misc
	import nipype.interfaces.utility as util     # utility
	import nipype.interfaces.io as nio           # Data i/o
	
		
	# connect components into a pipeline
	task = pe.Workflow(name="task")
	l1 = level1analysis();
	l2 = level1analysis(pppi_trim_dm,"level1_pppi");
	l2.name = "level1_pppi"

	fields=['subject','func','movement','design_matrix','contrasts','pppi_contrasts']
	inputnode = pe.Node(interface=util.IdentityInterface(fields=fields),name='input')
	task.connect([(inputnode,l1,[('func','input.func'),('design_matrix','input.design_matrix'),
				     ('movement','input.movement'),('contrasts','input.contrasts')])])
	task.connect([(inputnode,l2,[('func','input.func'),('design_matrix','input.design_matrix'),
				     ('movement','input.movement'),('contrasts','input.contrasts')])])

	fields = ['spm_mat_file','con_images']	
	for roi in pppi_rois:
		fields.append("pppi_"+roi[0]+"_con_images")

	outputnode = pe.Node(interface=util.IdentityInterface(fields=fields),name='output')	
	task.connect(l1,"output.con_images",outputnode,"con_images")
	task.connect(l1,"output.spm_mat_file",outputnode,"spm_mat_file")


	# now do gPPI analysis
	for roi in pppi_rois:
		pppi = pe.Node(interface=wrap.PPPI(), name="pppi_"+roi[0])
		pppi.inputs.voi_name = roi[0]
		pppi.inputs.voi_file = roi[1]
		task.connect(l2,'output.spm_mat_file',pppi,'spm_mat_file')
		task.connect(inputnode,'subject',pppi,'subject')
		
		contrast = pe.Node(interface = spm.EstimateContrast(), name="contrast"+roi[0])
		task.connect(inputnode,'pppi_contrasts',contrast,'contrasts')
		task.connect(pppi,'spm_mat_file',contrast,'spm_mat_file')
		task.connect(pppi,'beta_images',contrast,'beta_images')
		task.connect(pppi,'residual_image',contrast,'residual_image')
		task.connect(contrast,'con_images',outputnode,"pppi_"+roi[0]+"_con_images")
	
	return task


"""
add printing and saving of files to a workflow
workflow - where everything will be added
node     - node from which to get output files
files	 - files that need to be printed and saved
"""
def print_save_files(workflow,node,datasink,files):
	return save_files(workflow,node,datasink,files,True)
	
"""
add printing and saving of files to a workflow
workflow - where everything will be added
node     - node from which to get output files
files	 - files that need to be printed and saved
"""
def save_files(workflow,node,datasink,connectors,doPrint):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	# go over parameters
	for param in connectors:
		input_param=param
		output_param=param		
		if isinstance(param, tuple):		
			input_param=param[0]
			output_param=param[1]		
	
		workflow.connect(node,input_param,datasink,"data."+output_param)
		if doPrint:
			prnt= pe.Node(interface=wrap.Print(), name="print_"+output_param)
			prnt.inputs.out_file = output_param+".ps"
			workflow.connect(node,input_param,prnt,"in_file")	
			workflow.connect(prnt,"out_file",datasink,"ps.@par"+output_param)

"""
extract and save a set of ROIs
"""
def extract_save_rois(workflow,node,node_param,datasink,name,task_name,task_units,roi_list):
	import nipype.pipeline.engine as pe          # pypeline engine
	import wrappers as wrap
	from nipype.interfaces.utility import Function
	
	extract = pe.Node(interface=wrap.ROIExtractor(), name="extract_"+name)
	extract.inputs.roi_images = list(zip(*roi_list)[1])
	extract.inputs.average = 'none'
	extract.inputs.interpelation = 0
	#extract.inputs.source = source_image
	
	# save CSV file
	save_csv = pe.Node(name="save_csv_"+name,
		interface=Function(input_names=["task","units","names","ext","output"],
		output_names=["csv_file"],function=wrap.save_csv))
	save_csv.inputs.task = task_name
	save_csv.inputs.units = task_units
	save_csv.inputs.names = list(zip(*roi_list)[0])
	save_csv.inputs.output = name+".csv"
	
	workflow.connect(node,node_param,extract,"source")
	workflow.connect(extract,"mat_file",save_csv,"ext")
	workflow.connect(save_csv,"csv_file",datasink,"csv.@par"+name)

	return extract

