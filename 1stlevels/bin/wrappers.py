#
# Define ASL Class NODE to wrap matlab code
#
import sys
import os                                  
import re
import glob
import scipy.io as sio

from nipype.interfaces.matlab import MatlabCommand
from nipype.interfaces.base import TraitedSpec, BaseInterface, BaseInterfaceInputSpec, File, traits, InputMultiPath, StdOutCommandLine,StdOutCommandLineInputSpec
from string import Template

# Standard library imports
from copy import deepcopy

# Third-party imports
import numpy as np

# Local imports
from nipype.interfaces.base import (OutputMultiPath, TraitedSpec, isdefined,
                                    traits, InputMultiPath, File, CommandLineInputSpec,CommandLine)
from nipype.interfaces.spm.base import (SPMCommand, scans_for_fname,
                                        func_is_3d,
                                        scans_for_fnames, SPMCommandInputSpec)
from nipype.utils.filemanip import (fname_presuffix, filename_to_list,
                                    list_to_filename, split_filename)


class ASLInputSpec(BaseInterfaceInputSpec): 
	in_files = InputMultiPath(traits.List(File(exists=True)),desc="Input Functional 3D files",mandatory=True)
	first_image_type = traits.Int(desc="Is the first image Label (0) or Control (1)",default_value=0)
	TR	= traits.Int(desc="Time Repetition (TR)",default_value=4460)
	
class ASLOutputSpec(TraitedSpec):
	cbf_image = File(exists=True,desc="Mean CBF Image", mandatory=True)
	cbf_value = traits.Float(desc = "CBF Value", mandatory=True)

"""
	ASL Script Wrapper: runs asl_script matlab code
"""
class ASL(BaseInterface):
	"""
	Run matlab asl_script to comput CBF image
	Examples
	--------

	>>> import nipype.interfaces.spm as spm
	>>> asl = ASL()
	>>> asl.inputs.in_files = 'functional.nii'
	>>> asl.inputs.first_image_type = 0
	>>> asl.run() # doctest: +SKIP

	"""
	input_spec = ASLInputSpec
	output_spec = ASLOutputSpec
	
  	
	def _run_interface(self, runtime):
		from nipype.interfaces.spm.base import scans_for_fname,scans_for_fnames
		from nipype.utils.filemanip import filename_to_list,list_to_filename

		# setup parameters
		input_dir = "."
		in_files = "{"
		asl_first = str(self.inputs.first_image_type)
		TR = str(self.inputs.TR)
		# convert images to cell array string in matlab
		for f in sorted(scans_for_fnames(filename_to_list(self.inputs.in_files))):
			in_files += "'"+f+"',\n"
			input_dir = os.path.dirname(f)
		in_files = in_files[:-2]+"}"
		self.input_dir = input_dir

		d = dict(in_files=in_files,in_dir=input_dir,first_image_type=asl_first,TR =TR)
		myscript = Template("""
		warning('off','all');
		cd('$in_dir');
		input = char($in_files);
		asl_script(input,$first_image_type,0,$TR);
		exit;
		""").substitute(d)
		mlab = MatlabCommand(script=myscript,matlab_cmd="matlab -nodesktop -nosplash",mfile=True)
		result = mlab.run()
		return result.runtime

	def _list_outputs(self):
		import scipy.io as sp
		#setup subject
		input_dir = "."
		for f in scans_for_fnames(filename_to_list(self.inputs.in_files)):
			input_dir = os.path.dirname(f)
			break
		pt = re.compile(".*/embarc_CU_(.*)_mri_fmriraw.*")
		mt = pt.match(input_dir)
		if mt:
			subject = mt.group(1)

		# get the float value
		mat = sp.loadmat(os.path.abspath(input_dir+"/mean_CBF_spm_save.mat"),squeeze_me=True)
				
		outputs = self._outputs().get()
		outputs['cbf_image'] = os.path.abspath(input_dir+"/meanCBF_"+subject+".nii")
		outputs['cbf_value'] = mat['L']
		return outputs

class pCASLInputSpec(BaseInterfaceInputSpec): 
	in_file =  File(exists=True,desc="Input Functional 4D files",mandatory=True)
	ref_file = File(exists=True,desc="Input Reference file",mandatory=True)
	
class pCASLOutputSpec(TraitedSpec):
	cbf_image = File(exists=True,desc="Mean CBF Image", mandatory=True)
	#cbf_value = traits.Float(desc = "CBF Value", mandatory=True)

"""
	pCASL Script Wrapper: runs asl_script matlab code
"""
class pCASL(BaseInterface):
	"""
	Run matlab cbfmap_base_pCASL to comput CBF image
	Examples
	--------

	>>> import nipype.interfaces.spm as spm
	>>> asl = ASL()
	>>> asl.inputs.in_files = 'functional.nii'
	>>> asl.inputs.first_image_type = 0
	>>> asl.run() # doctest: +SKIP

	"""
	input_spec = pCASLInputSpec
	output_spec = pCASLOutputSpec
	
  	
	def _run_interface(self, runtime):
		# setup parameters
		in_file = "'"+str(self.inputs.in_file)+"'"
		ref_file = "'"+str(self.inputs.ref_file)+"'"
		
		d = dict(in_file=in_file,ref_file=ref_file)
		myscript = Template("""
		warning('off','all');
		cbf = cbfmap_base_pCASL($in_file,$ref_file);
		save([dirname($in_file) '/CBF_pCASL.mat'],'cbf');
		exit;
		""").substitute(d)
		mlab = MatlabCommand(script=myscript,matlab_cmd="matlab -nodesktop -nosplash",mfile=True)
		result = mlab.run()
		return result.runtime

	def _list_outputs(self):	
		import scipy.io as sp	
		outdir = os.path.dirname(self.inputs.in_file)
		#mat = sp.loadmat(os.path.abspath(outdir+"/CBF_pCASL.mat"),squeeze_me=True)

		outputs = self._outputs().get()
		outputs['cbf_image'] = os.path.join(outdir,"CBF_pCASL.img")
		#outputs['cbf_value'] = mat['cbf']
		
		return outputs
		
##############################################################################		
class PPPIInputSpec(BaseInterfaceInputSpec): 
	voi_file = InputMultiPath(File(exists=True),desc="VOI region file",field="VOI",mandatory=True)
	voi_name = traits.String(field='Region', mandatory=True,desc='VOI region name')
	subject = traits.String("subject",field='subject',desc='Subject Name')
	spm_mat_file = File(exists=True,field='directory', desc='absolute path to SPM.mat',copyfile=True,mandatory=True)
	comp_contrasts = traits.Int(0,desc="Compute Contrasts",field='CompContrasts', usedefault=True)	
class PPPIOutputSpec(TraitedSpec):
	beta_images = OutputMultiPath(File(exists=True), desc='beta images')
	residual_image = File(exists=True, desc='residual images')
	spm_mat_file = File(exists=True, desc='Updated SPM mat file')
	
"""
	PPPI Script Wrapper: runs PPPI matlab code
"""
class PPPI(BaseInterface):
	"""
	Run matlab PPPI to compute gPPI
	Examples
	--------

	>>> import nipype.interfaces.spm as spm
	>>> pppi = PPPI()
	>>> pppi.inputs.voi_file = 'BR3.nii'
	>>> pppi.inputs.voi_name = 'BR3'
	>>> pppi.inputs.subject = 'subject1'
	>>> pppi.inputs.spm_mat_file = 'SPM.mat'
	>>> pppi.run() 

	"""
	input_spec = PPPIInputSpec
	output_spec = PPPIOutputSpec
	
  	
	def _run_interface(self, runtime):
		from nipype.interfaces.spm.base import scans_for_fname,scans_for_fnames
		from nipype.utils.filemanip import filename_to_list,list_to_filename

		# setup parameters
		voi_file = str(self.inputs.voi_file)
		voi_name = "'"+str(self.inputs.voi_name)+"'"
		subject  = "'"+str(self.inputs.subject)+"'"
		spm_file = "'"+str(self.inputs.spm_mat_file)+"'"
		contrast = str(self.inputs.comp_contrasts)
		directory = os.path.dirname(re.sub("[\[\]']","",spm_file))
		
		d = dict(voi_file=voi_file,voi_name=voi_name,subject=subject,
				 spm_file=spm_file,directory=directory,contrast=contrast)
		myscript = Template("""
		warning('off','all');

		% copy input 
		spm = $spm_file;		
		load(spm);		
		SPM.VResMS.fname = [SPM.swd '/ResMS.img'];
		SPM.xVol.VRpv.fname = [SPM.swd '/RPV.img'];
		SPM.VM.fname = ls([SPM.swd '/../*/mask.img']);
		for i=1:length(SPM.Vbeta)
			if strcmp(SPM.Vbeta(i).fname(1),'/') == 0
				SPM.Vbeta(i).fname = [SPM.swd '/' SPM.Vbeta(i).fname];				
			end
		end		
		SPM.swd = '$directory';
		save(spm,'SPM');
		cd('$directory');
		%dos(['ln -sf ' path '/../*/*.[ih]* .']);
		load('ppi_master_template.mat')
		P.CompContrasts = $contrast;
		P.VOI=char($voi_file);
		P.Region=char($voi_name);
		P.subject=char($subject);
        	P.directory=char('$directory');
		P.FLmask =1;
		P.equalroi = 0;
		mat = [char($subject),'_analysis_',char($voi_name),'.mat'];
		save(mat,'P');
        	PPPI(mat);
       		exit;
		""").substitute(d)
		mlab = MatlabCommand(script=myscript, mfile=True)
		result = mlab.run()
		return result.runtime

	def _list_outputs(self):
		voi = str(self.inputs.voi_name)
		outdir = os.path.dirname(self.inputs.spm_mat_file)
		outdir = os.path.abspath(os.path.join(outdir,"PPI_"+voi))
		outputs = self._outputs().get()
		#outputs['spm_mat_file'] = self.inputs.spm_mat_file
		outputs['spm_mat_file'] = os.path.join(outdir,"SPM.mat")
		outputs['beta_images'] = sorted(glob.glob(outdir+'/beta_00*.img'))
		outputs['residual_image'] = os.path.join(outdir,"ResMS.img")
		return outputs

##############################################################################		
class ImageCalcInputSpec(SPMCommandInputSpec):
	#in_file =  File(exists=True, desc='Input Image',field='input',mandatory=True,copyFile=True)
	in_file =  InputMultiPath(File(exists=True), field='input',desc='Input Image',mandatory=True,copyfile=False)
	out_file = File(value='out.nii',desc='Output Image Name',field='output',usedefault=True, genfile=True, hash_files=False)
	out_dir =  File(value='', field='outdir', usedefault=True,desc='Output directory')
	expression = traits.String(field='expression', mandatory=True,desc='Expression for Calculation')
	dmtx = traits.Int(0, field='options.dmtx', usedefault=True,desc='DMTX')
	mask = traits.Int(0, field='options.mask', usedefault=True,desc='Mask')
	interp = traits.Int(0, field='options.interp', usedefault=True,desc='Interpalation')
	data_type = traits.Int(16, field='options.dtype', usedefault=True,desc='Datatype')

class ImageCalcOutputSpec(TraitedSpec):
	out_file = File(exists=True, desc='Output Image')

class ImageCalc(SPMCommand):
	"""Use spm_imcalc to do arbitrary arithmatic on images
	Examples
	--------

	>>> import nipype.interfaces.spm as spm
	>>> calc = ImageCalc()
	>>> calc.inputs.in_files = 'functional.nii'
	>>> calc.inputs.expression = 'i1/5'
	>>> calc.run() # doctest: +SKIP

	"""

	input_spec = ImageCalcInputSpec
	output_spec = ImageCalcOutputSpec
	
	_jobtype = 'util'
	_jobname = 'imcalc'

	def _generate_job(self, prefix='', contents=None):
		if isinstance(prefix,str):
			prefix = prefix.replace("imcalc{1}","imcalc")
		return super(ImageCalc, self)._generate_job(prefix,contents)
		
	def _format_arg(self, opt, spec, val):
		f = re.sub("[\[\]']","",str(self.inputs.in_file))
		if opt == 'out_file':
			val = "calc_"+os.path.basename(f)
			self.inputs.out_file = val
		if opt == 'out_dir':
			val = os.path.dirname(f)+"/"
			self.inputs.out_dir = val
		
		if opt == 'in_file' or opt == 'out_dir':
			if isinstance(val,list):
				return np.array(val,dtype=object)
			return np.array([val],dtype=object)
		return super(ImageCalc, self)._format_arg(opt, spec, val)

	def _list_outputs(self):
		outputs = self._outputs().get()
		out = self.inputs.out_file
		if os.path.isdir(self.inputs.out_dir):
			out = os.path.abspath(self.inputs.out_dir+"/"+out)
		outputs['out_file'] = out
		return outputs

##############################################################################		

class ROIExtractorInputSpec(SPMCommandInputSpec):
	source = InputMultiPath(File(exists=True), field='src.srcimgs', desc='Source Images',copyfile=False, mandatory=True)
	roi_images = traits.List(File(exists=True),field='roispec{*}.srcimg',desc='Lis of ROI images',mandatory=True)
	average = traits.Enum("none","vox",field="avg",usedefault=True,desc='Average')
	interpelation = traits.Int(0, field='interp', usedefault=True,desc='Interpelation')

class ROIExtractorOutputSpec(TraitedSpec):
	mat_file = File(exists=True, desc='Output Matlab File')

class ROIExtractor(SPMCommand):
	"""Use SPM Volume toolbox to extract ROIs from images
	Examples
	--------

	>>> import nipype.interfaces.spm as spm
	>>> roi = ROIExtractor()
	>>> roi.inputs.source = 'functional.nii'
	>>> roi.inputs.roi = [roi1.nii,roi2.nii']
	>>> roi.run() # doctest: +SKIP

	"""

	input_spec = ROIExtractorInputSpec
	output_spec = ROIExtractorOutputSpec
	
	_jobtype = 'tools'
	_jobname = 'vgtbx_Volumes{1}.Single_Volumes.tbxvol_extract'

	# fix the issue of {1} prefix, SPM won't run w/ it
	def _generate_job(self, prefix='', contents=None):
		if isinstance(prefix,str):
			prefix = prefix.replace("tbxvol_extract{1}","tbxvol_extract")
		return super(ROIExtractor, self)._generate_job(prefix,contents)
		
	def _format_arg(self, opt, spec, val):
		if opt == 'source' or opt == 'roi_images':
			return scans_for_fnames(val)
		return super(ROIExtractor, self)._format_arg(opt, spec, val)
	
	# overwrite parse_inputs, to make a list of ROIs
	def _parse_inputs(self, skip=()):
		spmdict = {}
		metadata = dict(field=lambda t: t is not None)
		for name, spec in self.inputs.traits(**metadata).items():
			if skip and name in skip:
				continue
			value = getattr(self.inputs, name)
			if not isdefined(value):
				continue
			field = spec.field
			if '{*}' in field and isinstance(value,list):
				for i,v in enumerate(value):
					spmdict[field.replace("{*}","{"+str(i+1)+"}")] = self._format_arg(name, spec,[v])	
			elif '.' in field:
				fields = field.split('.')
				dictref = spmdict
				for f in fields[:-1]:
					if f not in dictref.keys():
						dictref[f] = {}
					dictref = dictref[f]
				dictref[fields[-1]] = self._format_arg(name, spec, value)
			else:
				spmdict[field] = self._format_arg(name, spec, value)
		return [spmdict]
		
    # overwrite matlab command, to save ext
	def _make_matlab_command(self, contents, postscript=None):
		out = os.path.dirname(re.sub("[\[\]']","",str(self.inputs.source)))+"/"
		cmd = super(ROIExtractor,self)._make_matlab_command(contents, postscript)
		cmd = cmd + """
    			if ~exist('ext')
    				ext= evalin('base','ext');
				end
				save('"""+out+"ext.mat','ext');"
		return cmd   
	
	def _list_outputs(self):
		outputs = self._outputs().get()
		out = os.path.dirname(re.sub("[\[\]']","",str(self.inputs.source)))+"/"
		outputs['mat_file'] = out+"ext.mat"
		return outputs
		
###########################################################################		
# Save CSV file in EMBARC format
def save_csv(task,units,names,ext,output="output.csv"):
	# load .mat file
	import scipy.io as sp
	import numpy
	import os
	outfile = os.getcwd()+"/"+output
	values = sp.loadmat(ext,squeeze_me=True)
	raw = values.get('ext')['raw']
	with open(outfile, "w") as f:
		for i in range(0,len(raw)):
			mean = "N/A"
			std = "N/A"
			if raw[i].dtype.names == None:
				mean = numpy.mean(numpy.ma.masked_invalid(raw[i]))
				std  = numpy.std(numpy.ma.masked_invalid(raw[i]))
			elif "mean" in raw[i].dtype.names:
				mean = raw[i]["mean"].tolist()
				std = raw[i]["std"].tolist()
			f.write(task+","+names[i]+","+str(mean)+","+str(std)+","+units+"\n")
	return outfile
	
	
##############################################################################		
class NuisanceInputSpec(BaseInterfaceInputSpec): 
	source = InputMultiPath(File(exists=True),desc='Unsmoothed Functional 4D file',field="source",mandatory=True)
	brain_mask = InputMultiPath(File(exists=True),desc='Brain mask 3D nifti file',field="mask",mandatory=True)
	white_mask = InputMultiPath(File(exists=True),desc='White matter mask 3D nifti file',field="mask",mandatory=True)
	movement = InputMultiPath(File(exists=True,copyfile=True),desc='Realigned Movement .txt file',field="movement",mandatory=True)
	regressors = File(value='regressors.txt',desc='Output Regressors .txt  File',field='output',usedefault=True, genfile=True, hash_files=False)
	time_repetition=traits.Float(1.5, field='TR', usedefault=True, desc='Time Repetition')

class NuisanceOutputSpec(TraitedSpec):
	regressors = File(exists=True, desc='Generated regressors .txt file')
	
"""
	Nuisance Filtering
"""
class Nuisance(BaseInterface):
	"""
	Run nuisance code

	"""
	input_spec = NuisanceInputSpec
	output_spec = NuisanceOutputSpec
	
	def _run_interface(self, runtime):
		from nipype.interfaces.spm.base import scans_for_fname,scans_for_fnames
		from nipype.utils.filemanip import filename_to_list,list_to_filename

		# setup parameters
		d = dict()
		d['source'] = str(self.inputs.source)
		d['brain']  = str(self.inputs.brain_mask)
		d['white']  = str(self.inputs.white_mask)
		d['movement'] = str(self.inputs.movement)
		d['regressors'] = "'"+str(self.inputs.regressors)+"'"
		d['TR'] = str(self.inputs.time_repetition)
		#d['directory'] = os.path.dirname(re.sub("[\[\]']","",d['source']))
		myscript = Template("""
		warning('off','all');
		nuisance($source,$white,$brain,$movement,$TR,$regressors);
       	exit;
		""").substitute(d)
		mlab = MatlabCommand(script=myscript, mfile=True)
		result = mlab.run()
		return result.runtime

	def _list_outputs(self):
		outputs = self._outputs().get()
		#out = os.path.dirname(re.sub("[\[\]']","",str(self.inputs.source)))+"/"
		outputs['regressors'] = os.getcwd()+"/"+self.inputs.regressors
		return outputs	
		
##############################################################################		
class mCompCorInputSpec(BaseInterfaceInputSpec): 
	source = InputMultiPath(File(exists=True),desc='Unsmoothed Functional 4D file',field="source",mandatory=True)
	brain_mask = InputMultiPath(File(exists=True),desc='Brain mask 3D nifti file',field="mask",mandatory=True)
	white_mask = InputMultiPath(File(exists=True),desc='White matter mask 3D nifti file',field="mask",mandatory=True)
	movement = InputMultiPath(File(exists=True,copyfile=True),desc='Realigned Movement .txt file',field="movement",mandatory=True)
	regressors = File(value='regressors.txt',desc='Output Regressors .txt  File',field='output',usedefault=True, genfile=True, hash_files=False)

class mCompCorOutputSpec(TraitedSpec):
	regressors = File(exists=True, desc='Generated regressors .txt file')
	
"""
	mCompCore Filtering
"""
class mCompCor(BaseInterface):
	"""
	Run mCompCore code

	"""
	input_spec = mCompCorInputSpec
	output_spec = mCompCorOutputSpec
	
	def _run_interface(self, runtime):
		from nipype.interfaces.spm.base import scans_for_fname,scans_for_fnames
		from nipype.utils.filemanip import filename_to_list,list_to_filename

		# setup parameters
		d = dict()
		d['source'] = str(self.inputs.source)
		d['brain']  = str(self.inputs.brain_mask)
		d['white']  = str(self.inputs.white_mask)
		d['movement'] = str(self.inputs.movement)
		d['regressors'] = "'"+str(self.inputs.regressors)+"'"
		#d['directory'] = os.path.dirname(re.sub("[\[\]']","",d['source']))
		myscript = Template("""
		warning('off','all');
		mCompCor($source,$white,$brain,$movement,$regressors);
       	exit;
		""").substitute(d)
		mlab = MatlabCommand(script=myscript, mfile=True)
		result = mlab.run()
		return result.runtime

	def _list_outputs(self):
		outputs = self._outputs().get()
		#out = os.path.dirname(re.sub("[\[\]']","",str(self.inputs.source)))+"/"
		outputs['regressors'] = os.getcwd()+"/"+self.inputs.regressors
		return outputs	

		
##############################################################################		
##############################################################################		
class PrintInputSpec(BaseInterfaceInputSpec): 
	in_file = InputMultiPath(File(exists=True,copyfile=False),desc='Input File',mandatory=True)
	out_file = File(value='print.ps',desc='Output Postscript File',usedefault=True, genfile=True)
class PrintOutputSpec(TraitedSpec):
	out_file = File(exists=True, desc='Output PDF file')
	
"""
	Print NIFTI files OR SPM.mat design matrix OR Realign&Unwarp File
"""
class Print(BaseInterface):
	input_spec = PrintInputSpec
	output_spec = PrintOutputSpec
	
	def _run_interface(self, runtime):
		from nipype.interfaces.spm.base import scans_for_fname,scans_for_fnames
		from nipype.utils.filemanip import filename_to_list,list_to_filename

		# setup parameters
		d = dict()
		d['in_file'] = str(self.inputs.in_file).replace("[", "{").replace("]","}")
		d['out_file']  = str(self.inputs.out_file)
		myscript = Template("""
		warning('off','all');
		
		input = $in_file;
		output = '$out_file';
		
		if isstr(input)
			input = {input};	
		end			
		for in = input
			in = in{1};
			[pathstr,name,ext] = fileparts(in);
		
			if strcmp(ext,'.txt')
				plot_realignment_parameters(in,output);
			elseif strcmp(ext,'.mat')
				plot_design_matrix(in,output);
			elseif strcmp(ext,'.nii') || strcmp(ext,'.img')
				nifti2jpeg(in,'-axial -histogram');
				system(['jpeg2ps ' output ' ' pathstr '/*.jpg' ]);
			end
		end
	   	exit;
		""").substitute(d)
		mlab = MatlabCommand(script=myscript, matlab_cmd="matlab -nodesktop -nosplash",mfile=True)
		result = mlab.run()
		return result.runtime

	def _list_outputs(self):
		outputs = self._outputs().get()
		outputs['out_file'] = os.getcwd()+"/"+self.inputs.out_file
		return outputs

############################################################################
class FLTInputSpec(CommandLineInputSpec):
	in_file = traits.String(desc="Input Subject Directory", exists=True, mandatory=True, argstr="%s")
	out_file = traits.String(desc = "CSV file")
class FLTOutputSpec(TraitedSpec):
	out_file = traits.String(desc = "CSV file")

class FLT(CommandLine):
	_cmd = "create_FLT_summary.pl"
	input_spec = FLTInputSpec
	output_spec = FLTOutputSpec
	def _list_outputs(self):
		outputs = self.output_spec().get()
		outputs['out_file'] = self.inputs.out_file
		return outputs


############################################################################
class CorrelateROIsInputSpec(BaseInterfaceInputSpec): 
	in_files =  InputMultiPath(File(exists=True),desc='Lis of 1D text files with timeseries of averages',mandatory=True)
	roi_names = traits.List(traits.String(),desc='List of ROI names',mandatory=True)
	task_name = traits.String("",desc='Task Name')
	out_file = File(value='output.csv',desc='Generated CSV file with Z-Scores',usedefault=True, genfile=True, hash_files=False)

class CorrelateROIsOutputSpec(TraitedSpec):
	out_file = File(exists=True, desc='Generated CSV file with Z-Scores')
	
"""
	Nuisance Filtering
"""
class CorrelateROIs(BaseInterface):
	"""
	Run nuisance code

	"""
	input_spec = CorrelateROIsInputSpec
	output_spec = CorrelateROIsOutputSpec
	
	def _run_interface(self, runtime):
		from nipype.interfaces.spm.base import scans_for_fname,scans_for_fnames
		from nipype.utils.filemanip import filename_to_list,list_to_filename

		# setup parameters
		d = dict()
		d['in_files'] = str(self.inputs.in_files).replace("[", "{").replace("]","}")
		d['roi_names']  = str(self.inputs.roi_names).replace("[", "{").replace("]","}")
		d['task_name']  = str(self.inputs.task_name)
		d['out_file'] = str(self.inputs.out_file)
		
		myscript = Template("""
		warning('off','all');
		
		in_files = $in_files;
		roi_names = $roi_names;
		task_name = '$task_name';
		out_file = '$out_file';
		
		% sanity check
		if length(roi_names) ~= length(in_files)
			fprint('Error: number of 1D average files not the same as their names');
			exit;
		end
		
		% import 1D ROI average files into N vectors
		roi_img = cell(1,length(in_files));
		for i=1:length(in_files)          
			roi_img{i} = importdata(in_files{i});
		end
		
		% calculate Z-scores
		fid = fopen(out_file, 'wt');
		for i=1:length(roi_img)  
			for j=1:length(roi_img)   
				if i ~= j                   
 					z=atanh(corr(roi_img{i},roi_img{j}));
 					fprintf(fid,'%s,%s_%s,%d,N/A,Z_Score\\n',task_name,roi_names{i},roi_names{j},z);
				end
			end
		end
		fclose(fid);

       	
       


		exit;
		""").substitute(d)
		mlab = MatlabCommand(script=myscript, mfile=True)
		result = mlab.run()
		return result.runtime

	def _list_outputs(self):
		outputs = self._outputs().get()
		outputs['out_file'] = os.getcwd()+"/"+self.inputs.out_file
		return outputs	

##############################################################
class ColumnSelectInputSpec(StdOutCommandLineInputSpec):
	in_file = File(desc = "Input Delimited File", exists = True, mandatory = True, argstr="%s")
	delimeter = traits.String(desc="Delimter", mandatory = False, position = 0, argstr="-d %s")
	selection = traits.String(desc = "Select a set of columns",mandatory = True, position = 1, argstr="-f %s" )
	complement = traits.Bool(desc = "Invert column selection", mandatory = False, position = 2,argstr="--complement")	
class ColumnSelectOutputSpec(TraitedSpec):
	out_file = File(exists=True, desc='Output Delimeted file')

class ColumnSelect(StdOutCommandLine):
	input_spec = ColumnSelectInputSpec
	output_spec = ColumnSelectOutputSpec
	cmd = 'cut'
	
	def _gen_outfilename(self):
		return self.inputs.in_file + ".subset"

	def _list_outputs(self):
		outputs = self.output_spec().get()
		outputs['out_file'] = os.path.abspath(self.inputs.in_file + ".subset")
		return outputs
