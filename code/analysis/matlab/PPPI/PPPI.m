function SPM=PPPI(parameterstructure,structfile)
% SPM Toolbox (no GUI) for computing and executing PPIs (and PPPIs) analysis
% at the first-level. PPPI is simply an extension of the GLM that 
% includes psychological interactions and physiological interactions in 
% the same model. Requires SPM8. SPM files must be in the current matlab path.
% Dependencies not explicitly listed below: SPM.mat file that has been estimated;
% VOI.mat file can be generated from a single subject or can be created with
% the eigenvariate gui in SPM' OR a .img/.nii region of interest file may
% also be used. 
%
% WARNING: If you are using SPM2 files, be very careful with your VOI as
% SPM2 were analyze format files and are ussually flipped compared to NIFTI
% files.
%
% Parameter Structure/File Containing Parameter Structure with the following fields:
%      subject: the subject number
%    directory: either the first-level SPM.mat directory, or if you are
%               only estimating a PPI model, then the first-level PPI
%               directory.
%          VOI: name of VOI file ('.nii', '.img', '.mat'). Checks for
%               file before executing program. If you use a .mat file,
%               it should be 3 columns containing the ijk voxel
%               indices OR be a VOI.mat file from SPM.
%               **This can also be a structure. If it is a structure, there
%               are three/four fields:
%               VOI:    name of the VOI file ('.nii', '.img')
%               masks:  a cell array of statistic images ('.nii', '.img') to threshold to
%                       define subject specific ROI. Must be NxM array
%                       (e.g. use {}) where N is either 1 or the number of sessions and M is the number of statistical
%                       images to use to define subject specific ROI.
%               thresh: an NxM matrix of thresholds (e.g. use []) where N is either 1 or the number of sessions and 
%                       M is the number of statistical images to threshold;
%                       thresholds should be the statistic value (e.g. 3)
%                       and not the significance (e.g. .05). These must
%                       line up with the image names in the masks field.
%               exact:  if set to 1, will find a cluster of size VOImin; VOImin must be greater than 0;
%                       default is empty; anything other than 1 will cause this
%                       option to be ignored.
%               **NOTE: Using the exact option will only use the first
%                       image in the masks field. Peak_nii.m must be in the
%                       MATLAB path.
%               **NOTE: each session must use the same number of images to
%                       define the subject specific ROI.
%               **NOTE: if N=1, then the same ROI will be used for each
%                       session. This is recommended.
%       Region: name of output file(s), reqires two names for analysis
%               with two VOI, regions should be separated by a space
%               inside the ' '. Output directory will be Region. (if 2 regions,
%               then the two regions will be separated by a _ in the directory name.
%     contrast: contrast to adjust for. Adjustments remove the effect 
%               of the null space of the contrast. Set to 0 for no adjustment. Set to a
%               number, if you know the contrast number. Set to a contrast name, if you
%               know the name. The default is: 'Omnibus F-test for PPI Analyses'.
%     analysis: specifies psychophysiological interaction ('psy');
%               physiophysiological interaction ('phys'); or psychophysiophysiological 
%               interactions ('psyphy').
%      extract: specifies the method of ROI extraction, eigenvariate ('eig')
%               or mean ('mean')
%       method: specifies traditional SPM PPI ('trad') or generalized
%               condition-specific PPI ('cond')
%     equalroi: specifies the ROIs must be the same size in all subjects
%               NOTE: default=1 (true); set to 0 to lift the restriction
%       FLmask: specifies that the ROI should be restricted using the
%               mask.img from the 1L statistics. NOTE: default=0.
%         VOI2: name of 2nd VOI for physiophysiological interactions ('.nii', '.img', '.mat'). Checks for
%               file before executing program. If you use a .mat file,
%               it should be 3 columns containing the ijk voxel
%               indices OR be a VOI.mat file from SPM.
%               **This can also be a structure. If it is a structure, there
%               are three/four fields:
%               VOI:    name of the VOI file ('.nii', '.img')
%               masks:  a cell array of statistic images ('.nii', '.img') to threshold to
%                       define subject specific ROI. Must be NxM array
%                       (e.g. use {}) where N is either 1 or the number of sessions and M is the number of statistical
%                       images to use to define subject specific ROI.
%               thresh: an NxM matrix of thresholds (e.g. use []) where N is either 1 or the number of sessions and 
%                       M is the number of statistical images to threshold;
%                       thresholds should be the statistic value (e.g. 3)
%                       and not the significance (e.g. .05).These must
%                       line up with the image names in the masks field.
%               exact:  if set to 1, will find a cluster of size VOImin; VOImin must be greater than 0;
%                       default is empty; anything other than 1 will cause this
%                       option to be ignored.
%               **NOTE: Using the exact option will only use the first
%                       image in the masks field. Peak_nii.m must be in the
%                       MATLAB path.
%               **NOTE: each session must use the same number of images to
%                       define the subject specific ROI.
%               **NOTE: if N=1, then the same ROI will be used for each
%                       session. This is recommended.
%      Weights: for traditional PPI, you must specify weight vector for
%               each task. 
%        Tasks: In the generalized condition-specific PPI, you should specify the tasks to
%               include in the analyses, but put a 0 or 1 in front of them to specify if
%               they must exist in all sessions. 
%               For the trad. approach the task must appear in all runs to make the proper
%               contrast weighting, so no number is needed.
%               For the cond. approach the task has to occur in at least 1 run, which is 
%               why you have the option. Default is that it does not have to occur in each run.                
%               Examples:
%                P.Tasks = { '1' 'task1' 'task2' 'task3' 'task4' 'task5' 'task6'} %must exist in all sessions
%                P.Tasks = { '0' 'task1' 'task2' 'task3' 'task4' 'task5' 'task6'} %does not need to exist in all sessions
%               NOTE: In traditional PPI, specify the tasks that go with the weights.
%     Estimate: specifies whether or not to estimate the PPI design. 1 means to 
%               esimate the design, 2 means to estimate the design from already created
%               regressors (must be of the OUT structure), 0 means not to
%               estimate. Default is set to 1, so it will estimate.
%CompContrasts: 0 not to estimate any contrasts;
%               1 to estimate contrasts;  
%               2 to only use PPI txt file for 1st level (not recommended); 
%               3 to only use PPI txt file for 1st level and estimate contrasts (not recommended);
%               2&3 are not recommended as they potentially do not include
%               all tasks effects in the mode. Use at your own risk.
%               3 can not weight the contrasts based on the number of
%               trials
%               Default is 0.
%    Contrasts: cell array of tasks to create contrasts to evaluate OR it is a structure
%               with fields Left and Right that specify the tasks on each
%               side of the equation.
%                 left: tasks on left side of equation or 'none'
%                 right: tasks on right side of equation or 'none'
%                 Weighted: from Weighted above, default is 0.
%                 STAT: 'T' or 'F'
%                 c: contrast vector from createVec, automatically
%                    generated
%                 name: name of contrast, will be defined if left blank
%                 Prefix: prefix to the task name (optional), can be used
%                         to select each run
%                 Contrail: suffix after task name (e.g. parametric
%                           modulators, different basis function)
%                 MinEvents: must be specified and must be 1 or greater,
%                            this tells the program how many events you need to form a
%                            contrast. If there are fewer events, the contrast is not
%                            created.
%                 MinEventsPer: optional. This tells the program how many events per trial type you need to form a
%                               contrast. If there are fewer events, the contrast is not
%                               created. Default is MinEvents/number of
%                               trial type in the contrast.
%               **If left blank and CompContrasts=1, then it defines all
%               possible T contrasts for task components and across runs.
%    Weighted:  Default is not to weight tasks by number of trials (0); to
%               change this, specify which tasks should be weighted by trials.
%               If you want to weight trials, then specify a duration longer
%               than your events. If you have a mixed block event related
%               design, then you can average your events based on number of
%               trials and the blocks won't be averaged IF Weighted is set
%               to be a number that is shorter than the block duration and
%               longer than your events.
%   SPMver:     SPM version used to create SPM.mat files at the first
%               level.
%   maskdir:    location of mask to use (optional)
%   VOImin:     sets the minimum VOI size
%
% FOR DEFAULTS, see PPPI_checkstruct.m code.
%
%EXAMPLE:
%
%  PPPI('parameters.mat'); OR PPPI(parameters,structfile)
%
%Limitations: must compute results for all sessions at the same
%             time
%
%Output text file columns and MAT-files:
%   text files -- OUT.PPI.C OUT.P.C OUT.Y.C OUT.C.C
%   MAT-files  -- OUT structure with fields: PPI (interactions), P
%   (psychological), Y (seed region(s)), and C (covariates). Each of these
%   has a field C containing the values and a field name containing the
%   names of the columns in C.
%
%   PPPI.v9 -- Last modified on 10/03/2011 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)
%   Wisconsin Alzheimer's Disease Research Center - Imaging Core, Univ. of
%   Wisconsin - Madison
%   Neuroscience Training Program and Department of Medicine, Univ. of
%   Wisconsin - Madison
%   GRECC, William S. Middleton Memorial Veteren's Hospital, Madison, WI
%   GRECC, Bedford VAMC
%   Department of Neurology, Massachusetts General Hospital and Havard
%   Medical School
%
%   Much of the section on traditional ppi (psychophysio and physiophysio)
%   code was taken from spm_peb_ppi: 
%      Copyright (C) 2008 Wellcome Department of Imaging Neuroscience
%      Darren Gitelman
%      $Id: spm_peb_ppi.m
%   
%   License: This m-file is ditributed under the GNU General Public Licence as published by the
%   Free Software Foundation (either version 2, as given in file
%   spm_LICENCE.man, available in the SPM download) as a derivative work; 
%   however, m-file dependencies - developed separately -- may have their own license. 
%   See specific m-file for the license.
%
%   Version 4 modifications:
%   (1) Fix for the onset times so the first trial is not accidentally
%   excluded (33 offset was causing a problem with some experimental
%   designs
%   (2) Added slice-timing correction to be consistent with SPM8_r3684
%
%   Version 5 modifications:
%   (1) Changed Task specification from column number to column name
%   (2) Option to allow some subjects/session to have event types missing
%   OR force all subjects to have all events
%   (3) Will estimate the PPI design in addition to making the regressors
%   (4) Motion parameters will be used only if specified in 1st-level
%   model, parameters will be grabbed from SPM.mat file.
%   (5) Options related to 1st level contrast computations was also added
%
%   Version 6 modifications:
%   (1) Cleaned up the code.
%
%   Version 7 modifications:
%   (1) Added the ability to use gzip and bzip2 files (only for Y though)
%   (2) Fixed problem with parametric modulation
%
%   Version 8 modifications:
%   (1) Allows for subject specific masks based on subject statistical
%   images. Thanks to Bob Spunt for providing the code to make this possible.
%
%  Version 9 modifications:
%  (1) Fixes error messages if input options are not correct.
%  (2) Allows the user to specify an exact VOI size (P.VOI.extact/P.VOI2.exact and P.VOImin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
toolboxdir=fileparts(which('PPPI'));
addpath(toolboxdir);%%
ind= find(toolboxdir == filesep);
if max(ind)==numel(toolboxdir)
    addpath(toolboxdir(1:ind(end-1)-1));
else
    addpath(toolboxdir(1:ind(end)-1));
end
if ~exist('createVec','file')
    errorval='createVec is not in path. createVec should be one in the same folder as the PPPI directory.'
end
if ~exist('defContrasts','file')
    errorval='defContrasts is not in path. defContrasts should be one in the same folder as the PPPI directory.'
end

%% Input Parser -- Will create variables for the specified options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(spm('Ver'),'SPM12')
    disp('PROGRAM ABORTED:')
    disp('  You must use SPM8 to process your data; however, you can use SPM.mat files')
    disp('  generated with SPM2 or SPM5. In these cases, simply specify the option SPMver')
    disp('  in single qoutes followed by a comma and the version number.')
    disp(' ')
    disp('Make sure to add SPM8 to your MATLAB path before re-running.')
    return
end
if exist('parameterstructure','var') && (~isstruct(parameterstructure) && exist(parameterstructure,'file')) 
    [path structfile ext]=fileparts(parameterstructure);
    structfile=[structfile ext];
    parameterstructure=load(parameterstructure);
    if ~isfield(parameterstructure,'correct') || parameterstructure.correct~=1
        try
            [P, errorval]=PPPI_checkstruct(parameterstructure,structfile);
        catch
            try 
                [P, errorval]=PPPI_checkstruct(parameterstructure);
            catch
                disp('Program will exit. Input structure is not correct.')
                return;
            end
        end
    else
        P=parameterstructure;
        errorval=[];
    end
elseif exist('parameterstructure','var') && isstruct(parameterstructure)
    if ~isfield(parameterstructure,'correct') || parameterstructure.correct~=1
        try
            [P, errorval]=PPPI_checkstruct(parameterstructure,structfile);
        catch
            try 
                [P, errorval]=PPPI_checkstruct(parameterstructure);
            catch
                disp('Program will exit. Input structure is not correct.')
                return;
            end
        end
    else
        P=parameterstructure;
        errorval=[];
    end
else
    disp('Program will exit. Required input is not a structure or a file.')
    help PPPI
    return;
end

if ~isempty(errorval)
    disp('One or more inputs are not correct.')
    for ee=1:numel(errorval)
        disp(['ERROR ' num2str(ee) ':' errorval{ee}])
    end
    disp('errorvals saved to errorvals.mat')
    save errorvals errorval
    return;
else
    cd(P.directory)
    try
       [path structfile ext]=fileparts(structfile);
       if isempty(path)
           path=pwd;
       end
       save([path filesep P.subject '_' strtok(P.Region) '_' P.analysis P.method '_' structfile '_' date ext],'P');
    catch
       save([pwd filesep P.subject '_' strtok(P.Region) '_' P.analysis P.method '_PPI_structure_ ' date '.mat'],'P'); 
    end
end

%% Check for zipped files.
try
    load SPM.mat
    for ii=1:numel(SPM.xY.VY)
        a{ii}=SPM.xY.VY(ii).fname;
    end
    a=unique(a);
    filesgz={}; filesbz={};
    for ii=1:numel(a)
        if exist([a{ii} '.gz'],'file')==2
            gunzip([a{ii} '.gz'])
            filesgz{end+1}=a{ii};
        end
        if ~exist(a{ii},'file')==2
            error('SPM.xY.VY points to the wrong locations')
        end
    end
    for ii=1:numel(a)
        if exist([a{ii} '.bz2'],'file')==2
            try
                system(sprintf('bunzip2 %s',[a{ii} '.bz2']));
                filesbz{end+1}=a{ii};
            catch
                error('bunzip2 is a valid command. This was setup for *NIX.')
            end
        end
        if ~exist(a{ii},'file')==2
            error('SPM.xY.VY points to the wrong locations')
        end
    end
catch
   error('Checking for zipped files failed') 
end

%% Check for estimate only

if P.Estimate==2
       SPM=spm_estimate_PPI(P.subject,'SPM.mat',P.Region,P.method,P.analysis,P.CompContrasts);
       try 
           if P.CompContrasts==1
               SPM=spm_contrasts_PPI(SPM,P.subject,P.Contrasts,P.Weighted,P.method);
           end
       catch
            disp('PPI Contrasts were not estimated for some reason.')
       end
       zipfiles(filesgz,filesbz)
       return;
end

%% Program
[region1,region2]=strtok(P.Region);
region2=region2(2:end);
load SPM.mat
Sessions=numel(SPM.Sess);

%Begin Session Loop
%==========================================================================
for z=1:Sessions
    timeseries_extract;
    if errorchk==0
    elseif errorchk==-2
        disp('Cluster was smaller than P.VOImin')
        disp('Now exiting.')
        return;
    elseif errorchk==-3
        disp('Statistical Images (P.VOI/VOI2.masks Used in Masking do not intersect P.VOI.VOI and/or P.VOI2.VOI')
        disp('Now exiting.')
        return;
    else
        return;
    end

    % Computation of PPI vectors
	%======================================================================
	RT=SPM.xY.RT;
	dt=SPM.xBF.dt;
	NT=round(RT/dt);
    fMRI_T0 = SPM.xBF.T0; % Corrects for slice-timing offset
    sess=num2str(z);

    if strcmp(P.analysis,'psy')
    % Psychophysiological Interactions
        %==================================================================
	    timeseries =load(strcat(P.subject,'_session',sess,'_',region1,'.mat'), 'xY');
        xY = timeseries.xY;
        Sess = SPM.Sess(z); 

    else
	% Physiophysiological or Psychophysiophysiological interactions
        %==================================================================
		timeseries = load(strcat(P.subject,'_session',sess,'_',region1,'.mat'),'xY');
        xY(1) = timeseries.xY;
        clear timeseries
        timeseries = load(strcat(P.subject,'_session',sess,'_',region2,'.mat'),'xY');
        xY(2) = timeseries.xY;numel(P.Tasks)
        Sess  = SPM.Sess(z);  
    end

  % Setup variables
  %------------------------------------------------------------------------
  N     = length(xY(1).u);
  k     = 1:NT:N*NT;  			% microtime to scan time indices   
  
    if strcmp(P.analysis,'psy') || strcmp(P.analysis,'psyphy')
      u       = length(Sess.U);
      Sess.U=spm_get_ons(SPM,z); %Get onsets in microtime (SPM8 was updated)
	  U.name  ={};
	  U.u     =[];
	  U.w     =[];
          
      if strcmp(P.method,'trad')
        % Traditional Method  
        %==============================================================
        T=zeros(numel(P.Tasks),1);
        I=zeros(numel(P.Tasks),1);
        Match=[];
        for j=1:numel(P.Tasks);
                for jj=1:numel(Sess.U)
                    match=0;
                    if strcmp(P.Tasks{j},Sess.U(jj).name(1))
                        T(j)=jj;
                        I(j)=j; %index of useable P.Tasks
                        match=1;
                        Match(end+1)=1;
                        break
                    end
                end
            if ~match
                    disp([P.Tasks{j} ' does not exist. Program will exit'])
                    Match(end+1)=0;
                    return
            end
        end
        if ~sum(Match)
            disp('Tasks do not exist OR were poorly defined. Program will exit')
            return
        end
        I=I(I~=0); 
        T=T(T~=0);
        for i=1:numel(T)
            for j=1:length(Sess.U(T(i)).name)
                if any(Sess.U(T(i)).u(33:end,j))
                    U.u             = [U.u Sess.U(T(i)).u(33:end,j)];
                    U.name{end + 1} = Sess.U(T(i)).name{j};
                    U.w             = [U.w P.Weights(i)];
                end
            end
        end
        PSY     = zeros(N*NT,1);
        for i = 1:size(U.u,2)
            PSY = PSY + full(U.u(:,i)*U.w(:,i));
            %         PSY = spm_detrend(PSY);  <- removed centering of psych
            %         variable
            %         prior to multiplication with xn. Based on discussion
            %         with Karl
            %         and Donald McLaren.
        end
      else           
        % Condition-Specific Method
        %=================================
        T=zeros(numel(P.Tasks)-1,1);
        I=zeros(numel(P.Tasks)-1,1);
        Match=[];
        for j=2:numel(P.Tasks)
            for jj=1:numel(Sess.U)
                    match=0;
                    if strcmp(P.Tasks{j},Sess.U(jj).name(1))
                        T(j)=jj;
                        I(j)=j; %index of useable P.Tasks
                        match=1;
                        Match(end+1)=1;
                        break
                    end
             end
             if ~match && str2double(P.Tasks{1})~=0
                    disp([P.Tasks{j} ' does not exist. Program will exit'])
                    Match(end+1)=0;
                    return
             end
        end
        if ~sum(Match)
            disp('Tasks do not exist OR were poorly defined. Program will exit')
            return
        end
        I=I(I~=0); 
        T=T(T~=0);
        for i=1:numel(T)
            for j=1:length(Sess.U(T(i)).name)
                if any(Sess.U(T(i)).u(33:end,j))
                    U.u             = [U.u Sess.U(T(i)).u(33:end,j)];
                    U.name{end + 1} = Sess.U(T(i)).name{j};
                    U.w             = 1;
                end
            end
        end 
        for i = 1:size(U.u,2)
                PSY(:,i)     = zeros(N*NT,1);
                PSY(:,i)	 = PSY(:,i) + full(U.u(:,i)*U.w);
                %         PSY = spm_detrend(PSY);  <- removed centering of psych
                %         variable
                %         prior to multiplication with xn. Based on discussion
                %         with Karl
                %         and Donald McLaren.
        end
      end
    end
    
    % create basis functions and hrf in scan time and microtime
	%----------------------------------------------------------------------
	hrf   = spm_hrf(dt);

	% create convolved explanatory {Hxb} variables in scan time
	%-------------------------------------------------------------------------
	xb    = spm_dctmtx(N*NT + 128,N);
	Hxb   = zeros(N,N);
        for i = 1:N
            Hx       = conv(xb(:,i),hrf);
            Hxb(:,i) = Hx(k + 128);
        end
	xb    = xb(129:end,:);

	% get confounds (in scan time) and constant term
	%----------------------------------------------------------------------
	X0    = xY(1).X0;
	M     = size(X0,2);

	% get response variable (timeseries),
	%-------------------------------------------------------------------------
	l=xY.u;
    S=size(xY,2);
    Y=zeros(length(l),S);
    if strcmp(P.extract,'eig')
        for i = 1:S
            Y(:,i) = xY(i).u;
        end
    else
        for i = 1:S
            Y(:,i) = xY(i).yy;
        end
    end

	% remove confounds and save Y in ouput structure
	%----------------------------------------------------------------------
	Yc    = Y - X0*inv(X0'*X0)*X0'*Y;
	PPI.Y = Yc;

	% specify covariance components; assume neuronal response is white
	% treating confounds as fixed effects
	%----------------------------------------------------------------------
	Q      = speye(N,N)*N/trace(Hxb'*Hxb);
	Q      = blkdiag(Q, speye(M,M)*1e6  );

	% get whitening matrix (NB: confounds have already been whitened)
	%----------------------------------------------------------------------
	W      = SPM.xX.W(Sess.row,Sess.row);

	% create structure for spm_PEB0 not to estimate any contrasts;
	%----------------------------------------------------------------------
    PEBP=cell(2,1);
    PEBP{1}.X = [W*Hxb X0];		% Design matrix for lowest level
	PEBP{1}.C = speye(N,N)/4;		% i.i.d assumptions
	PEBP{2}.X = sparse(N + M,1);	% Design matrix for parameters (0's)
	PEBP{2}.C = Q;

    switch P.analysis
        case{'psy'}
        % Psychophysiological interactions
        %==================================================================
    	C       = spm_PEB(Y,PEBP);
    	xn      = xb*C{2}.E(1:N);
        xn      = spm_detrend(xn);
        
        PSYxn=zeros(size(PSY,1),size(PSY,2));
        PSYHRF=zeros(numel((k-1) + fMRI_T0),size(PSY,2));
        for j=1:size(PSY,2)
                % multiply psychological variable by neural signal
                %----------------------------------------------------------
                PSYxn(:,j)   = PSY(:,j).*xn;
                % convolve and resample at each scan for bold signal
                %----------------------------------------------------------
                ppit	    = conv(PSYxn(:,j),hrf);
                ppit        = ppit((k-1) + fMRI_T0);
                ppi(:,j)    = spm_detrend(ppit);
                
                % similarly for psychological effect
                %----------------------------------------------------------
                PSYHRFtmp   = conv(PSY(:,j),hrf);
                PSYHRF(:,j) = PSYHRFtmp((k-1) + fMRI_T0);
        end

        % save psychological variables
    	%------------------------------------------------------------------
    	PPI.P   = PSYHRF(:,any(PSYHRF));
    	PPI.ppi = ppi(:,any(ppi));

        case{'phy'}
        % Physiophysiological Interactions
		%==================================================================
		C       = spm_PEB(Y(:,1),PEBP);
    	xn1     = xb*C{2}.E(1:N);
    	C       = spm_PEB(Y(:,2),PEBP);
    	xn2     = xb*C{2}.E(1:N);
    	xn1     = spm_detrend(xn1);
    	xn2     = spm_detrend(xn2);
    	xnxn    = xn1.*xn2;

    
    	% convolve and resample at each scan for bold signal
    	%------------------------------------------------------------------
    	ppi     = conv(xnxn,hrf);
    	ppi     = ppi((k-1) + fMRI_T0);

    	% save variables
    	%------------------------------------------------------------------
    	PPI.xn  = [xn1 xn2];
    	PPI.ppi = spm_detrend(ppi);
                
        case{'psyphy'}
        % Psychophysiophysiological Interactions
        %==================================================================
        C       = spm_PEB(Y(:,1),PEBP);
        xn1     = xb*C{2}.E(1:N);
    	C       = spm_PEB(Y(:,2),PEBP);
    	xn2     = xb*C{2}.E(1:N);
    	xn1     = spm_detrend(xn1);
    	xn2     = spm_detrend(xn2);
    	xnxn    = xn1.*xn2;
    
    	% convolve and resample at each scan for bold signal
    	%------------------------------------------------------------------
        ppi         = zeros(length(k),3*length(Sess.U)+1); %Initialize ppi
        ppitmp      = conv(xnxn,hrf);
        ppi(:,1)    = spm_detrend(ppitmp((k-1) + fMRI_T0));
 
        PSYxn1=zeros(size(PSY,1),size(PSY,2));
        PSYxn2=zeros(size(PSY,1),size(PSY,2));
        PSYxnxn=zeros(size(PSY,1),size(PSY,2));
        PSYHRF=zeros(numel((k-1) + fMRI_T0),size(PSY,2));
        for j=1:size(PSY,2)
            % multiply psychological variable by neural signal
    		%----------------------------------------------------------
            PSYxn1(:,j)   = PSY(:,j).*xn1;
            PSYxn2(:,j)   = PSY(:,j).*xn2;
            PSYxnxn(:,j)  = PSY(:,j).*xnxn;

            % convolve and resample at each scan for bold signal
            %----------------------------------------------------------
            ppi1                    =   conv(PSYxn1(:,j),hrf);
            ppi(:,(j-1)*3+2)              = spm_detrend(ppi1((k-1) + fMRI_T0));
            ppi2                    = conv(PSYxn2(:,j),hrf);
            ppi(:,(j-1)*3+3)    = spm_detrend(ppi2((k-1) + fMRI_T0));
            ppixnxn                 = conv(PSYxnxn(:,j),hrf);
            ppi(:,(j-1)*3+4)  = spm_detrend(ppixnxn((k-1) + fMRI_T0));

            % similarly for psychological effect
            %--------------------------------------------------------------
            PSYHRFtmp  = conv(PSY(:,j),hrf);
            PSYHRF(:,j)  = PSYHRFtmp((k-1) + fMRI_T0);
        end
    
     
        % save variables
        %---------------------------------------------------------------------
        PPI.xn  = [xn1 xn2];
    	PPI.P   = PSYHRF(:,any(PSYHRF));
    	PPI.ppi = ppi(:,any(ppi)); % xnxninteraction+ntasks*3interactions
    end
       
    % Define Output
    %=====================================================================
     if strcmp(P.analysis,'phy')
        OUT.P.C = [];
        OUT.P.name={};
        OUT.PPI.C=PPI.ppi;
        OUT.Y.C=PPI.Y;
        OUT.Y.name={[region1 '_seedtc'] [region2 '_seedtc']};
        OUT.C.C=SPM.Sess(z).C.C;
        OUT.C.name=SPM.Sess(z).C.name;
    elseif strcmp(P.analysis,'psy')
        OUT.P.C = PPI.P;
        OUT.PPI.C=PPI.ppi;
        if str2double(P.Tasks{1})==0 || str2double(P.Tasks{1})==1
            OUT.P.name={};
            OUT.PPI.name={};
            for i=1:numel(T)
                for nn=1:length(Sess.U(T(i)).name)
                    if any(Sess.U(T(i)).u(33:end,nn))
                        OUT.P.name{end+1}=Sess.U(T(i)).name(nn);
                        OUT.PPI.name{end+1}=['PPI_' cell2mat(Sess.U(T(i)).name(nn))];
                    end
                end
            end
        else
            OUT.P.name={'PSY'};
            OUT.PPI.name={'PPI'};
        end
        OUT.Y.C=PPI.Y;
        OUT.Y.name={[region1 '_seedtc']};
        OUT.C.C=SPM.Sess(z).C.C;
        OUT.C.name=SPM.Sess(z).C.name;
    else
        OUT.P.C = PPI.P;
        OUT.PPI.C=PPI.ppi;
        if str2double(P.Tasks{1})==0 || str2double(P.Tasks{1})==1
            OUT.PPI.name={['PPI_' region1 '_' region2]};
            OUT.P.name={};
            for i=1:numel(T)
                for nn=1:numel(Sess.U(T(i)).name)
                    if any(Sess.U(T(i)).u(33:end,nn))
                        OUT.P.name{end+1}=Sess.U(T(i)).name(nn);
                    end
                end
                suffix={region1 region2 [region1 '_' region2]};
                for nn=1:numel(Sess.U(T(i)).name)
                    if any(Sess.U(T(i)).u(33:end,nn))
                        for jj=1:numel(suffix)
                            OUT.PPI.name{end+1}=['PPI_' cell2mat(Sess.U(T(i)).name(nn)) '_' cell2mat(suffix(jj))];
                        end
                    end
                end
            end
        else
            OUT.P.name={'PSY'};
            OUT.PPI.name={['PPI_' region1 '_' region2] ['PPI_' region1] ['PPI_PSY_' region2] ['PPI_PSY_' region1 '_' region2]};
        end
        OUT.Y.C=PPI.Y;
        OUT.Y.name={[region1 '_seedtc'] [region2 '_seedtc']};
        OUT.C.C=SPM.Sess(z).C.C;
        OUT.C.name=[SPM.Sess(z).C.name];  
    end  

    % Save Output
    %======================================================================
    if strcmp(P.analysis,'psy')
        dlmwrite([P.subject '_' region1 '_session' num2str(z) '_' P.method '_PPI_regressors.txt'],[OUT.PPI.C OUT.P.C OUT.Y.C OUT.C.C],' ');
        save([P.subject '_' region1 '_session' num2str(z) '_' P.method '_PPI_regressors.mat'],'OUT','-v7.3');
    elseif strcmp(P.analysis,'phy')
        dlmwrite([P.subject '_' region1 '_and_' region2 '_session' num2str(z) '_PPI_regressors.txt'],[OUT.PPI.C OUT.P.C OUT.Y.C OUT.C.C],' ');
        save([P.subject '_' region1 '_and_' region2 '_session' num2str(z) '_PPI_regressors.mat'],'OUT','-v7.3');
    else
        dlmwrite([P.subject '_' region1 '_and_' region2 '_session' num2str(z) '_' P.method '_PPPI_regressors.txt'],[OUT.PPI.C OUT.P.C OUT.Y.C OUT.C.C],' ')
        save([P.subject '_' region1 '_and_' region2 '_session' num2str(z) '_' P.method '_PPPI_regressors.mat'],'OUT','-v7.3')
    end
    clear PPI PSY ppi OUT
end

if P.Estimate==1
   try
        SPM=spm_estimate_PPI(P.subject,'SPM.mat',P.Region,P.method,P.analysis,P.CompContrasts);
   catch
       disp('Estimation Failed')
   end
   try 
      if P.CompContrasts==1
         SPM=spm_contrasts_PPI(SPM,P.subject,P.Contrasts,P.Weighted,P.method);
      end
   catch
      disp('PPI Contrasts were not estimated for some reason.')
   end
else
   SPM={'PPI was not estimated. Please estimate the PPI model now.'};
end
zipfiles(filesgz,filesbz)
end
%% Check for zipped files.
function zipfiles(filesgz,filesbz)
try
    if ~isempty(filesgz)
        for ii=1:numel(filesgz)
            gzip(files{ii});
        end
    end
    if ~isempty(filesbz)
        for ii=1:numel(filesbz)
            system(sprintf('bzip2 %s',a{ii}));
        end
    end
catch
   error('Zipping files failed') 
end
end
