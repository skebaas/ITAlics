%% Extract Timeseries Data from a VOI(s)
% Due to a large number of dependencies from PPPI, this CANNOT be run
% independently. Much of this script mirrors the contents of spm_regions.
% It has been included here to automate the process inside of the PPPI
% function. Also, computes the regional mean in addition to the
% eigenvariate for those who want it.
% OUTPUT FILENAME: [P.subject '_session' sess '_' region# '.mat']
% OUTPUT FILE CONTAINS xY (same as xY from spm_regions.m)
%
%   License: This m-file is ditributed under the GNU General Public Licence as published by the
%   Free Software Foundation (either version 2, as given in file
%   spm_LICENCE.man, available in the SPM download) as a derivative work; 
%   however, m-file dependencies - developed separately -- may have their own license. 
%   See specific m-file for the license.
%
%   Version 3 -- Modfied on 09/27/2011 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)
%   Neuropsychology Neuroimaging Laboratory, Univ. of Wisconsin - Madison
%   Neuroscience Training Program and Department of Medicine
%   GRECC, William S. Middleton Memorial Veteren's Hospital, Madison, WI
%   GRECC, Edith Norse Rogers Memorial Veterans Hospital, Bedford, MA
%   Department of Neurology, Massachusetts General Hospital and Harvard
%       Medical School
%
%   Allows the use of subject specific masks.
%

%% Begin Program Here.
errorchk=0;
if ~strcmp(P.analysis,'psy')
	vois=2;
else
	vois=1;
end
clear xY
sess=num2str(z);
for c=1:vois
    if (c==1 && isstruct(P.VOI) && ~P.VOI.exact) || (c==2 && isstruct(P.VOI2) && ~P.VOI2.exact)
        if ~isnumeric(P.contrast)
            xY.Ic=str2double(P.contrast);
        else
            xY.Ic=P.contrast;
        end
        if c==1
            matlabbatch{1}.spm.util.voi.spmmat = {[SPM.swd filesep 'SPM.mat']}; %#ok<*SAGROW>
            matlabbatch{1}.spm.util.voi.adjust = xY.Ic;
            matlabbatch{1}.spm.util.voi.session = z;
            [path file ext]=fileparts(P.VOI.VOI);
            if size(P.VOI.masks,1)~=1
                matlabbatch{1}.spm.util.voi.name = [file '_' P.subject '_Sess' sess];
            else
                matlabbatch{1}.spm.util.voi.name = [file '_' P.subject '_Sess'];
            end
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(P.VOI.VOI);
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            if ~isempty(fileparts(SPM.VM.fname))
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr(SPM.VM.fname);
            else
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr([SPM.swd filesep SPM.VM.fname]);
            end
            matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
            expression = 'i1>0 & i2>0';
            if size(P.VOI.masks,1)~=1
                for m = 1:size(P.VOI.masks,2)
                    matlabbatch{1}.spm.util.voi.roi{m+2}.mask.image = cellstr(P.VOI.masks{z,m});
                    matlabbatch{1}.spm.util.voi.roi{m+2}.mask.threshold = P.VOI.thresh(z,m);
                    addexp = sprintf(' & i%d>0',m+2);
                    expression = [expression addexp]; %#ok<*AGROW>
                end
                matlabbatch{1}.spm.util.voi.expression = expression;
                spm_jobman('run',matlabbatch);
                clear matlabbatch
                if exist(deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess' sess '_' sess '.mat']),'file')==2
                    load(deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess' sess]),'xY');
                else
                    errorchk=-3;
                    return
                end
                xY.mask=spm_vol([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess' sess '.img']);
            else
                for m = 1:size(P.VOI.masks,2)
                    matlabbatch{1}.spm.util.voi.roi{m+2}.mask.image = cellstr(P.VOI.masks{m});
                    matlabbatch{1}.spm.util.voi.roi{m+2}.mask.threshold = P.VOI.thresh(m);
                    addexp = sprintf(' & i%d>0',m+2);
                    expression = [expression addexp];
                end
                matlabbatch{1}.spm.util.voi.expression = expression;
                spm_jobman('run',matlabbatch);
                clear matlabbatch
                if exist(deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess_' sess '.mat']),'file')==2
                    load(deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess_' sess]),'xY');
                else
                    errorchk=-3;
                    return
                end
                xY.mask=spm_vol([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess.img']);
            end
        else
            matlabbatch{1}.spm.util.voi.spmmat = {[SPM.swd filesep 'SPM.mat']};
            matlabbatch{1}.spm.util.voi.adjust = xY.Ic;
            matlabbatch{1}.spm.util.voi.session = z;
            [path file ext]=fileparts(P.VOI2.VOI);
            if size(P.VOI2.masks,1)~=1
                matlabbatch{1}.spm.util.voi.name = [file '_' P.subject '_Sess' sess];
            else
                matlabbatch{1}.spm.util.voi.name = [file '_' P.subject '_Sess'];
            end
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(P.VOI2.VOI);
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            if ~isempty(fileparts(SPM.VM.fname))
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr(SPM.VM.fname);
            else
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr([SPM.swd filesep SPM.VM.fname]);
            end
            matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
            expression = 'i1>0 & i2>0';
            if size(P.VOI2.masks,1)~=1
                for m = 1:size(P.VOI2.masks,2)
                    matlabbatch{1}.spm.util.voi.roi{m+2}.mask.image = cellstr(P.VOI2.masks{z,m});
                    matlabbatch{1}.spm.util.voi.roi{m+2}.mask.threshold = P.VOI2.thresh(z,m);
                    addexp = sprintf(' & i%d>0',m+2);
                    expression = [expression addexp];
                end
                matlabbatch{1}.spm.util.voi.expression = expression;
                spm_jobman('run',matlabbatch);
                clear matlabbatch
                if exist(deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess' sess '_' sess '.mat']),'file')==2
                    load(deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess' sess]),'xY');
                else
                    errorchk=-3;
                    return
                end
                xY.mask=spm_vol([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess' sess '.img']);
            else
                for m = 1:size(P.VOI2.masks,2)
                    matlabbatch{1}.spm.util.voi.roi{m+2}.mask.aimage = cellstr(P.VOI2.masks{m});
                    matlabbatch{1}.spm.util.voi.roi{m+2}.mask.threshold = P.VOI2.thresh(m);
                    addexp = sprintf(' & i%d>0',m+2);
                    expression = [expression addexp];
                end
                matlabbatch{1}.spm.util.voi.expression = expression;
                spm_jobman('run',matlabbatch);
                clear matlabbatch
                if exist(deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess_' sess '.mat']),'file')==2
                    load(deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess_' sess]),'xY');
                else
                    errorchk=-3;
                    return
                end
                xY.mask=spm_vol([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess.img']);
            end
        end
    elseif (c==1 && isstruct(P.VOI)) || (c==2 && isstruct(P.VOI2))
        if c==1
            if size(P.VOI.masks,2)~=1
                zz=z;
            else
                zz=1;
            end
            [path file ext]=fileparts(P.VOI.VOI);
            mapparameters.out=deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess_' sess]);
            mapparameters.mask=P.VOI.VOI;
            mapparameters.STAT='none';
            mapparameters.df1=[];
            mapparameters.thresh=P.VOI.thresh(zz,1);
            mapparameters.sign='pos';
            mapparameters.exact=1;
            mapparameters.cluster=P.VOImin;
            peak_nii(P.VOI.masks{zz,1},mapparameters)
            [xY,errorchk]=set_mask(SPM,[mapparameters.out '_thresh' num2str(mapparameters.thresh) '_extent'  num2str(mapparameters.cluster) '_clusters.nii'],region1,P.maskdir);
        else
            if size(P.VOI2.masks,2)~=1
                zz=z;
            else
                zz=1;
            end
            [path file ext]=fileparts(P.VOI2.VOI);
            mapparameters.out=deblank([SPM.swd filesep 'VOI_' file '_' P.subject '_Sess_' sess]);
            mapparameters.mask=P.VOI2.VOI;
            mapparameters.STAT='none';
            mapparameters.df1=[];
            mapparameters.thresh=P.VOI2.thresh(zz,1);
            mapparameters.sign='pos';
            mapparameters.exact=1;
            mapparameters.cluster=P.VOImin;
            peak_nii(P.VOI2.masks{zz,1},mapparameters)
            [xY,errorchk]=set_mask(SPM,[mapparameters.out '_thresh' num2str(mapparameters.thresh) '_extent'  num2str(mapparameters.cluster) '_clusters.nii'],region1,P.maskdir);
        end
    else
        if c==1
            [xY,errorchk]=set_mask(SPM,P.VOI,region1,P.maskdir);
        else
            [xY,errorchk]=set_mask(SPM,P.VOI2,region2,P.maskdir);
        end
    end
    if errorchk==-1 || errorchk==-3
        return;
    end
    numvox=size(xY.XYZmm,2);
    if isfield(P,'VOImin') && numvox<P.VOImin
        fprintf('\n\nWARNING! VOI has %d voxels! Continuing...\n\n',numvox)
        errorchk=-2;
        return;
    else
        fprintf('\nVOI has %d voxels\n',numvox)
    end
    
    clear xY.y xY.yy xY.u xY.v xY.s xY.X0
    xY.Sess=z;
    if ~isnumeric(P.contrast)
        xY.Ic=str2double(P.contrast);
    else
        xY.Ic=P.contrast;
    end
    
    %-Get voxel numbers that are in VOI
    %======================================================================
    xSPM.hdr=spm_vol(SPM.VM.fname);
    [xSPM.Y, xSPM.XYZmm]=spm_read_vols(xSPM.hdr);
    mXYZ=inv(xY.mask.mat)*[xSPM.XYZmm;ones(1,size(xSPM.XYZmm,2))]; % pixel coordinates
    tmpQ = spm_sample_vol(xY.mask,mXYZ(1,:),mXYZ(2,:),mXYZ(3,:),0);
    tmpQ(~isfinite(tmpQ)) = 0;
    if P.FLmask==1 && P.equalroi==0 
       tmpQ=tmpQ.*reshape(xSPM.Y,1,numel(xSPM.Y));
    end
    Q = find(tmpQ);
    if isempty(find(xSPM.Y(Q)))
        try
            xY.Y=spm_read_vols(xY.mask);
            tmpQv=xSPM.Y.*xY.Y; % fails if image data and mask are different dimensions.
            tmpQ(find(xY.Y>0))=xY.Y(find(xY.Y>0));
            Q = find(tmpQ);
            if isempty(find(xSPM.Y(Q)))
                invokecatchstatement
            end
        catch
            errorchk=-1;
            diary on
            disp('Voxels do not exist; image and VOI do not overlap entirely; VOI must be within the image.')
            disp('OR No voxels in VOI. Program will exit')
            diary off
            return;
        end
    elseif sum(isfinite(xSPM.Y(Q)./xSPM.Y(Q)))~=length(Q) && P.equalroi==1
        errorchk=-1;
        diary on
        disp('VOI is larger than dataset. Program will exit');
        diary off
        return;
    end
    [R,C,S]=ndgrid(1:xSPM.hdr.dim(1),1:xSPM.hdr.dim(2),1:xSPM.hdr.dim(3));
	xY.XYZ = [R(:)';C(:)';S(:)'];
    xY.XYZ=xY.XYZ(:,Q);
	clear R C S
    
    %% spm_regions.m was the source of most of the code below.
    %-Get Data
    %======================================================================
    y=spm_get_data(SPM.xY.VY,xY.XYZ);
    y=spm_filter(SPM.xX.K,SPM.xX.W*y);

    %-Computation
	%=======================================================================
    %-Parameter estimates: beta = xX.pKX*xX.K*y
	%---------------------------------------------------------------
	beta=spm_get_data(SPM.Vbeta,xY.XYZ);

    %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
	%---------------------------------------------------------------
	if xY.Ic~=0
       y=y-spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);
    end

    % confounds
    %-----------------------------------------------------------------------
    xY.X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);

	% extract session-specific rows from data and confounds
	%-----------------------------------------------------------------------
	try
		i     = SPM.Sess(xY.Sess).row;
		y     = y(i,:);
		xY.X0 = xY.X0(i,:);
    catch
	end

	% and add session-specific filter confounds
	%-----------------------------------------------------------------------
	try
		xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).X0];
    catch
	end

	%=======================================================================
	try
		xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).KH]; % Compatibility check
    catch
	end

	%=======================================================================
	% Remove null space of X0
	%-----------------------------------------------------------------------
	xY.X0   = xY.X0(:,~~any(xY.X0));

    % compute regional response in terms of first eigenvariate
	%-----------------------------------------------------------------------
	[m n]   = size(y);
	if m > n
		[v s v] = svd(y'*y);
		s       = diag(s);
		v       = v(:,1);
		u       = y*v/sqrt(s(1));
	else
		[u s u] = svd(y*y');
		s       = diag(s);
		u       = u(:,1);
		v       = y'*u/sqrt(s(1));
	end
	d       = sign(sum(v));
	u       = u*d;
	v       = v*d;
	Y       = u*sqrt(s(1)/n);

	% set in structure
	%-----------------------------------------------------------------------
    xY.y     = y;
    xY.yy    = transpose(mean(transpose(y))); %average (not in spm_regions)
	xY.u    = Y; %eigenvariate
	xY.v    = v;
	xY.s    = s;
    
    if c==1
       save([P.subject '_session' sess '_' region1 '.mat'],'xY', '-v6')
    else
       save([P.subject '_session' sess '_' region2 '.mat'],'xY', '-v6')
    end
    clear xY XYZ hdr img beta y d u v Y s
end
return;




