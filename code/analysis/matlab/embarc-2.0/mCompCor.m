%
% Nuisance inputs: unsmoothed, but processed 4D func, brain_mask, white_test3.nii 
% ROI, rp motion params, fsl motion outliers text files
% output: confounds
%
%
function out = mCompCor(unsmoothed_nii_path,white_mask_roi_path,brain_mask_file_path,rp_movement_file_path,out)

	%dir = '/data/DIAMOND/HC/diamond_20140723.10001/';
	%unsmoothed_nii_path = [dir '/analysis/dynamic_faces/preprocess/despike/rmeas_MID171_Dynamic_Faces_twix_FID29429_MB_flirt_warp_despike.nii'];
	%white_mask_roi_path = '/usr/local/software/data/ROI_EMBARC_scaled2mm/white_test3.nii';
	%brain_mask_file_path = [dir '/output/dynamic_faces/data/mask/rmeas_MID171_Dynamic_Faces_twix_FID29429_MB_flirt_warp_despike_smooth_brain_mask.nii'];
	%rp_movement_file_path = [dir '/output/dynamic_faces/data/movement/rp_meas_MID171_Dynamic_Faces_twix_FID29429_MB.txt'];
	
	if ~exist('out')
		out = 'confounds.txt';
	end
	
	% white matter ROI
	mask_hdr = spm_vol(white_mask_roi_path);
	mask_img = spm_read_vols(mask_hdr);

	% read in movement parameters
	mov_params = importdata(rp_movement_file_path);

	% we need to use data before smoothing in preprocessing, so I'll need to add another output
	% to preprocess() function for unsmoothed data and feed it to this script
	uhdr = spm_vol(unsmoothed_nii_path);
	uimg = spm_read_vols(uhdr);
	NumScans = length(uimg);
	ScanDim  = mask_hdr.dim;
	
	% constants being generated for hi-pass/lo-pass filter
	%f = (0:NumScans-1)*(0.5/NumScans);
	%BPF = ((0.01<abs(f))& (abs(f) < 0.08));
	%BPF2 = ((0.009<abs(f))& (abs(f) < 0.6));

	% Read in brain mask (BET output from preprocess())
	mask3_hdr = spm_vol(brain_mask_file_path);
	mask3_img = spm_read_vols(mask3_hdr);

	% check mask vs white dimensions
	
	n = 0;
	nn = 0;
	nnn = 0;


	% preallocate space
	store = zeros(3,prod(ScanDim));
	list = zeros(prod(ScanDim),4);
	
	%% better initialize?
	list2 = 0;
	list3 = 0;
	lijst = zeros(1000, 3);

	for x = 1:ScanDim(1)
	    for y = 1:ScanDim(2)
	        for z = 1:ScanDim(3)
	            if (mask_img(x, y, z) >= 1)
	                % create a list of coordinates of white matter
	                n= n+1;
	                store(1, n) = x;
	                store(2, n) = y;
	                store(3, n) = z;
	                
	            elseif x<mask3_hdr.dim(1) & y<mask3_hdr.dim(2) & z<mask3_hdr.dim(3) & mask3_img(x,y,z) == 1
	                % create a list of standard deviations per each
	                % voxel and voxel locations
	                nn = nn+1;
	                list(nn, 1) = std(uimg(x, y, z, :));
	                list(nn, 2) = x;
	                list(nn, 3) = y;
	                list(nn, 4) = z;
	            end
	            
% 	            if x<mask3_hdr.dim(1) & y<mask3_hdr.dim(2) & z<mask3_hdr.dim(3) & mask3_img(x,y,z) == 1 
%                 	% global signal 
%                 	nnn = nnn+1;
%                 	lijst(nnn, 1) = x;
%                 	lijst(nnn, 2) = y;
%                 	lijst(nnn, 3) = z;
%             	end
	        end
	    end
	end
	


	%truncate
	store = store(:,1:n);
	list = list(1:nn,:);

	% create an array of timeseries in white matter
	array = zeros(n, NumScans);
	for q = 1:n
		array(q, :) = uimg(store(1, q), store(2, q), store(3, q), :);
	end

% 	% create an array of timeseries in global signal
% 	array2 = zeros(nnn, NumScans);
% 	for q = 1:nnn
% 	        array2(q, :) = uimg(lijst(q, 1), lijst(q, 2), lijst(q, 3), :);
% 	end
% 	global_sig = mean(array2, 1);


	% output of above
	% a spreadsheet with each column being a voxel in ROI
	% and each row is a timeseries


	% resort gray matter matrix by sd
	list2 = sortrows(list, 1);

	% top 2% of voxels with hightest SD
	% in list 3
	len = length(list2);
	lenny = int32(len - (len/50));
	l = len-lenny;
	list3 = zeros(l,3);
	for q = lenny:len
	    list3(q-lenny+1,1) = list2(q,2);
	    list3(q-lenny+1,2) = list2(q,3);
	    list3(q-lenny+1,3) = list2(q,4);
	end

	% timeseries of the top 2% of voxels (new_array)
	new_array = zeros(l, NumScans);
	for q = 1:l
	  new_array(q, :) = uimg(list3(q, 1), list3(q, 2), list3(q, 3), :);
	end

	% do principal component analysis (PCA)

	%%% better initialize?
	third_array = [array' new_array'];
	third_array = third_array';
    mean_array1 = nanmean(third_array, 1);
    
f = (0:NumScans-1)*(0.67/NumScans);

BPF_high = ((0.15<abs(f))& (abs(f) < 1));

spektrum =(fft(mean_array1(:)-(mean(mean_array1(:)))))/NumScans;

spektrum1 = BPF_high'.*spektrum;

mean_array = real(ifft((spektrum1)));


% 	[M, N] = size(third_array);
% 	mn = mean(third_array, 2);
% 	third_array = third_array - repmat(mn, 1, N);
% 	ta = fft(third_array)/NumScans;
% 	bp = repmat(BPF2, M, 1);
% 	ta2 = bp.*ta;
% 	third_array1 = real(ifft(ta2));
% 
% 
% 	Y = third_array1/(sqrt(N-1));
% 	[u, S, PC] = svd(Y);
% 
% 	p_confounds = PC(:, 1:5);
% 	p_confounds(:,6) = global_sig;

% 	m_confounds = mov_params;

	% analyze motion params and do FFT (fourier) 
% 	for k = 1:6
% 	   spektrum =(fft(m_confounds(:,k)'-mean(m_confounds(:,k)')))/NumScans;
% 	   spektrum1 = BPF2.*spektrum;
% 	   m_confounds(1:NumScans, k) = real(ifft((spektrum1)));
% 	end
% 
% 	for c = 1:6
% 	    m_confounds(1, c+6) = 0;
% 	    for t = 2:NumScans 
% 	    m_confounds(t, c+6) = m_confounds(t, c) - m_confounds(t-1, c);
% 	    end
% 	    
% 	    p_confounds(1, c+6) = 0;
% 	    for t = 2:NumScans 
% 	    p_confounds(t, c+6) = p_confounds(t, c) - p_confounds(t-1, c);
% 	    end
% 
% 	end


	confounds = [mean_array, mov_params];


	% now creates an output file
	dlmwrite(out,confounds,'\t');
