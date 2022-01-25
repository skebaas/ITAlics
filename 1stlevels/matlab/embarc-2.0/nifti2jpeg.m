% nifti2jpeg - create JPEG files for NIFTI file
% author: Eugene Tseytlin (University of Pittsburgh)
% file : filename of nii image
% flags: -axial  -sagital -coronal -3d -histogram 
function nifti2jpeg(file, flags)
	warning('off','all');
    % setup flags
	if ~exist('flags') || isempty(flags)
		flags = '-axial -sagital -coronal -histogram';
    %elseif ~isempty(findstr(flags,'-3d'))
	%	flags = [flags ' -axial -sagital -coronal']; 
	end

	%
	% strip .gz suffix
	%
	if strcmp(file(length(file)-2:length(file)),'.gz') 
		dos(['gunzip ' file]);
		file = file(1:length(file)-3);
	end
	
	% load nifti file
    out=file(1:length(file)-4);
	nii=load_untouch_nii(file);
	
	%axial: 
	if ~isempty(findstr(flags,'-axial'))
		figure
		montage(imrotate(permute(nii.img(:,:,:,1),[1 2 4 3]), 90),'DisplayRange',[])
		colormap gray
		colorbar
		title(['axial of ' file]);
        saveas(gcf,[out '_axial.jpg']);
	end
	
	%sagital: 
	if ~isempty(findstr(flags,'-sagital'))
		figure
		montage(imrotate(permute(nii.img(:,:,:,1),[2 3 4 1]), 90),'DisplayRange',[])
		colormap gray
		colorbar
        title(['sagital of ' file]);
		saveas (gcf, [out '_sagital.jpg']);
	end
	
	%coronal: 
	if ~isempty(findstr(flags,'-coronal'))
		figure
		montage(imrotate(permute(nii.img(:,:,:,1),[1 3 4 2]), 90),'DisplayRange',[])
		colormap gray
		colorbar
        title(['coronal of ' file]);
		saveas (gcf,[ out '_coronal.jpg']);
	end
	
	%3d image
	%if ~isempty(findstr(flags,'-3d'))
	%	[axial, map1]=imread([out '_axial.jpg']);
	%	[sagital, map2]=imread([out '_sagital.jpg']);
	%	[coronal, map3]=imread([out '_coronal.jpg']);
	%	image1=cat(1, axial, sagital, coronal);
	%	figure
	%	montage(image1, 'DisplayRange',[]);
    %   title(['3D of ' file]);
	%	imwrite(image1, [out '_3D.jpg']);
	%end
	% create histogram image
	if ~isempty(findstr(flags,'-histogram'))
		hdr = spm_vol(file);
		[n, x] = histvol(hdr, 100);
		figure;
		bar(x,n);
		title(['histogram of ' file]);
		saveas (gcf,[out '_histogram.jpg']);
    end
    close all;
    clear nii;
 
