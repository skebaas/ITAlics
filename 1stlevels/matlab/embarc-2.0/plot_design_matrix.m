function plot_design_matrix(SPM_file,PS_file)
	load(SPM_file);
	spm_DesRep('DesMtx',SPM.xX,reshape({SPM.xY.VY.fname},size(SPM.xY.VY)),SPM.xsDes);
	spm_print(PS_file);
	fg=spm_figure('FindWin','Graphics');
	if ~isempty(fg)
		close(fg);
	end
return;
