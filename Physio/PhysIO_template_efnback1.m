%% Example script using PhysIO with Matlab only (no SPM needed)
%  For documentation of the parameters, see also tapas_physio_new (e.g., via edit tapas_physio_new)

%% Create default parameter structure with all fields
physio = tapas_physio_new();

%% Individual Parameter settings. Modify to your need and remove default settings
%physio.save_dir = {'physio_out'};
physio.log_files.vendor = 'Siemens';
physio.log_files.cardiac = {'pulse.puls'};
physio.log_files.respiration = {'resp.resp'};
physio.log_files.scan_timing = {'efnback1.dcm'};
physio.log_files.sampling_interval = [];
physio.log_files.align_scan = 'first';
physio.scan_timing.sqpar.Nslices = 70;
physio.scan_timing.sqpar.TR = 1.5;
physio.scan_timing.sqpar.Ndummies = 0;
physio.scan_timing.sqpar.Nscans = 280;
physio.scan_timing.sqpar.time_slice_to_slice = 0.021428571;
physio.scan_timing.sqpar.onset_slice = 1;
physio.scan_timing.sync.method = 'scan_timing_log';
physio.preproc.cardiac.modality = 'PPU';
physio.preproc.cardiac.filter.include = 1;
physio.preproc.cardiac.filter.type = 'butter';
physio.preproc.cardiac.filter.passband = [0.3 9];
physio.preproc.cardiac.initial_cpulse_select.method = 'auto_matched';
physio.preproc.cardiac.initial_cpulse_select.max_heart_rate_bpm = 120;
physio.preproc.cardiac.initial_cpulse_select.file = 'initial_cpulse_kRpeakfile.mat';
physio.preproc.cardiac.initial_cpulse_select.min = 0.4;
physio.preproc.cardiac.posthoc_cpulse_select.method = 'off';
physio.preproc.cardiac.posthoc_cpulse_select.percentile = 80;
physio.preproc.cardiac.posthoc_cpulse_select.upper_thresh = 60;
physio.preproc.cardiac.posthoc_cpulse_select.lower_thresh = 60;
physio.model.output_multiple_regressors = 'multiple_regressors_efnback1.txt';
physio.model.output_physio = 'physio_efnback1.mat';
physio.model.retroicor.include = 1;
physio.model.retroicor.order.c = 3;
physio.model.retroicor.order.r = 4;
physio.model.retroicor.order.cr = 0;
physio.model.rvt.include = 1;
physio.model.rvt.delays = 0:5:20;
physio.model.hrv.include = 1;
physio.model.hrv.delays = 0:6:24;
physio.verbose.level = 0;
physio.verbose.process_log = cell(0, 1);
physio.verbose.fig_handles = zeros(1, 0);
physio.verbose.fig_output_file = 'output_efnback1.jpg'
physio.verbose.use_tabs = false;
physio.verbose.show_figs = false;
physio.verbose.save_figs = true;
physio.ons_secs.c_scaling = 1;
physio.ons_secs.r_scaling = 1;

%% Run physiological recording preprocessing and noise modeling
physio = tapas_physio_main_create_regressors(physio);
