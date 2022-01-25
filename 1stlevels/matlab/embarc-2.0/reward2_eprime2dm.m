function reward_eprime2dm(eprime_file)
	% extract subject if possible
	directory = fileparts(eprime_file);
	participant = 'subject';

	[st data]= system(['eprime2csv ', eprime_file, ' 1']);
	[head content] = parse_csv(data);
	[st2 data2]= system(['eprime2csv ', eprime_file, ' 2']);
	[head2 content2] = parse_csv(data2);
	[st3 data3]= system(['eprime2csv ', eprime_file, ' 3']);
	[head3 content3] = parse_csv(data3);
	[st4 data4]= system(['eprime2csv ', eprime_file, ' 4']);
	[head4 content4] = parse_csv(data4);

	anticipation_onset=0;
	outcome_trial = 0;
	outcome_onset=0;
	ant_duration = 0;
	
	trial_type = 0;
	range = 0;
	re_array = 0;
	pe_array = 0;
    

	GamStim_RESP_c1 = 2;
	GamStim_RESP_c2 = 3;
	

	resp1 = 0;
	response_onset=0;
	response = zeros(24,1)-1;
	ant_win = 0;
	ant_win_onset=0;
	ant_loss = 0;
	ant_loss_onset=0;
	out_winwin = 0;
	out_winwin_onset=0;
	out_winloss = 0;
	out_winloss_onset=0;
	out_losswin = 0;
	out_losswin_onset=0;
	out_lossloss =0;
	out_lossloss_onset=0;
	
    	ant_mix = 0;
	ant_mix_onset=0;
	ant_neu = 0;
	ant_neu_onset=0;
	out_mixwin = 0;
	out_mixwin_onset=0;
	out_mixloss = 0;
	out_mixloss_onset=0;
	out_neuwin = 0;
	out_neuwin_onset=0;
	out_neuloss =0;
	out_neuloss_onset=0;
    
    
   	baseline = 0;
	baseline_onset=0;
	
    	error_trial = 0;
	error_trial_onset=0;

	%%offset
	Procedure =         find(ismember(head3{1}, 'Procedure')==1);	       %1
	GamStim_OnsetTime = find(ismember(head3{1}, 'GamStim.OnsetTime')==1);  %27 
	GamStim_RESP  =     find(ismember(head3{1}, 'GamStim.RESP')==1);       %28 
	GamStim_RT    =     find(ismember(head3{1}, 'GamStim.RT')==1);         %29
	RewardShuffleImage_OnsetTime = find(ismember(head3{1}, 'RewardShuffleImage.OnsetTime')==1);  %46
	FeedbackR_OnsetTime          = find(ismember(head3{1}, 'FeedbackR.OnsetTime')==1);     %22 %16
	FeedbackRN_OnsetTime         = find(ismember(head3{1}, 'FeedbackRN.OnsetTime')==1);          %19
	LossShuffleImage_OnsetTime   = find(ismember(head3{1}, 'LossShuffleImage.OnsetTime')==1);    %40 
	FeedbackL_OnsetTime          = find(ismember(head3{1}, 'FeedbackL.OnsetTime')==1);           %10
	FeedbackLN_OnsetTime         = find(ismember(head3{1}, 'FeedbackLN.OnsetTime')==1);          %13

	% fix
	MixedShuffleImage_OnsetTime = find(ismember(head3{1}, 'MixedShuffleImage.OnsetTime')==1);  %46
	
	NeutralShuffleImage_OnsetTime   = find(ismember(head3{1}, 'NeutralShuffleImage.OnsetTime')==1);    %40 
	
	%
    GotResp = find(ismember(head3{1}, 'GotRespDuration')==1);
    
    	first_scan = str2num(char(content3{1,1}(GamStim_OnsetTime)));


	for ii=1:48 % each block separately
        
	   trial_scan_time = str2num(char(content3{1,1}(GamStim_OnsetTime)));
	    
	    flag = 0;
	    if isempty(str2num(char(content3{ii}(GamStim_RESP)))) == 1
	        
	        flag = 1;
	        error_trial = error_trial +1;
	        error_trial_onset(error_trial) = str2num(char(content3{ii}(GamStim_OnsetTime)))  - first_scan;

	    %change the error length    
        end
	    
        if (isempty(str2num(char(content3{ii}(GotResp)))) == 1) && (flag == 0)
            flag = 1;
	        error_trial = error_trial +1;
	        error_trial_onset(error_trial) = str2num(char(content3{ii}(GamStim_OnsetTime)))  - first_scan;
            
        end
        
        
        if flag == 0;
        
	    if str2num(char(content3{ii,1}(GamStim_RESP))) == GamStim_RESP_c1
	        resp1 = resp1 + 1;
	        response_onset(resp1) = str2num(char(content3{ii}(GamStim_OnsetTime))) - first_scan;
	        RT(resp1) = str2num(char(content3{ii}(GamStim_RT)));
            
	    end
	    if str2num(char(content3{ii}(GamStim_RESP))) == GamStim_RESP_c2
	        resp1 = resp1 + 1;
	        response_onset(resp1) = str2num(char(content3{ii}(GamStim_OnsetTime))) - first_scan;
	        RT(resp1) = str2num(char(content3{ii}(GamStim_RT)));
            
	    end
	    
	    
	    
	        if (strcmp(content3{ii}(Procedure), 'RewardWinTrialProc') == 1) 
	            
              
	            ant_win = ant_win + 1;
	            out_winwin = out_winwin+1;
	            ant_win_onset(ant_win) = str2num(char(content3{ii}(RewardShuffleImage_OnsetTime))) - first_scan;
		       	out_winwin_onset(out_winwin) = str2num(char(content3{ii}(FeedbackR_OnsetTime ))) - first_scan;
	       		anticipation_onset(resp1) = str2num(char(content3{ii}(RewardShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{ii}(FeedbackR_OnsetTime ))) - first_scan;
                
                ant_duration(resp1) = outcome_onset(resp1) - anticipation_onset(resp1);
                
                trial_type(resp1) = 1;
                
                range(resp1) = 1;
                re_array(resp1) = 0.5;
                pe_array(resp1) = 0.5;
                
			 end
	        
	        if strcmp(content3{ii}(Procedure), 'RewardNeuTrialProc') == 1
	           

	            ant_win = ant_win + 1;
	            out_winloss = out_winloss+1;
	            ant_win_onset(ant_win) = str2num(char(content3{ii}(RewardShuffleImage_OnsetTime))) - first_scan;
	            out_winloss_onset(out_winloss) = str2num(char(content3{ii}(FeedbackRN_OnsetTime  ))) - first_scan;
	            anticipation_onset(resp1) = str2num(char(content3{ii}(RewardShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{ii}(FeedbackRN_OnsetTime ))) - first_scan;
	            
                ant_duration(resp1) = outcome_onset(resp1) - anticipation_onset(resp1);
                
                trial_type(resp1) = 2;
                
                range(resp1) = 1;
                re_array(resp1) = 0.5;
                pe_array(resp1) = -0.5;
                
	        end
	        
	        if strcmp(content3{ii}(Procedure), 'LossNeuTrialProc') == 1
	            
                

	            ant_loss = ant_loss + 1;
	            out_losswin = out_losswin+1;
	            ant_loss_onset(ant_loss) = str2num(char(content3{ii}(LossShuffleImage_OnsetTime))) - first_scan;
	            out_losswin_onset(out_losswin) = str2num(char(content3{ii}(FeedbackLN_OnsetTime)))- first_scan;
	            anticipation_onset(resp1) = str2num(char(content3{ii}(LossShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{ii}(FeedbackLN_OnsetTime ))) - first_scan;
                
                ant_duration(resp1) = outcome_onset(resp1) - anticipation_onset(resp1);
                
                trial_type(resp1) = 3;
                
                range(resp1) = 0.75;
                re_array(resp1) = -0.375;
                pe_array(resp1) = 0.375;
                
	        end
	        
	        if strcmp(content3{ii}(Procedure), 'LossLoseTrialProc') == 1
	                            

	            ant_loss = ant_loss + 1;
	            out_lossloss = out_lossloss+1;
	            ant_loss_onset(ant_loss) = str2num(char(content3{ii}(LossShuffleImage_OnsetTime))) - first_scan;
	            out_lossloss_onset(out_lossloss) = str2num(char(content3{ii}(FeedbackL_OnsetTime)))- first_scan;
	            anticipation_onset(resp1) = str2num(char(content3{ii}(LossShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{ii}(FeedbackL_OnsetTime ))) - first_scan;
                
                ant_duration(resp1) = outcome_onset(resp1) - anticipation_onset(resp1);
                
                trial_type(resp1) = 4;
                
                range(resp1) = 0.75;
                re_array(resp1) = -0.375;
                pe_array(resp1) = -0.375;
                
            end
            
            if strcmp(content3{ii}(Procedure), 'MixedWin') == 1
	            
                
 	            ant_mix = ant_mix + 1;
 	            out_mixwin = out_mixwin+1;
 	            ant_mix_onset(ant_mix) = str2num(char(content3{ii}(MixedShuffleImage_OnsetTime))) - first_scan;
 	            out_mixwin_onset(out_mixwin) = str2num(char(content3{ii}(FeedbackR_OnsetTime)))- first_scan;
 	            anticipation_onset(resp1) = str2num(char(content3{ii}(MixedShuffleImage_OnsetTime))) - first_scan;
	            outcome_onset(resp1) = str2num(char(content3{ii}(FeedbackR_OnsetTime ))) - first_scan;
                
                ant_duration(resp1) = outcome_onset(resp1) - anticipation_onset(resp1);
                
                trial_type(resp1) = 5;
                
                range(resp1) = 1.75;
                re_array(resp1) = 0.125;
                pe_array(resp1) = 0.875;
                
	        end
	        
	        if strcmp(content3{ii}(Procedure), 'MixedLose') == 1
	            
                
 	            ant_mix = ant_mix + 1;
 	            out_mixloss = out_mixloss+1;
 	            ant_mix_onset(ant_mix) = str2num(char(content3{ii}(MixedShuffleImage_OnsetTime))) - first_scan;
 	            out_mixloss_onset(out_mixloss) = str2num(char(content3{ii}(FeedbackL_OnsetTime)))- first_scan;
 	            anticipation_onset(resp1) = str2num(char(content3{ii}(MixedShuffleImage_OnsetTime))) - first_scan;
 	            outcome_onset(resp1) = str2num(char(content3{ii}(FeedbackL_OnsetTime ))) - first_scan;
                
                ant_duration(resp1) = outcome_onset(resp1) - anticipation_onset(resp1);
                
                trial_type(resp1) = 6;
                
                range(resp1) = 1.75;
                re_array(resp1) = 0.125;
                pe_array(resp1) = -0.875;
                
            end
            
            if strcmp(content3{ii}(Procedure), 'NeutralCorrect') == 1
	            
                

 	            ant_neu = ant_neu + 1;
 	            out_neuwin = out_neuwin+1;
 	            ant_neu_onset(ant_neu) = str2num(char(content3{ii}(NeutralShuffleImage_OnsetTime))) - first_scan;
 	            out_neuwin_onset(out_neuwin) = str2num(char(content3{ii}(FeedbackLN_OnsetTime)))- first_scan;
 	            anticipation_onset(resp1) = str2num(char(content3{ii}(NeutralShuffleImage_OnsetTime))) - first_scan;
 	            outcome_onset(resp1) = str2num(char(content3{ii}(FeedbackLN_OnsetTime ))) - first_scan;
                
                ant_duration(resp1) = outcome_onset(resp1) - anticipation_onset(resp1);
                
                trial_type(resp1) = 7;
                
                range(resp1) = 0;
                re_array(resp1) = 0;
                pe_array(resp1) = 0;
	        end
	        
	         if strcmp(content3{ii}(Procedure), 'NeutralIncorrect') == 1
	            
                

 	            ant_neu = ant_neu + 1;
 	            out_neuloss = out_neuloss + 1;
 	            ant_neu_onset(ant_neu) = str2num(char(content3{ii}(NeutralShuffleImage_OnsetTime))) - first_scan;
 	            out_neuloss_onset(out_neuloss) = str2num(char(content3{ii}(FeedbackRN_OnsetTime)))- first_scan;
 	            anticipation_onset(resp1) = str2num(char(content3{ii}(NeutralShuffleImage_OnsetTime))) - first_scan;
 	            outcome_onset(resp1) = str2num(char(content3{ii}(FeedbackRN_OnsetTime ))) - first_scan;
                
		        ant_duration(resp1) = outcome_onset(resp1) - anticipation_onset(resp1);
		        
		        trial_type(resp1) = 8;
		        
		        range(resp1) = 0;
		        re_array(resp1) = 0;
		        pe_array(resp1) = 0;
	        end
            
            
		end
	end

    

	% reward behav
	clear st;
	clear data;
	clear head;
	clear content;

	%if jj == 1
		error_trial1 = error_trial;
	%end
	    

	 if error_trial > 0

	     
	names=cell(4,1);
	onsets=cell(4,1);
	durations=cell(4,1);


	pmod = struct('name', {' '}, 'param', {}, 'poly',{});
	pmod(2).name{1} = 'reward_expectancy';
	pmod(2).param{1} = re_array;
	pmod(2).poly{1} = 1;

	pmod(2).name{2} = 'uncertainty';
	pmod(2).param{2} = range;
	pmod(2).poly{2} = 1;

	pmod(3).name{1} = 'prediction_error';
	pmod(3).param{1} = pe_array;
	pmod(3).poly{1} = 1;

	names{1}='questionmark';
	%names{2}='response'; %?
	names{2}='anticipation';
	names{3}='outcome';
	names{4}='errors';

	onsets{1} = response_onset/1000;
	onsets{2} = anticipation_onset/1000;
	onsets{3} = outcome_onset/1000;
	onsets{4} = error_trial_onset/1000;





	durations{1}=[4];
	durations{2}= ant_duration/1000; %[6]; % trial specific
	durations{3}=[1];
	durations{4}=[5];
	 end

	if error_trial == 0
	 
		names=cell(3,1);
		onsets=cell(3,1);
		durations=cell(3,1);


		pmod = struct('name', {' '}, 'param', {}, 'poly',{});
		pmod(2).name{1} = 'reward_expectancy';
		pmod(2).param{1} = re_array;
		pmod(2).poly{1} = 1;

		pmod(2).name{2} = 'uncertainty';
		pmod(2).param{2} = range;
		pmod(2).poly{2} = 1;

		pmod(3).name{1} = 'prediction_error';
		pmod(3).param{1} = pe_array;
		pmod(3).poly{1} = 1;

		names{1}='questionmark';
		%names{2}='response'; %?
		names{2}='anticipation';
		names{3}='outcome';
		%names{4}='errors';

		onsets{1} = response_onset/1000;
		onsets{2} = anticipation_onset/1000;
		onsets{3} = outcome_onset/1000;
		%onsets{4} = error_trial_onset/1000;


		durations{1}=[4];
		durations{2}= ant_duration/1000; %[6]; % trial specific
		durations{3}=[1];
		%durations{4}=[5];
	 
	end

	      
	% design_matrix = char([directory2, '/nDM_check_reward.mat']);
	% save(design_matrix,'names','onsets','durations', 'pmod');
	 	design_matrix2 = char([directory, '/nDM_RE_PE_reward2.mat']);
	save(design_matrix2,'names','onsets','durations', 'pmod');
	% 
	%error_array(i, jj) = error_trial;
