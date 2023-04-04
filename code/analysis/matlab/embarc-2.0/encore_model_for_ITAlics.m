function encore_model_for_ITAlics(eprime_file)
   %subname = 'sub-21001'
   % Get the basename of eprime_file and set it equal to subname
   [output, subname, ~] = fileparts(eprime_file)
   %Split subname based on "_" and pull second element
   subname = strsplit(subname, '_');
   subname = ['sub-',subname{2}]
   subj{1,1} = ([subname, 'a']);
   subj{1,2} = ([subname, 'b']);
   subj{1,3} = (['BHV_e2_','*.', subname(5:9)]);   

   ss = 1;

   err_n = 0;
   % if ss == 3
   % err_n = 1;    
   % end
   % if ss == 7
   % err_n = 1;    
   % end
    ss
   %other?

   bhv = dir(eprime_file)
    bhv
   load([bhv.folder,'/', bhv.name]);


   av_value = -1; %*(p{1,1}.muPhi(1));
   ap_value = 1;


   
   %clear/initialize new variables
   
   flag = 0;
   %   mkdir(['/projects/ENCORE_2/data/First_levels/first_level_', char(subj{ss,3}(end-4:end)), '/parametric/']); %FIX
   bct = 0;
   choose_safe = 0;
   design_mat = 0;
   for sq = 1:6
      ct=0;
      win_ct = 0;
   nowin_ct = 0;
   neu_ct = 0;
   neg_ct = 0;

   max_av_value=0;
   min_av_value=0;
   max_ap_value=0;
   min_ap_value=0;
   spin_value = 0;
   onset_choice = 0;
   duration_choice =0;
   onset_neg= 0;
   duration_neg =0;
   onset_neu = 0;
   duration_neu =0;
   onset_nowin = 0;
   duration_nowin =0;
   onset_spin = 0;
   duration_spin =0;
   onset_win = 0;
   duration_win =0;

   for iii = 1:20 %ntrials

   ii = (20*(sq-1)) +iii;

   if isnan(results.choiceResp(ii)) == 0

   ct = ct + 1;    
      
   max_av_value(ct) = max(results.negProb1(ii), results.negProb2(ii));
   min_av_value(ct) = min(results.negProb1(ii), results.negProb2(ii));
   max_ap_value(ct) = max(results.posProb1(ii), results.posProb2(ii));
   min_ap_value(ct) = min(results.posProb1(ii), results.posProb2(ii));
   flag(ct) = 0;
   %choose_safe(ct) = 2;

   if results.isFCtrial(ii) == 0
   bct = bct + 1;
   design_mat(bct,1) = max_av_value(ct);
   design_mat(bct,2) = min_av_value(ct);
   design_mat(bct,3) = max_ap_value(ct);
   design_mat(bct,4) = min_ap_value(ct);

   if (results.choiceResp(ii)==1) && (results.negProb1(ii)<results.negProb2(ii))
   choose_safe(bct,1) = 1;
   flag(ct) = 1;
   end

   if (results.choiceResp(ii)==1) && (results.negProb1(ii)>results.negProb2(ii))
   choose_safe(bct,1) = 0;
   flag(ct) = 1;
   end
   if (results.choiceResp(ii) == 2) && (results.negProb1(ii)>results.negProb2(ii))
   choose_safe(bct,1) = 1;
   flag(ct) = 1;
   end

   if (results.choiceResp(ii) == 2) && (results.negProb1(ii)<results.negProb2(ii))
   choose_safe(bct,1) = 0;
   flag(ct) = 1;
   end
   end

   choice_trial(ct) = results.isFCtrial(ii);


   onset_choice(ct) = results.timing.relative.onsetChoicePhase(ii); %results.timing.absolute.onsetFixation(ii) - results.timing.absolute.scanBlockStartTTLglobal(1);
   duration_choice(ct) = results.timing.relative.offsetChoicePhase(ii) - results.timing.relative.onsetChoicePhase(ii);



   onset_spin(ct) = results.timing.relative.onsetSpinning(ii); % - results.timing.absolute.scanBlockStartTTLglobal(1);
   duration_spin(ct) = results.timing.relative.offsetSpinning(ii) - results.timing.relative.onsetSpinning(ii);



   if results.finalSpinValence(ii) == -1
   spin_value(ct) = av_value*results.finalSpinProb(ii);   
      
   if results.outcome(ii) == -1
   neg_ct = neg_ct + 1;
   onset_neg(neg_ct) = results.timing.relative.onsetOutcomeImage(ii); % - results.timing.absolute.scanBlockStartTTLglobal(1);
   duration_neg(neg_ct) = results.timing.relative.offsetOutcomeImage(ii) - results.timing.relative.onsetOutcomeImage(ii);
   end
   if results.outcome(ii) == 0
   neu_ct = neu_ct + 1;
   onset_neu(neu_ct) = results.timing.relative.onsetOutcomeImage(ii); % - results.timing.absolute.scanBlockStartTTLglobal(1);
   duration_neu(neu_ct) = results.timing.relative.offsetOutcomeImage(ii) - results.timing.relative.onsetOutcomeImage(ii);
   end     
   end

   if results.finalSpinValence(ii) == 1
   spin_value(ct) = ap_value*results.finalSpinProb(ii);       
      
   if results.outcome(ii) == 1
   win_ct = win_ct + 1;
   onset_win(win_ct) = results.timing.relative.onsetOutcomeImage(ii); % - results.timing.absolute.scanBlockStartTTLglobal(1);
   duration_win(win_ct) = results.timing.relative.offsetOutcomeImage(ii) - results.timing.relative.onsetOutcomeImage(ii);
   end
   if results.outcome(ii) == 0
   nowin_ct = nowin_ct + 1;
   onset_nowin(nowin_ct) = results.timing.relative.onsetOutcomeImage(ii); % - results.timing.absolute.scanBlockStartTTLglobal(1);
   duration_nowin(nowin_ct) = results.timing.relative.offsetOutcomeImage(ii) - results.timing.relative.onsetOutcomeImage(ii);
   end    
   end




   end

   end        

      
   names=cell(6,1);
   onsets=cell(6,1);
   durations=cell(6,1);

   names{1} ='choice';
   names{2}='anticipation';
   names{3}='win';
   names{4}='nowin';
   names{5}='neutral';
   names{6}='negative';

   onsets{1} = onset_choice;
   onsets{2} = onset_spin;
   onsets{3} = onset_win;
   onsets{4} = onset_nowin;
   onsets{5} = onset_neu;
   onsets{6} = onset_neg;

   durations{1} = duration_choice;
   durations{2}=duration_spin; %[6]; % trial specific
   durations{3}=duration_win;
   durations{4}=duration_nowin;
   durations{5}=duration_neu;
   durations{6}=duration_neg;

   pmod = struct('name', {' '}, 'param', {}, 'poly',{});
   pmod(1).name{1} = 'max_av';
   pmod(1).param{1} = max_av_value;
   pmod(1).poly{1} = 1;

   pmod(1).name{2} = 'min_av';
   pmod(1).param{2} = min_av_value;
   pmod(1).poly{2} = 1;

   pmod(1).name{3} = 'max_ap';
   pmod(1).param{3} = max_ap_value;
   pmod(1).poly{3} = 1;

   pmod(1).name{4} = 'min_ap';
   pmod(1).param{4} = min_ap_value;
   pmod(1).poly{4} = 1;

   pmod(2).name{1} = 'av_ant';
   pmod(2).param{1} = spin_value;
   pmod(2).poly{1} = 1;


   % 
      design_matrix1 = char([output,'/DM_', char(num2str(sq)), '_parametric.mat'])
       save(design_matrix1,'names','onsets','durations', 'pmod');

   clear names onsets duration pmod

   end
   
   YY = categorical(choose_safe);
    
    [B, dev, stats] = mnrfit(design_mat, YY);
   
   behav_output(ss, :) = B;
   behav_mat = char([output,'/Data_', subname, '_behave.mat'])
   save(behav_mat, 'behav_output');
   clear max_ap_value min_ap_value max_av_value min_av_value spin_value
   clear neg_ct neu_ct nowin_ct win_ct
   clear u y X
   clear onset_choice onset_neg onset_neu onset_win onset_nowin onset_spin
   clear onset_win vector duration_choice duration_neg duration_neu
   clear duration_nowin duration_win duration_spin

   end
   
%end
