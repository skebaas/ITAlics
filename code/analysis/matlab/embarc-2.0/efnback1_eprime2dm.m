% Generate a SPM design matrix file from E-Prime file
% author: Eugene Tseytlin (University of Pittsburgh

function ndm_file = efnback_eprime2dm(eprime_file)

% process eprime file
ep = read_eprime(eprime_file);

% iterate through log3
names = {'zeroblank';'zeroneutral';'zerofear';'zerohappy';'twoblank';
         'twoneutral';'twofear';'twohappy';'instructions';'error'};
durations = cell(length(names),1);
durations(1:length(names)) = {4}; % hard -coded for now to 4s
onsets = cell(length(names),1);
accuracy = cell(length(names),1); %%%AS/TB%%%
instructions = zeros(1,ep.LogFrame.Level3.N);
%instructions = {0;52;104;156;208;260;312;364};
error = [];
map = {};
ready = intmax;

% go over all conditions 
for i=1:ep.LogFrame.Level3.N
    face = ep.LogFrame.Level3.PicList{i};
    back = ep.LogFrame.Level3.ListType{i};
    code = ep.LogFrame.Level3.LetList{i};
    fprintf('\n')
    
    % construct name
    nm = '';
    if strcmp(back,'0-Back')
        nm = 'zero';
    else
        nm = 'two';
    end
    if strcmp(face,'PicListNeu')
        nm = [nm 'neutral'];
    elseif  strcmp(face,'PicListBlank')
        nm = [nm 'blank'];
    elseif  strcmp(face,'PicListPos')
        nm = [nm 'happy'];
    elseif  strcmp(face,'PicListNeg')
        nm = [nm 'fear'];
    end
    
    % save in map

    map.(code).name = nm;
    map.(code).ready = str2double(ep.LogFrame.Level3.GetReady_OnsetTime{i});
    map.(code).onsets = {};
    ready =  min(ready,map.(code).ready);
    instructions(i) = double(round((map.(code).ready-ready)/1000));
end


% construct onsets
n = ep.LogFrame.Level4.N;
codes = fieldnames(map);
for i=1:n
    code = '';
    % find zero/two + face code that onset in question
    % is related to
    for j=1:length(codes)
        if ~isempty(ep.LogFrame.Level4.(codes{j}){i})
            code = codes{j};
            break;
        end
    end
    if ~isempty(code)
	% maybe need to use the first ready?        
	onset = str2double(ep.LogFrame.Level4.Slide3_OnsetTime{i});
        onset = double(round((onset - ready)/1000));
%%%AS/TB%%% Start 
        acc = str2double(ep.LogFrame.Level4.Slide3_ACC{i})
        rt = str2double(ep.LogFrame.Level4.Slide3_RT{i})
        try
            map.(code).accuracy;
        catch
            map.(code).accuracy = {};
        end
        map.(code).accuracy{end+1} = acc; 

%%%AS/TB%%% End
	if str2double(ep.LogFrame.Level4.Slide3_ACC{i}) == 1
            map.(code).onsets{end+1} = onset;
        else
            error{end+1} = onset;
        end
    end
end

% go in pre-determined order of names
for i=1:8
    % find an appropriate onset data
    for j=1:length(codes)
        if strcmp(names{i},map.(codes{j}).name)
            onsets{i} = cell2mat(map.(codes{j}).onsets);
            accuracy{i} = cell2mat(map.(codes{j}).accuracy); %%%AS/TB%%%
            break;
        end
    end
end
% add instructions and errors   
onsets{9} = instructions;
if isempty(error)
	onsets{10} = [0];
else
	onsets{10} = cell2mat(error);
end

% save file
subject = [ep.Header.Subject{1} '.' ep.Header.Session{1}];
[path name ext] = fileparts(eprime_file);
if ~isempty(path)
	path = [path '/'];
end
ndm_file = [path 'nDM_efnback1_' subject '.mat'];
accuracy_file = [path 'acc_efnback1_' subject '.mat']; %%%AS/TB%%%
save(accuracy_file, 'names', 'accuracy');%%%AS/TB%%%
save(ndm_file,'names','onsets','durations');


% check if low accuracy (<75%)
if length(error) > n*.25
	display(['Warning: very low accuracy, error rate is ' num2str(length(error)*n/100) '%']);	
end

%ListType -> 0-Back vs 2-Back
%PicList  -> PicListBlank (blank), PicListNeu (neutral), PicListNeg (fear),PicListPos (happy)
%LetList  -> ListType+PicList combo for LogFrame 4
%GetReady.OnsetTime ->
%Slide3.ACC -> (1,0) if the anwer was correct, if < 70% is wrong, print warning
%Slide3.OnsetTime - > calulate onset time for LetList: onset = if correct? (Slide3-GetReady)/1000
%Slide3.RT -> reaction time, don't need for now, might need to print warning to exclude studies
% Calculate duration: either hard-code 4s or add CrossTime+Slide3.Duration (round to seconds)
% names vector: zero(blank|neutral|fear|happy), two(blank|neutral|fear|happy), instructions,error
% durations: 4s
% onsets: instructions: [0,52,104,156,208,260,312,364]
% onsets for error are concatonation of onsets if ACC was 0

