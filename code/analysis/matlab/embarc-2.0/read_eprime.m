%
% This is a MATLAB rewrite of a generic E-Prime TXT file
% reader that creates a ePrime stucture as output
% rewrite of Perl eprime2csv file
% author: Eugene Tseytlin (University of Pittsburgh)
%

function  eprime = read_eprime(eprime_file)
% read in the input file
fid = fopen(eprime_file);
tline = fgets(fid);
while ischar(tline)
    tline = strtrim(tline);
    
    % look at cues for start and stop of frame or header
    hd = regexp(tline,'\*\*\* ([A-Za-z]+) (Start|End) \*\*\*','tokens'); %"tokens" breaks the matched values into an array
    %%% (EX 1) *** Logframe Start *** would be tokenized as hd{1}{1} = 'Logframe' and hd{1}{2} = 'Start' %%%
    kv = regexp(tline,'([\w\.]+): (.*)','tokens'); %Find key value pairs of the form "Key: Value"
    
    if ~isempty(hd) && strcmp(hd{1}{2},'Start') % Look for either 'Header Start' or 'Logframe Start'
        % start new block
        block = strtrim(hd{1}{1}); % Look at (EX 1)
        if ~strcmp(block,'Header') %%% *** LogFrame Start ***
            try % try to increment value by 1 
                eprime.(block).(level).N = 1 + eprime.(block).(level).N;
            catch % If value cannot be incremented, set it equal to 1
                eprime.(block).(level).N = 1;
            end
        end
    elseif ~isempty(hd) && strcmp(hd{1}{2},'End') && ~strcmp(block,'Header') % If line is equal to '*** LogFrame End ***'
        % add blank values to key/val that has
        % not been mentioned
        f = fieldnames(eprime.(block).(level)); %fieldnames will give keys that have already been found
        for i=1:length(f) % Iterate through each key 
            %f{i}  %%%TESTING%%%
            %eprime.(block).(level).(f{i})   %%%TESTING%%%
            x = length(eprime.(block).(level).(f{i})); % How many values exists for key=f{i}
            y = eprime.(block).(level).N; % How many values should exist
            if ~strcmp(f{i},'N') && x < y % If less values exist than there should be, then execute this block
                for j=x:y-1 %Iterate through elements starting with most recent and 
			    %finishing on y as specified by number of levels present
                    eprime.(block).(level).(f{i}){end+1} = ''; % append empty values 
                end
            end
        end

    %Block to set (level) for eprime object
    elseif ~isempty(kv) && strcmp(kv{1}{1},'Level')
        % read Level: 1
        level = ['Level' strtrim(kv{1}{2})];

%%% THIS IS WHERE THE ACTUAL WORK IS DONE %%%   
    % We know this block won't have any lines of the following forms ['Header Start|End', 'LogFrame Start|End', 'Level: N']
    elseif ~isempty(kv)
        % while in block, add to current map
        key = strrep(strtrim(kv{1}{1}),'.','_'); %Set key to be the first value 
        val = strtrim(kv{1}{2});
        %%% (Ex 2) 'ENBackList.Cycle: 1' will have the following values: key='ENBackList_Cycle' and value=1

        % Execute this block if we are in the 'Header' block

        if strcmp(block,'Header')
            try
                eprime.(block).(key); % Has this key been set yet?
            catch
                eprime.(block).(key) = {}; % Well set it then!
            end
            eprime.(block).(key){end+1} = val; %append val to the list of values for the key
 
        else % We are in a 'LogFrame' block

            try
                eprime.(block).(level).(key); % What about this key? Has it been set?
            catch
                eprime.(block).(level).(key) = {}; %Initialize value
                if eprime.(block).(level).N > 1 % If this the first iteration of the Level, then we can skip this.
                    for j=1:eprime.(block).(level).N-1 % Otherwise, we must fill in the previous levels for the key with blank values
                        eprime.(block).(level).(key){end+1} = '';
                    end
                end
            end
            eprime.(block).(level).(key){end+1} = val; %append the value to the list of values for this key
        end
    end
    tline = fgets(fid);
end
fclose(fid);
