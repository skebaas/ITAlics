%
% Read DynFaces generated task file
%
%
function ndm_file = dynfaces_task2dm(task_file)
%task_file = 'subject20140723.10001_2_20140723160824.txt';

% read in the input file
fid = fopen(task_file);
tline = fgets(fid);
categories = {};
while ischar(tline)
	tline = strtrim(tline);
	h1 = regexp(tline,'(\w+)=\s*([^=\s]+)','tokens');
	h2 = regexp(tline,'(\w+)','tokens');
	h3 = regexp(tline,'([^\t]+)\t([^\t]+)','tokens');
	% this is a header line with one attr=value
	if length(h1) == 1
		task.(strtrim(h1{1}{1})) = strtrim(h1{1}{2});
	% button sdefinition	
	elseif length (h1) == 3
		for i=1:length(h1)
			task.buttons.(strtrim(h1{i}{1})) = strtrim(h1{i}{2});
		end
	% category definitions	
	elseif length(h2) == 8 && strcmp(h2{1},'Trial')
		for i=1:length(h2)
			category = char(strtrim(h2{i}));
			categories{end+1} = category;
			task.(category) = {};
		end
	elseif length(h2) == 8
		for i=1:length(h2)
			task.(categories{i}){end+1} = strtrim(h2{i});
		end
	elseif length(h3) == 3
		task.started = strtrim(h1{1}{2});
		task.closed = strtrim(h1{2}{2});
		task.length = strtrim(h1{3}{2});
	end
	tline = fgets(fid);
end
fclose(fid);

% define presets
names = {'Anger';'Fear';'Sad';'Happy';'IDmorph';'Shape'}; %'errors'
durations = cell(length(names),1);
durations(1:length(names)) = {1}; % hard -coded for now to 4s
onsets = cell(length(names),1);

% generate offsets
for i=1:length(task.Face)
	% classify face
	j = 0;
	if regexp(char(task.Face{i}),'.+A[A-Z]$')
		j = 1; % Anger
	elseif regexp(char(task.Face{i}),'.+F[A-Z]$')
		j = 2; % Fear
	elseif regexp(char(task.Face{i}),'.+S[A-Z]$')
		j = 3; % SAD
	elseif regexp(char(task.Face{i}),'.+H[A-Z]$')
		j = 4; % Happy
	elseif regexp(char(task.Face{i}),'.+to[0-9][A-Z]$')
		j = 5; % IDmorhp
	elseif regexp(char(task.Face{i}),'^LL.+')
		j = 6; % Shape
	end
	% check for correct response
	%if ~strcmp(char(task.Crect{i}),char(task.Resp{i}))
	%	j = 7; % error
	%end	
	onsets{j}{end+1} = str2double(task.OnSet{i}) / 1000.0;
end


% save file
subject = char(task.SubID);
[path name ext] = fileparts(task_file);
if ~isempty(path)
	path = [path '/'];
end
ndm_file = [path 'nDM_dynface_' subject '.mat'];
save(ndm_file,'names','onsets','durations');



