clear all

directory1 = {'/raid/r4/p3/Phillips/LAMS_reward_gPPI'};

addpath('/raid/r4/p2/Phillips/projects/Henry/PPPI/')
addpath('/raid/r3/p1/Phillips/projects/EMBARC/software/spm8/')

Subjects(1)={'401.20110718'};  
Subjects(2)={'407.20111017'}; 
Subjects(3)={'402.20120409'};  
Subjects(4)={'409.20111108'};  
Subjects(5)={'404.20110825'};  
Subjects(6)={'410.20111111'};
Subjects(7)={'405.20110915'};  
Subjects(8)={'411.20120130'};
Subjects(9)={'406.20110929'};  
Subjects(10)={'412.20120210'};



%User input required (region files)
regionfile={[char(directory1), '/VS_3x3x3.img']};

%User input required (region names)
region={'bilateralVS'};

% or individualised file - create VOI for each subject?



load('ppi_master_template.mat')
P.CompContrasts = 0;
P.VOI=char(regionfile);
P.Region=char(region);
P.Tasks = { '1' 'reward' 'control'};
P.Contrasts.left = 'none';
P.Contrasts.right = 'none';
save([char(directory1), '/', char(region), '.mat'],'P');


for i=1:10
    try
        %User input required directory of SPM.mat files
        SDirectory = {[char(directory1), '/', char(Subjects(i)), '/']};
        cd(char(SDirectory))
        
        % mkdir PPI directory for each subject
        %mkdir('PPI')
        
        %User input required (directory same as line 23 above)
        load([char(directory1), '/', char(region), '.mat']);
        
        P.subject=char(Subjects(i));
        P.directory=char(SDirectory);
        
        %User input required (change analysis to be more specific)
        save([char(Subjects(i)), '_analysis_', char(region), '.mat'],'P');
        PPPI([char(Subjects(i)), '_analysis_', char(region), '.mat']);
    catch
        disp(['Failed: ' Subjects(i)])
    end
end

