clear all

directory1 = {'/raid/r4/p2/Phillips/projects/Matt/Matt_Design_Files/MDD_vs_HC_Project/FIrstLevels_and_PPI/PreProcessing_2/FL_PreProcessing_2_Despike_MDx1_NewGN'};

addpath('/raid/r4/p2/Phillips/projects/Henry/PPPI/')
addpath('/raid/r3/p1/Phillips/projects/EMBARC/software/spm8/')


Subjects(1)={'fl_V2_HC_20080407_210315p_c4'};
Subjects(2)={'fl_V2_MDD_20080918_210825p_c3'};
Subjects(3)={'fl_V2_HC_20080501_210335_c2'};   
Subjects(4)={'fl_V2_MDD_20081008_203141p_c4'};
Subjects(5)={'fl_V2_HC_20080610_210460_c2'};   
Subjects(6)={'fl_V2_MDD_20081028_210911_c3'};
Subjects(7)={'fl_V2_HC_20080617_210446_c3'};    
Subjects(8)={'fl_V2_MDD_20081114_211026_c3'};
Subjects(9)={'fl_V2_HC_20080618_210445_c4'};    
Subjects(10)={'fl_V2_MDD_20081204_210366p_c4'};
Subjects(11)={'fl_V2_HC_20080701_209489p_c1'};   
Subjects(12)={'fl_V2_MDD_20081219_211131_c3'};
Subjects(13)={'fl_V2_HC_20080722_210600p_c4'};   
Subjects(14)={'fl_V2_MDD_20090203_211205p_c3'};
Subjects(15)={'fl_V2_HC_20080814_210632p_c2'};   
Subjects(16)={'fl_V2_MDD_20090316_211289_c4'};
Subjects(17)={'fl_V2_HC_20080903_210627p_c1'};   
Subjects(18)={'fl_V2_MDD_20090319_211266_c1'};
Subjects(19)={'fl_V2_HC_20081203_211114p_c3'};   
Subjects(20)={'fl_V2_MDD_20090330_211370_c1'};
Subjects(21)={'fl_V2_HC_20090113_211132_c3'};    
Subjects(22)={'fl_V2_MDD_20090415_204042_c1'};
Subjects(23)={'fl_V2_HC_20090211_211174p_c1'};   
Subjects(24)={'fl_V2_MDD_20090422_211522p_c2'};
Subjects(25)={'fl_V2_HC_20090302_211040p_c2'};   
Subjects(26)={'fl_V2_MDD_20090423_204848_c3'};
Subjects(27)={'fl_V2_HC_20090313_211218_c3'};    
Subjects(28)={'fl_V2_MDD_20090427_211279p_c4'};
Subjects(29)={'fl_V2_HC_20090403_211342p_c2'};   
Subjects(30)={'fl_V2_MDD_20090603_211725p_c3'};
Subjects(31)={'fl_V2_HC_20090527_211643p_c1'};   
Subjects(32)={'fl_V2_MDD_20090611_211766_c1'};
Subjects(33)={'fl_V2_HC_20090826_207407p_c1'};   
Subjects(34)={'fl_V2_MDD_20090715_210938_c4'};
Subjects(35)={'fl_V2_HC_20090921_211947_c3'};    
Subjects(36)={'fl_V2_MDD_20090724_211945_c3'};
Subjects(37)={'fl_V2_HC_20091027_212306p_c2'};   
Subjects(38)={'fl_V2_MDD_20090915_212168_c1'};
Subjects(39)={'fl_V2_HC_20100217_212656p_c1'};   
Subjects(40)={'fl_V2_MDD_20090918_212130_c2'};
Subjects(41)={'fl_V2_HC_20100708_213493p_c1'};   
Subjects(42)={'fl_V2_MDD_20090930_212244_c4'};
Subjects(43)={'fl_V2_HC_20100817_213911_c3'};   
Subjects(44)={'fl_V2_MDD_20091123_211973p_c3'};
Subjects(45)={'fl_V2_HC_20100922_214212_c4'};    
Subjects(46)={'fl_V2_MDD_20100225_212759p_c2'};
Subjects(47)={'fl_V2_HC_20101002_214307_c1'};    
Subjects(48)={'fl_V2_MDD_20100325_212817p_c2'};
Subjects(49)={'fl_V2_HC_20101006_214248p_c2'};  
Subjects(50)={'fl_V2_MDD_20100401_212891p_c3'};
Subjects(51)={'fl_V2_HC_20101110_214469_c1'};    
Subjects(52)={'fl_V2_MDD_20100422_212997p_c3'};
Subjects(53)={'fl_V2_HC_20101119_214312_c2'};    
Subjects(54)={'fl_V2_MDD_20100811_213473p_c2'};
Subjects(55)={'fl_V2_HC_20101124_213955_c4'};    
Subjects(56)={'fl_V2_MDD_20100824_213198p_c1'};
Subjects(57)={'fl_V2_HC_20101201_214611_c2'};   
Subjects(58)={'fl_V2_MDD_20101019_210719_c2'};
Subjects(59)={'fl_V2_HC_20110111_214794p_c3'};  
Subjects(60)={'fl_V2_MDD_20101118_214482_c1'};
Subjects(61)={'fl_V2_HC_20110209_214834p_c1'};   
Subjects(62)={'fl_V2_MDD_20101209_214708p_c3'};
Subjects(63)={'fl_V2_MDD_20080508_208521p_c1'};  
Subjects(64)={'fl_V2_MDD_20110128_214858p_c3'};
Subjects(65)={'fl_V2_MDD_20080519_210392_c2'};
Subjects(66)={'fl_V2_MDD_20080728_210500_c3'};
Subjects(67)={'fl_V2_MDD_20080904_202962_c2'};


 spm_jobman('initcfg');

for i=44:67


SDirectory = {[char(directory1), '/', char(Subjects(i)), '/PPI_BR3_only/SPM.mat']};
matlabbatch{1}.spm.stats.con.spmmat = SDirectory;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'AngerShapePPI';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [0 0 0 0 0 0 1 0 0 0 0 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'FearShapePPI';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [0 0 0 0 0 0 0 1 0 0 0 -1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'HappyShapePPI';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [0 0 0 0 0 0 0 0 1 0 0 -1];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'SadShapePPI';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.convec = [0 0 0 0 0 0 0 0 0 0 1 -1];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'BR3regressor';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.convec = [0 0 0 0 0 0 0 0 0 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;


%spm('defaults', 'FMRI');
%spm_jobman('initcfg');
spm_jobman('serial', matlabbatch);
clear matlabbatch


end


