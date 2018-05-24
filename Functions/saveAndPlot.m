function saveAndPlot(PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, foldername)
% save data in txt files and make plots to visualize

mkdir Data
mkdir('Data/', foldername)
path = ['Data/',foldername,'/'];

Name = [string('SUC2'), string('GAL1'), string('GAL4'), string('MAL61'), string('MAL62'), string('MAL63'), string('HXT1'), string('HXT2'), string('HXT3'), string('HXT4'), string('HXK2'), string('CTT1'),  string('DDR2'), string('HSP12'),string('RHR2'),string('PDC5'),string('PDC6'),string('LSC1'), string('NCR')]';
active = ones(length(Name),1);
Output = table(Name, active);

Output{1,2} = PromSite{1,2}; 
Output{2,2} =PromSite{2,2};
Output{3,2} =PromSite{2,2};
Output{4,2} =PromSite{3,2};
Output{5,2} =PromSite{3,2};
Output{6,2} =PromSite{3,2};
Output{7,2} =PromSite{7,2};
Output{8,2} =PromSite{4,2};
Output{9,2} =PromSite{4,2};
Output{10,2} =PromSite{4,2};
Output{11,2} =PromSite{8,2}; 
Output{12,2} =PromSite{5,2}; 
Output{13,2} =PromSite{5,2}; 
Output{14,2} =PromSite{5,2}; 
Output{15,2} =PromSite{6,2}; 
Output{16,2} =PromSite{6,2}; 
Output{17,2} =PromSite{6,2}; 
Output{18,2} =PromSite{6,2}; 
Output{19,2} =PromSite{9,2};

% save all outputs
writetable(Metabolites, [path,'Metabolites.txt'],'Delimiter','\t');
writetable(Miscl, [path,'Miscl.txt'],'Delimiter','\t');
writetable(Output, [path,'Output.txt'],'Delimiter','\t');
writetable(PKApw, [path,'PKApw.txt'],'Delimiter','\t');
writetable(PromSite, [path,'PromSite.txt'],'Delimiter','\t');
writetable(R2S3pw, [path,'R2S3pw.txt'],'Delimiter','\t');
writetable(Snf1pw, [path,'Snf1pw.txt'],'Delimiter','\t');
writetable(TORpw, [path,'TORpw.txt'],'Delimiter','\t');

% create and save picture of the corresponding output
disp('Create pictures...')
disp('  ')
dataNames = {'PKApw.txt', 'Snf1pw.txt', ...
    'R2S3pw.txt', 'Miscl.txt', 'TORpw.txt'};
mkdir(path,'Pictures'); 
for n = 1:length(dataNames)
    makeModelPicture(path, dataNames{n});
end

end