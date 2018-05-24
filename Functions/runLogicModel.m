function runLogicModel(glucoseLevel, nitrogenLevel, knockouts, crosstalks, foldername)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% glucoseLevel is a vector giving a sequence of glucose levels, 
% e.g. [0 1 0 2], with
% 0: no glucose
% 1: low glucose
% 2: high glucose

% nitrogenLevel is a vector giving a sequence of nitrogen levels, 
% e.g. [0 1 1 0], with
% 0: no nitrogen
% 1: nitrogen
% it has to have the same length as glucoseLevel

% knockouts is a cell of strings which contain the names of proteins that
% are knocked out, e.g. {'Snf3', 'Std1'}

% filename specifies the folder where the results are saved. If it does
% not exist it will be created, e.g. 'MyFoldername'

% After initialization the model will loop through the vector and runs the
% model. The output of each round will be the input of the next round.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALISATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, ...
    placeholders] = initalizeModel();

% change knockouts
[PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw] = ...
    knockout(PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, knockouts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

levels = {'zero', 'low', 'high'}; % only for displaying
levels2 = {'zero', 'high'}; % only for displaying
nSteps = length(glucoseLevel);

for i = 1:nSteps
    
    disp(['STEP ', num2str(i),': GLUCOSE LEVEL ', levels{glucoseLevel(i)+1}, ', NITROGEN LEVEL ', levels2{nitrogenLevel(i)+1}])
    Metabolites{1,2} = glucoseLevel(i);
    Metabolites{4,2} = nitrogenLevel(i);
    [PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders] = ...
        runUntilSteadyState(PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders, crosstalks);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% OUTPUT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nSteps == 1
        saveAndPlot(PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, ...
            foldername);
    else
        saveAndPlot(PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, ...
            [foldername, '/step', num2str(i)]);
    end
    
end


end
