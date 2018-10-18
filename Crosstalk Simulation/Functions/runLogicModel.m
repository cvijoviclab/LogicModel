function output = runLogicModel(glucoseLevel, nitrogenLevel, knockouts, crosstalks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set glucose- and nitrogenlevel and knockouts. For each crosstalk
% configuration in crosstalks (matrix with each row being a sequence of 13
% binary numbers that indicate which crosstalk is active and which not

% all gene expression outputs are saved in a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nConfigurations, ~] = size(crosstalks);
output = size(nConfigurations, 5);

for i = 1:nConfigurations
    i
    
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

    tmpCrosstalks = crosstalks(i,:);
   
    nSteps = length(glucoseLevel);

    for j = 1:nSteps
        Metabolites{1,2} = glucoseLevel(j);
        Metabolites{4,2} = nitrogenLevel(j);
        [PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders] = ...
            runUntilSteadyState(PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, placeholders, tmpCrosstalks);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% OUTPUT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    output(i,1) = PromSite{1,2}; % SUC GAL MAL
    output(i,2) = PromSite{4,2}; % HXT
    output(i,3) = PromSite{8,2}; % HXK
    output(i,4) = PromSite{5,2}; % STRE
    output(i,5) = PromSite{6,2}; % PDS

end

end
