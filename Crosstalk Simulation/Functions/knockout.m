function [PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw] = ...
    knockout(PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw, knockouts)
% modifies the tables and sets all protein/genes in knockouts to presence =
% 0
% knockout is a cell of strings which specifies the genes/proteins that
% are knocked out

nKnockouts = length(knockouts);

for i = 1:nKnockouts
    proteinName = knockouts{i};
    
    % find position of protein in data files
    for j = {PromSite, Metabolites, Miscl, Snf1pw, R2S3pw, PKApw, TORpw}
        players = j{1}.Name;
        pos = find(proteinName == players, 1);
        if ~isempty(pos)
            break;
        end
    end
    
    % change presence/activity in corresponding table to 0
    if isequal(j{1},PromSite)
        PromSite{pos,2} = 0;
    elseif isequal(j{1},Metabolites)
        Metabolites{pos,2} = 0;
    elseif isequal(j{1},Miscl)
        Miscl{pos,2} = 0;
    elseif isequal(j{1},Snf1pw)
        Snf1pw{pos,2} = 0; 
    elseif isequal(j{1},R2S3pw)
        R2S3pw{pos,2} = 0; 
    elseif isequal(j{1},PKApw)
        PKApw{pos,2} = 0; 
    elseif isequal(j{1}, TORpw)
        TORpw{pos,2} = 0;
    end
end

end