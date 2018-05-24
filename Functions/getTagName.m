function [tagName, DNA_binding] = getTagName(data, proteinIndex)

% reads the tag name out of the given protein if there is one

    tagName = '';

    if data.phosphorylation(proteinIndex) > 0
        tagName = 'P';
    elseif data.GXP(proteinIndex) > 0 
        if data.GXP(proteinIndex) == 1
            tagName = 'GDP';
        elseif data.GXP(proteinIndex) == 2
            tagName = 'GTP';
        end
    end
    
    if data.DNA_binding(proteinIndex) == 1
        DNA_binding = true;
    else
        DNA_binding = false;
    end

end