function [cNucleus, cNucleusDNAb, cCytosol, cMembrane] = ...
    getProteinPositions(nNucleus, nNucleusDNAb, nCytosol, nMembrane)

% distribute positions depending on how many proteins are in the
% nucelus/cytosol/membrane
    
    % the ones that are membrane bound are just distributed along the
    % circular sector with angular difference of dPhi
    % depending on how many boxes are drawn
    % max approximately 8 membrane bound proteins possible
    if nMembrane > 8
        disp('too many membrane bound proteins: cannot be displayed correctly')
    end
    
    if nMembrane < 4
        dPhi = pi/25; 
    else
        dPhi = pi/(25+(nMembrane-4)*5);
    end
    
    if mod(nMembrane,2) == 0
        vec = pi/2 + dPhi * [-fliplr(1:2:nMembrane), 1:2:nMembrane];
        cMembrane = fliplr(3 * [cos(vec); sin(vec)]);
    else
        vec = pi/2 + dPhi * [-fliplr(2:2:nMembrane), 0:2:nMembrane];
        cMembrane = fliplr(3 * [cos(vec); sin(vec)]);
    end
    
    
    % proteins that are in the nucleus and DNA bound are just positionned
    % along the line in the nucleus
    % maximal 4 proteins find space
    maxN = 4;
    if nNucleusDNAb > maxN
        disp('too many DNA bound proteins in the nucleus: cannot be displayed correctly')
    end
   
    cNucleusDNAb = [];
    for i = 1:nNucleusDNAb
        newPosition = [-0.65; 0.3] + (i-1) * [0.65; 0];
        cNucleusDNAb = [cNucleusDNAb newPosition];
    end
    
    
    % proteins that are in the nucleus and not DNA bound are positionned by 
    % defined coordinates (brute force)
    % maximal 5 proteins find space
    maxN = 5;
    if nNucleus > maxN
        disp('too many proteins in the nucleus: cannot be displayed correctly')
        positions{maxN+1:nNucleus} = deal([[-0.25 1; 0.25 1; -0.5 0.65; ...
            0 0.65; 0.5 0.65]'; repmat([0 0.75]',1,nNucleus-maxN)]);
        
    elseif nNucleus > 0 && nNucleus <= maxN
        positions{1} = [0 0.75]';
        positions{2} = [-0.3 0.75; 0.3 0.75]';
        positions{3} = [0 1; -0.25 0.65; 0.25 0.65]';
        positions{4} = [0 1; -0.5 0.65; 0 0.65; 0.5 0.65]'; 
        positions{5} = [-0.25 1; 0.25 1; -0.5 0.65; 0 0.65; 0.5 0.65]'; 

        cNucleus = positions{nNucleus};
        
    else
        cNucleus = [];
        
    end
    
    clear positions;
    
    % proteins that are in the cytosol are positionned by defined
    % coordinates (brute force)
    % maximal 17 proteins find space
    maxN = 17;
    if nCytosol > maxN
        disp('too many proteins in the cytosol: cannot be displayed correctly')
        positions{15:nCytosol} = deal([[-0.25 1; 0.25 1; -0.5 0.65; ...
            0 0.65; 0.5 0.65]'; repmat([0 2]',1,nCytosol-maxN)]);
        
    elseif nCytosol > 0 && nCytosol <= maxN
        positions{1} = [0 2]';
        positions{2} = [-0.5 2; 0.5 2]';
        positions{3} = [0 2.3; -0.5 1.8; 0.5 1.8]';
        positions{4} = [-0.5 2.3; 0.5 2.3; -0.5 1.8; 0.5 1.8]'; 
        positions{5} = [-0.3 2.3; 0.3 2.3; -0.6 1.8; 0 1.8; 0.6 1.8]'; 
        positions{6} = [-0.6 2.3; 0 2.3; 0.6 2.3; -0.6 1.8; 0 1.8; 0.6 1.8]';
        positions{7} = [0 2.5; -0.6 2.1; 0 2.1; 0.6 2.1; -0.6 1.7; 0 1.7; 0.6 1.7]';
        positions{8} = [-0.3 2.5; 0.3 2.5; -0.6 2.1; 0 2.1; 0.6 2.1; -0.6 1.7; ...
            0 1.7; 0.6 1.7]';
        positions{9} = [-0.6 2.5; 0 2.5; 0.6 2.5; -0.6 2.1; 0 2.1; 0.6 2.1; ...
            -0.6 1.7; 0 1.7; 0.6 1.7]';
        positions{10} = [-0.3 2.5; 0.3 2.5; -0.9 2.1; -0.3 2.1; 0.3 2.1; ...
             0.9 2.1; -0.9 1.7; -0.3 1.7; 0.3 1.7; 0.9 1.7]';
        positions{11} = [-0.6 2.5; 0 2.5; 0.6 2.5; -0.9 2.1; -0.3 2.1; 0.3 2.1; ...
             0.9 2.1; -0.9 1.7; -0.3 1.7; 0.3 1.7; 0.9 1.7]';
        positions{12} = [-0.3 2.5; 0.3 2.5; -1.2 2.1; -0.6 2.1; 0 2.1; 0.6 2.1; ...
             1.2 2.1; -1.2 1.7; -0.6 1.7; 0 1.7; 0.6 1.7; 1.2 1.7]';
        positions{13} = [-0.6 2.5; 0 2.5; 0.6 2.5; -1.2 2.1; -0.6 2.1; 0 2.1; ...
            0.6 2.1; 1.2 2.1; -1.2 1.7; -0.6 1.7; 0 1.7; 0.6 1.7; 1.2 1.7]';
        positions{14} = [-0.9 2.5; -0.3 2.5; 0.3 2.5; 0.9 2.5; -1.2 2.1; -0.6 2.1; ...
            0 2.1; 0.6 2.1; 1.2 2.1; -1.2 1.7; -0.6 1.7; 0 1.7; 0.6 1.7; 1.2 1.7]';
        positions{15} = [-0.3 2.65; 0.3 2.65; -0.6 2.3; 0 2.3; 0.6 2.3; -1.2 1.95; ...
            -0.6 1.95; 0 1.95; 0.6 1.95; 1.2 1.95; -1.2 1.6; -0.6 1.6; 0 1.6; ...
            0.6 1.6; 1.2 1.6]';
        positions{16} = [-0.3 2.65; 0.3 2.65; -0.9 2.3; -0.3 2.3; 0.3 2.3; 0.9 2.3; ...
            -1.2 1.95; -0.6 1.95; 0 1.95; 0.6 1.95; 1.2 1.95; -1.2 1.6; -0.6 1.6; ...
            0 1.6; 0.6 1.6; 1.2 1.6]';
        positions{17} = [-0.6 2.65; 0 2.65; 0.6 2.65; -0.9 2.3; -0.3 2.3; 0.3 2.3; ...
             0.9 2.3; -1.2 1.95; -0.6 1.95; 0 1.95; 0.6 1.95; 1.2 1.95; -1.2 1.6; ...
             -0.6 1.6; 0 1.6; 0.6 1.6; 1.2 1.6]';

        cCytosol = positions{nCytosol};   
        
    else
        cCytosol = [];
        
    end
   
end