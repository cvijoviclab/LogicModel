function [protein, noLow, lowHigh, noHigh] = processRNAData()

% reads the data in RNAData and produces output that is an array containing a
% -1 0 or 1 for each protein (rows), indicating if its up or downregulated
% when changing the glucose concentration from nno to low, from low to high
% and from no to high
% if there is no data NaN

RNAdata = readtable('RNAData.txt','Delimiter','\t');
RNAdata = table2struct(RNAdata, 'ToScalar', true);

protein = RNAdata.Protein; % just the names of the proteins/genes
noGlucoseExp = RNAdata.noGl;
lowGlucoseExp = RNAdata.lowGl;
highGlucoseExp = RNAdata.highGl;

for i = 1:size(protein,1)
    % for every protein check if it is up- or downregulated by changing the
    % glucose concentration from 0 to low, from low to high and from 0 to high
    no = noGlucoseExp(i);
    low = lowGlucoseExp(i);
    high = highGlucoseExp(i);
    
    noLow(i) = -sign(no-low);
    lowHigh(i) = -sign(low-high);
    noHigh(i) = -sign(no-high);
end

end