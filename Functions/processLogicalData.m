function [protein, logicalNoLow, logicalLowHigh, logicalNoHigh] = processLogicalData()

logicalDataNo = readtable('Data/NoGlucose/Output.txt','Delimiter','\t');
logicalDataNo = table2struct(logicalDataNo, 'ToScalar', true);
logicalDataLow = readtable('Data/LowGlucose/Output.txt','Delimiter','\t');
logicalDataLow = table2struct(logicalDataLow, 'ToScalar', true);
logicalDataHigh = readtable('Data/HighGlucose/Output.txt','Delimiter','\t');
logicalDataHigh = table2struct(logicalDataHigh, 'ToScalar', true);

logicalNoLow = - logicalDataNo.active + logicalDataLow.active ;
logicalLowHigh = - logicalDataLow.active  + logicalDataHigh.active ;
logicalNoHigh = - logicalDataNo.active  + logicalDataHigh.active;

% name is the same for all three types
protein = logicalDataNo.Name;

end