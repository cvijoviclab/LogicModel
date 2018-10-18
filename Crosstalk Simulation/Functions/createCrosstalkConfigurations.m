function configs = createCrosstalkConfigurations(nCrosstalks)

% create vector with all possible combinations of active(1) or inactive(0) crosstalks
% (2^nCrosstalks possibilities) 

n = 2^nCrosstalks-1;
configs = zeros(n, nCrosstalks);

for i = 1:n
    x = dec2bin(i);
    for j = 1:length(x)
        % take i+1 such that first row is zeros only
        configs(i+1, nCrosstalks+1-j) = str2double(x(length(x)+1-j));
    end
end

end