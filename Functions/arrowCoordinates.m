function coordinates = arrowCoordinates(center, direction)

% draws a arrow with the top at center
% direction = 1: arrow points up (upregulation)
% direction = -1: arrow pointsdown (downregulation)
% direction = 0: arrow is a square around the center (no change)

headWidth = 0.2;
headLength = 0.15;
f = 0.3;

if direction == 1
    original = [0 -headWidth -f*headWidth -f*headWidth f*headWidth f*headWidth headWidth 0; ...
        0 -headLength -headLength -1 -1 -headLength -headLength 0];
    coordinates = repmat(center',1,8) + original;
elseif direction == -1
    original = [0 headWidth f*headWidth f*headWidth -f*headWidth -f*headWidth -headWidth 0; ...
        0 +headLength +headLength 1 1 +headLength +headLength 0];
    coordinates = repmat(center',1,8) + original;
else
    original = [-2*f*headWidth 2*f*headWidth 2*f*headWidth -2*f*headWidth -2*f*headWidth; ...
        -f*headWidth -f*headWidth f*headWidth f*headWidth -f*headWidth];
    coordinates = repmat(center',1,5) + original;
end

end
    
