function includeProtein(center, name, addition, DNAb)

% adds a rectangle with a proteins name ('MIG1', 'SNF1', ..) in the current figure
% addition can be an empty string
% but if it is not there will be a small ellipse on the top right corner
% with the string ('P', 'GTP', ..) in it
% if DNAb is true a line beneath the box is plotted as well

    % define size of boxes for proteins
    length = 0.4*1.3;
    width = 0.16*1.3;
    fontsize = 28;

    % define size and fontsize of addition
    l = strlength(addition);
    if l == 1
        semiaxisMajor = 0.07;
        semiaxisMinor = 0.07;
    else
        semiaxisMajor = 0.13;
        semiaxisMinor = 0.07;
    end
    fontsizeAdd = 18;

    % plot 
    edges = rectangleMaker(center, length, width);
    hold on
    plot(edges(:,1), edges(:,2), 'k')
    fill(edges(:,1), edges(:,2), 'w', 'LineWidth', 1);
    text(center(1), center(2), ['\textbf{',name,'}'],'FontSize', fontsize, ...
        'HorizontalAlignment', 'center', 'Interpreter' , 'latex')
    
    % add additional tag is added on the top right corner 
    % if it is not an empty string
    if ~isempty(addition) 
        upperCorner = edges(1,:);
        [xAdd, yAdd] = ellipseMaker(upperCorner + [-0.04 0.04], semiaxisMajor, semiaxisMinor);
        plot(xAdd, yAdd, 'k')
        fill(xAdd, yAdd, [0.9 0.9 0.9], 'LineWidth', 1)
        text(upperCorner(1) -0.04, upperCorner(2) +0.03, ['\textbf{',addition,'}'],'FontSize', fontsizeAdd, ...
        'HorizontalAlignment', 'center', 'Interpreter' , 'latex')
    end
    
    % add additional line if DNA bound
    if DNAb
        plot([edges(3,1), edges(3,1)], [edges(3,2), edges(3,2) - 0.09], 'k', 'LineWidth', 1)
    end
    
    hold off
end