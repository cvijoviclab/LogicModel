function fig = plotEmptyCell()

% creates a figure with an empty cell (nucelus and membrane only basically)

    % Nucleus plotting
    cNuc = [0; 0]; % center of Nucleus
    rNuc = 1.25; % Radius of Nucleus
    % Define the coordinates of the nucleus
    [xNuc,yNuc] = ellipseMaker(cNuc,rNuc,rNuc);% Coordinates of nucleus, circular


    % Membrane plotting
    cMem = [0; 0]; % center of membrane
    rMem = 3; % Radius of membrane
    % Define the coordinates of the membrane
    [xMem,yMem] = ellipseMaker(cMem,rMem,rMem);% Coordinates of membrane, circular
    
    
    % The figure itself
    fig = figure('Name','Cell model','units','normalized','outerposition',...
        [0 0 1 1], 'Visible', 'off');
    plot(xMem, yMem, 'k')
    cytosolPlot = fill(xMem,yMem,[0.6 0.8 0.6], 'LineWidth', 1);
    
    hold on
    plot(xNuc, yNuc);
    nucleusPlot = fill(xNuc,yNuc,[0.6 0.6 0.6], 'LineWidth', 1);

    membranePlot = plot(xMem, yMem, 'Color', 'k', 'LineWidth', 7);
    
    % add line for the DNA bounded proteins
    plot([-1, 1], [0.1, 0.1], 'k', 'LineWidth', 2)
    
    % adjust and add legend
    axis([-1.6,1.6,0,3.7])
    daspect([1 1 1])
    set(gca,'xtick',[],'ytick',[])
%     legCell = legend([nucleusPlot,cytosolPlot,membranePlot],'  \textbf{Nucleus}', ...
%         '  \textbf{Cytosol}','  \textbf{Membrane}');
%     set(legCell, 'Interpreter', 'latex', 'FontSize', 24, 'Location', 'NorthEast',...
%         'AutoUpdate','off');
%     legend boxoff

end
