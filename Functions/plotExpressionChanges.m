function plotExpressionChanges(path, gene, changesExperimental, changesLogical, name)

% plots the expression changes for two data sets 
%
% changesExperimental: experimental data
% changesLogical: data from logical model 
% (the two vectors should have the same length)
%
% name: name of the file that is saved

nData = length(changesExperimental);

% create figures
fig = figure('units','normalized','outerposition',...
    [0 0 1 0.4], 'Visible', 'off');
plot(linspace(0,nData+1), zeros(1,100), 'k--')
axis([0 nData+1 -1 1]);
xticks(1:nData+1)
set(gca,'TickLabelInterpreter', 'tex')
xticklabels(gene);
set(gca,'YTick',[])
set(gca,'fontsize',17)

hold on

% go through all proteins and add arrows to the plot for both models
for i = 1:nData
    arrow = arrowCoordinates([i,changesExperimental(i)],changesExperimental(i));
    f = fill(arrow(1,:), arrow(2,:),[0.6 0.8 0.6], 'LineWidth', 1,...
        'EdgeColor', [0.6 0.8 0.6]);
    arrowLogical = arrowCoordinates([i,changesLogical(i)],changesLogical(i));
    p = plot(arrowLogical(1,:),arrowLogical(2,:),'k','LineWidth', 1.5); 
end

%add legend and title
leg = legend([f,p],'   experiment', ...
        '  logical model');
legend boxoff 
set(leg, 'Interpreter', 'latex', 'FontSize', 20, 'Location', 'EastOutside');
title(name, 'FontSize', 20,'Interpreter' , 'latex')  
hold off

% save the picture
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto', 'PaperSize',[pos(3),pos(4)]);
print(gcf,[path,'/Data/', name], '-depsc', '-r0')

close(fig)

end