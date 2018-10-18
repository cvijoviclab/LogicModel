%% LOGIC MODEL CROSSTALK SIMULATION WITH TOR
%% 2018

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
addpath('Functions'); % path to functions that are used 
warning('off'); % Don't show all the warnings
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

%%% the file consists of 2 parts:
%     - Figure 4 B: comparison of the gene expression for the WT, expected
%       gene expression vs. the model output with crosstalks and gapfilling
%       vs. the model without crosstalks and gapfilling
%       under no nitrogen and nitrogen and no glucose and glucose
%     - Simulation of all possible combinations and test which can
%       recapture gene expression after deletion of certain proteins 

% find all output pictures in the folder Data/TOR


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare gene expressions with crosstalk and gapfilling and without
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foldername = 'Data/TOR';
mkdir(foldername)

result = [];

% no nitrogen
glucoseLevel = [0 0];
expressions1 = 0.1 + runLogicModelTOR(glucoseLevel, [1 0], {'Xxx1','Xxx2','Xxx3','Xxx4','Xxx5','Xxx6'}, zeros(1,13));
expressions2 = 0.1 + runLogicModelTOR(glucoseLevel, [1 0], {}, ones(1,13));

result = [result; [1.1 expressions1(6) expressions2(6)]];

glucoseLevel = [1 1];
expressions1 = 0.1 + runLogicModelTOR(glucoseLevel, [1 0], {'Xxx1','Xxx2','Xxx3','Xxx4','Xxx5','Xxx6'}, zeros(1,13));
expressions2 = 0.1 + runLogicModelTOR(glucoseLevel, [1 0], {}, ones(1,13));
result = [result; [1.1 expressions1(6) expressions2(6)]];

% nitrogen
glucoseLevel = [0 0];
expressions1 = 0.1 + runLogicModelTOR(glucoseLevel, [0 1], {'Xxx1','Xxx2','Xxx3','Xxx4','Xxx5','Xxx6'}, zeros(1,13));
expressions2 = 0.1 + runLogicModelTOR(glucoseLevel, [0 1], {}, ones(1,13));

result = [result; [1.1 expressions1(6) expressions2(6)]];

glucoseLevel = [1 1];
expressions1 = 0.1 + runLogicModelTOR(glucoseLevel, [0 1], {'Xxx1','Xxx2','Xxx3','Xxx4','Xxx5','Xxx6'}, zeros(1,13));
expressions2 = 0.1 + runLogicModelTOR(glucoseLevel, [0 1], {}, ones(1,13));
result = [result; [0.1 expressions1(6) expressions2(6)]];

f = figure('Name', 'expressions2', 'Position', [0, 0, 700, 500]);
colors = [83 81 84;
    128 133 133;
    211 94 96] / 256;

s1 = subplot(1,1,1, 'Position', [0.13 0.63, 0.8, 0.3]);
set(s1, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar(result(1:2,:),...
    'edgecolor','none','BarWidth', 1);
ylim([0 1.2]);
yticks([0.1 1.1]);
s1.YTickLabels = {'0','1'};
s1.XTickLabels = {};
set(s1, 'FontSize', 20);
title('no nitrogen', 'FontSize', 20);
box off;

s2 = subplot(2,1,2, 'Position', [0.13 0.25, 0.8, 0.3]);
set(s2, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar(result(3:4,:),...
    'edgecolor','none','BarWidth', 1);
ylim([0 1.2]);
yticks([0.1 1.1]);
s2.YTickLabels = {'0','1'};
s2.XTickLabel = {'NCR -Glc', 'NCR +Glc'};
set(s2, 'FontSize', 20);
title('nitrogen', 'FontSize', 20);
box off;

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1); 
text(0.08, 0.59, 'expression', 'FontSize', 25, 'Interpreter', 'Latex', ...
        'Rotation', 90, 'HorizontalAlignment', 'Center');
text(0.531, 0.16, 'genes', 'FontSize', 25, 'Interpreter', 'Latex', ...
         'HorizontalAlignment', 'Center');

legend(s2,{'Expected',['WT without' newline 'crosstalk/gap filling'],['WT with' newline 'crosstalk/gap filling']},...
    'Position', [0.169 -0.054 0.719 0.149], 'FontSize', 20, 'Interpreter','Latex',...
    'Orientation', 'horizontal');
legend(s2,'boxoff')
    
print(f,[foldername, '/geneExpression'], '-depsc', '-r0');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Produce gene expression data for crosstalk combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chose a deletion and a specific gene that we want to look at
% only the ones where deletion changes the gene expression are of interest

knockouts = {'Snf1'};
%knockouts = {'Rgt2', 'Snf3'};
%knockouts = {'Tpk1', 'Tpk2', 'Tpk3', 'Bcy1'};

gene = 1;
geneName = 'SUC-GAL-MAL';

crosstalks = createCrosstalkConfigurations(11);
crosstalks = [crosstalks(:,1:6) ones(2^11,1) crosstalks(:,7) ones(2^11,1) crosstalks(:,8:11)];
% crosstalk 7 and 9 are always on, and the other combinations are tested
% crosstalk 10-13 are from TOR are included now

glucoseLevel = 0;
nitrogenLevel = 0;

expressionsWT = runLogicModel(glucoseLevel, nitrogenLevel, {}, [0 0 0 0 0 0 1 0 1 0 0 0 0]);
expressionsDeletion = runLogicModel(glucoseLevel, nitrogenLevel, knockouts, [0 0 0 0 0 0 1 0 1 0 0 0 0]);
expressions = runLogicModel(glucoseLevel, nitrogenLevel, knockouts, crosstalks);

% get the combinations that are similar to the wt and to deletion without
% crosstalk
% TOR crosstalks INCLUDED
% makes only sense when something changes
nActive = 2^11/2;

indWT = expressions(:,gene) == expressionsWT(gene);
crosstalksWT = crosstalks(indWT,[1,2,3,4,5,6,8,10,11,12,13]);
nWT = size(crosstalksWT,1);
activityWT = sum(crosstalksWT == 1)/nActive;

indDeletion = expressions(:,gene) == expressionsDeletion(gene);
crosstalksDeletion = crosstalks(indDeletion,[1,2,3,4,5,6,8,10,11,12,13]);
nDeletion = size(crosstalksDeletion,1);
activityDeletion = sum(crosstalksDeletion == 1)/nActive;

f = figure('Name', 'activity', 'Position', [0, 0, 700, 500]);
colors = [
    211 94 96;
%    132 186 91;
%    114 147 203;
    144 103 167
    ] / 256;

s1 = subplot(1,1,1);
set(s1, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar([activityWT' activityDeletion'], 'stacked', ...
    'edgecolor','none','BarWidth', 0.8);
yticks([0 1]);
xticks(1:11);
xlabel('active crosstalk')
ylabel('\% of crosstalk combinations','Interpreter', 'Latex')
s1.YTickLabel = {'0', '100'};
s1.XTickLabel = {'1','2','3','4','5','6','8','10','11','12','13'};

if glucoseLevel == 0
    titlename = [geneName, ' for no glucose'];
else
    titlename = [geneName, ' for glucose'];
end
title(titlename);
name = [];
name2 = [];
for i = 1:length(knockouts)
    name = [name, knockouts{i}, '$\Delta$'];
    name2 = [name2, knockouts{i}];
end
legend('WT', name, 'Location', 'NorthEast')
set(s1,'FontSize',24)
box off;

mkdir Data TOR
if glucoseLevel == 0
    savename = [foldername, '/', geneName, '-Gluc-',name2];
else
    savename = [foldername, '/', geneName, '+Gluc-',name2];
end

print(f,savename, '-depsc', '-r0');