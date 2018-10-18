%% LOGIC MODEL CROSSTALK SIMULATION
%% 2018

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
addpath('Functions'); % path to functions that are used 
warning('off'); % Don't show all the warnings
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

%%% the file consists of 4 parts:
%     - Figure 2 B: comparison of the gene expression for deletion
%       strains under no glucose and glucose
%     - Figure 2 A: comparison of the gene expression for the WT, expected
%       gene expression vs. the model output with crosstalks and gapfilling
%       vs. the model without crosstalks and gapfilling
%       under no glucose and glucose
%     - Figure 2 C: compare the expression of HXK/HXT under Snf3/Rgt2
%       deletion for the model with crosstalk and gapfilling and without
%       to see how crosstalk can recapture the gene expression
%     - Simulation of all possible combinations and test which can
%       recapture gene expression after deletion of certain proteins 

% find all output pictures in the folder Data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 B: Compare gene expressions for different deletions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foldername = 'Data';
mkdir(foldername)

% define all knockouts we want to look at
knockouts1 = {};
knockouts2 = {'Snf1'};
knockouts3 = {'Rgt2', 'Snf3'};
knockouts4 = {'Tpk1', 'Tpk2', 'Tpk3', 'Bcy1'};

crosstalk = [0 0 0 0 0 0 1 0 1 0 0 0 0];
nitrogenLevel = 0;

% no glucose
% add 0.1 to make it more visibile in the plots
glucoseLevel = 0;
expressions1 = 0.1 + runLogicModel(glucoseLevel, nitrogenLevel, knockouts1, crosstalk);
expressions2 = 0.1 + runLogicModel(glucoseLevel, nitrogenLevel, knockouts2, crosstalk);
expressions3 = 0.1 + runLogicModel(glucoseLevel, nitrogenLevel, knockouts3, crosstalk);
expressions4 = 0.1 + runLogicModel(glucoseLevel, nitrogenLevel, knockouts4, crosstalk);

f = figure('Name', 'expressions', 'Position', [0, 0, 1000, 500]);
colors = [211 94 96;
    132 186 91;
    114 147 203;
    144 103 167] / 256;

s1 = subplot(2,1,1, 'Position', [0.13 0.63, 0.8, 0.3]);
set(s1, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar([expressions1' expressions2' expressions3' expressions4'],...
    'edgecolor','none','BarWidth', 1);
ylim([0 1.2]);
yticks([0.1 1.1]);
s1.YTickLabels = {'0','1'};
xticks([]);
set(s1, 'FontSize', 20);
title('no glucose', 'FontSize', 20);
box off;

% repeat everything for glucose
glucoseLevel = 1;
expressions1 = 0.1 + runLogicModel(glucoseLevel, nitrogenLevel, knockouts1, crosstalk);
expressions2 = 0.1 + runLogicModel(glucoseLevel, nitrogenLevel, knockouts2, crosstalk);
expressions3 = 0.1 + runLogicModel(glucoseLevel, nitrogenLevel, knockouts3, crosstalk);
expressions4 = 0.1 + runLogicModel(glucoseLevel, nitrogenLevel, knockouts4, crosstalk);

s2 = subplot(2,1,2, 'Position', [0.13 0.25, 0.8, 0.3]);
set(s2, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar([expressions1' expressions2' expressions3' expressions4'],...
    'edgecolor','none','BarWidth', 1);
ylim([0 1.2]);
yticks([0.1 1.1]);
s2.YTickLabels = {'0','1'};
s2.XTickLabels = {'SUC/GAL/MAL','HXT','HXK','STRE','PDS'};
set(s2, 'FontSize', 20);
title('glucose', 'FontSize', 20);
box off;

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1); 
text(0.08, 0.59, 'expression', 'FontSize', 25, 'Interpreter', 'Latex', ...
        'Rotation', 90, 'HorizontalAlignment', 'Center');
text(0.531, 0.15, 'genes', 'FontSize', 25, 'Interpreter', 'Latex', ...
         'HorizontalAlignment', 'Center');

legend(s2,{'WT', 'Snf1$\Delta$', 'Snf3$\Delta$Rgt2$\Delta$', ['Tpk1$\Delta$Tpk2$\Delta$' newline 'Tpk3$\Delta$Bcy1$\Delta$']},...
    'Position', [0.169 -0.054 0.719 0.149], 'FontSize', 20, 'Interpreter','Latex',...
    'Orientation', 'horizontal');
legend(s2,'boxoff')

hold off

print(f,[foldername, '/geneExpressions'], '-depsc', '-r0');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 A: Compare expected gene expressions with those from the simulation with or without gapfilling and crosstalk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foldername = 'Data';
mkdir(foldername)

% we don't look at nitrogen now
nitrogenLevel = [0 0];

% no glucose
% expectation from literature
expressionsExp = 0.1 + [1 0 0 1 1];
% model with knockout of all gapfillers and crosstalks
expressions1 = 0.1 + runLogicModel([1 0], nitrogenLevel, {'Xxx1','Xxx2','Xxx3','Xxx4','Xxx5','Xxx6'}, zeros(1,13));
% model with all gapfillers and crosstalks included
expressions2 = 0.1 + runLogicModel([1 0], nitrogenLevel, {}, ones(1,13));

f = figure('Name', 'expressions2', 'Position', [0, 0, 1000, 500]);
colors = [83 81 84;
    128 133 133;
    211 94 96] / 256;

s1 = subplot(2,1,1, 'Position', [0.13 0.63, 0.8, 0.3]);
set(s1, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar([expressionsExp' expressions1' expressions2'],...
    'edgecolor','none','BarWidth', 1);
ylim([0 1.2]);
yticks([0.1 1.1]);
s1.YTickLabels = {'0','1'};
s1.XTickLabels = {};
set(s1, 'FontSize', 20);
title('no glucose', 'FontSize', 20);
box off;

% repeat same experiment for glucose
expressionsExp = 0.1 + [0 1 1 0 0];
expressions1 = 0.1 + runLogicModel([0 1], nitrogenLevel, {'Xxx1','Xxx2','Xxx3','Xxx4','Xxx5'}, zeros(1,13));
expressions2 = 0.1 + runLogicModel([0 1], nitrogenLevel, {}, ones(1,13));

s2 = subplot(2,1,2, 'Position', [0.13 0.25, 0.8, 0.3]);
set(s2, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar([expressionsExp' expressions1' expressions2'],...
    'edgecolor','none','BarWidth', 1);
ylim([0 1.2]);
yticks([0.1 1.1]);
s2.YTickLabels = {'0','1'};
s2.XTickLabel = {'SUC/GAL/MAL','HXT','HXK','STRE','PDS'};
set(s2, 'FontSize', 20);
title('glucose', 'FontSize', 20);
box off;

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1); 
text(0.08, 0.51, 'expression', 'FontSize', 25, 'Interpreter', 'Latex', ...
        'Rotation', 90, 'HorizontalAlignment', 'Center');
text(0.531, 0.15, 'genes', 'FontSize', 25, 'Interpreter', 'Latex', ...
         'HorizontalAlignment', 'Center');

legend(s2,{'Expected',['WT without' newline 'crosstalk/gap filling'],['WT with' newline 'crosstalk/gap filling']},...
    'Position', [0.169 -0.054 0.719 0.149], 'FontSize', 20, 'Interpreter','Latex',...
    'Orientation', 'horizontal');
legend(s2,'boxoff')
    
print(f,[foldername, '/geneExpression2'], '-depsc', '-r0');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2 C: compare the expression of HXK/HXT under Snf3/Rgt2 deletion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foldername = 'Data';
mkdir(foldername)

nitrogenLevel = [0 0];

% glucose
expressionsExp = 0.1 + runLogicModel([0 1], nitrogenLevel, {}, [0 0 0 0 0 0 1 0 1 0 0 0 0]);
expressions1 = 0.1 + runLogicModel([0 1], nitrogenLevel, {'Snf3','Rgt2'}, [0 0 0 0 0 0 1 0 1 0 0 0 0]);
expressions2 = 0.1 + runLogicModel([0 1], nitrogenLevel, {'Snf3','Rgt2'}, [1 0 1 0 0 0 1 0 1 0 0 0 0]);

f = figure('Name', 'expressions3', 'Position', [0, 0, 1000, 500]);
colors = [211 94 96;
    114 147 203;
    57 106 177] / 256;

s2 = subplot(1,1,1,'Position', [0.13 0.329 0.8 0.55]);
set(s2, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar([(expressionsExp(2:3))' (expressions1(2:3))' (expressions2(2:3))'],...
    'edgecolor','none','BarWidth', 1);
ylim([0 1.2]);
yticks([0.1 1.1]);
s2.YTickLabels = {'0','1'};
s2.XTickLabel = {'HXT','HXK'};
set(s2, 'FontSize', 20);
title('glucose', 'FontSize', 20);
box off;
xlabel('genes', 'FontSize', 25, 'Interpreter', 'Latex', ...
         'HorizontalAlignment', 'Center');

ax1 = axes('Position',[0 0 1 1],'Visible','off');
axes(ax1); 
text(0.08, 0.6, 'expression', 'FontSize', 25, 'Interpreter', 'Latex', ...
        'Rotation', 90, 'HorizontalAlignment', 'Center');


legend(s2,{' WT','Snf3$\Delta$Rgt2$\Delta$',['Snf3$\Delta$Rgt2$\Delta$ with' newline 'additional crosstalk 1 and 3']},...
    'Position', [0.27 0 0.568 0.157], 'FontSize', 20, 'Interpreter','Latex',...
    'Orientation', 'horizontal');
legend(s2,'boxoff')
    
print(f,[foldername, '/geneExpressionHXT'], '-depsc', '-r0');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation of all possible combinations and test which can recapture gene expression after deletion of certain proteins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chose a deletion and a specific gene that we want to look at
% only the ones where deletion changes the gene expression are of interest

%knockouts = {'Snf1'};
knockouts = {'Rgt2', 'Snf3'};
%knockouts = {'Tpk1', 'Tpk2', 'Tpk3', 'Bcy1'};

gene = 3;
geneName = 'HXK';

crosstalks = createCrosstalkConfigurations(7);
crosstalks = [crosstalks(:,1:6) ones(2^7,1) crosstalks(:,7) ones(2^7,1) zeros(2^7,4)];
% crosstalk 7 and 9 are always on, and the other combinations are tested
% crosstalk 10-13 are from TOR are not included now

glucoseLevel = 1;
nitrogenLevel = 0;

expressionsWT = runLogicModel(glucoseLevel, nitrogenLevel, {}, [0 0 0 0 0 0 1 0 1 0 0 0 0]);
expressionsDeletion = runLogicModel(glucoseLevel, nitrogenLevel, knockouts, [0 0 0 0 0 0 1 0 1 0 0 0 0]);
expressions = runLogicModel(glucoseLevel, nitrogenLevel, knockouts, crosstalks);

% get the combinations that are similar to the wt and to deletion without
% crosstalk
% TOR crosstalks are neglected
% makes only sense when the gene expression changes due to the deletion
nActive = 2^7/2;

indWT = expressions(:,gene) == expressionsWT(gene);
crosstalksWT = crosstalks(indWT,[1,2,3,4,5,6,8]);
nWT = size(crosstalksWT,1);
activityWT = sum(crosstalksWT == 1)/nActive;

indDeletion = expressions(:,gene) == expressionsDeletion(gene);
crosstalksDeletion = crosstalks(indDeletion,[1,2,3,4,5,6,8]);
nDeletion = size(crosstalksDeletion,1);
activityDeletion = sum(crosstalksDeletion == 1)/nActive;

f = figure('Name', 'activity', 'Position', [0, 0, 700, 500]);
colors = [
    211 94 96;
%    132 186 91;
    114 147 203;
%    144 103 167
    ] / 256;

s1 = subplot(1,1,1,'Position',[0.13 0.306 0.777 0.617]);
set(s1, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
bar([activityWT' activityDeletion'], 'stacked', ...
    'edgecolor','none','BarWidth', 0.8);
yticks([0 1]);
xticks(1:7);
xlabel('active crosstalk')
ylabel('\% of crosstalk combinations','Interpreter', 'Latex')
s1.YTickLabel = {'0', '100'};
s1.XTickLabel = {'1','2','3','4','5','6','8'};

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

legend({'WT', name}, 'Position', [0.133 0.059 0.769 0.089], 'Orientation', 'horizontal')
legend('boxoff')
set(s1,'FontSize',24)
box off;

if glucoseLevel == 0
    savename = [foldername, '/', geneName, '-Gluc-',name2];
else
    savename = [foldername, '/', geneName, '+Gluc-',name2];
end

print(f,savename, '-depsc', '-r0');

    