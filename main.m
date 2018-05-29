%% LOGIC MODEL 
%% 2018

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
addpath('Functions'); % path to functions that are used 
warning('off'); % Don't show all the warnings

% RULES 
% metabolite(glucose level: no(0) high(1), external on/off, nitrogen level: no(0) high(1))
% proteinname(presence, localization nucleus(2); cytosol(1); membrane(0); phosphorylation, guanylation non(0); GDP(1);GTP(2), dna bound)
% genes(promotor on/off)

%% Produce data for a defined glucose, nitrogen combination

% choose the folder name where the data is saved
foldername = '-NitrTo+Nitr-Gluc+Crosstalk';
% define all knockouts, i.e. {'Snf1', 'Tor1'}
knockouts = {};
% define all active crosstalks out of 13, i.e. 
% crosstalks = [0 0 0 0 0 0 0 0 0 0 0 0 0]; 
crosstalks = [1 1 1 1 1 1 1 1 1 1 1 1 1]; 

% set glucose and nitrogen level sequence
% both have to have the same length, and each entry corresponds to one
% 'time step'
glucoseLevels = [0 0];
nitrogenLevels = [0 1];

% run logic model and create txt files and pictures with data in the
% defined folder
runLogicModel(glucoseLevels, nitrogenLevels, knockouts, crosstalks, foldername);

    