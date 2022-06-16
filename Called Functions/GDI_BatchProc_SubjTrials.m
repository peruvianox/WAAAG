function [Representative] = GDI_BatchProc_SubjTrials(files2get)
% -------------------------------------------------------------------------
% GAIT DEVIATION INDEX - CONTROL DATABASE BATCH PROCESSING
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: July 14, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This function batch processes all the subject's trials and
%               determines most representative gait cycle.
% 
% Input:        C3Dfile(s)      user defined C3D file(s) to be analyzed
%               
%
% Output:       Representative  structure with representative gait cycle
%                               and corresponding trial name.
%   
warning off
% %DEBUG LOOP
% close all
% clear all
% clc
% 
% [files2get, pathname] = uigetfile('*.c3d', 'Please select all C3D file(s) to be analyzed',...
%                                        'MultiSelect', 'on');
% 
% % %%% END DEBUG LOOP %%%%%

%% Batch Processes All Files
% files2get = files2get';  % must be enabled when running with GUI!
if iscell(files2get) == 0
    NumTrials = 1;
else
    NumTrials = length(files2get);
end

if NumTrials > 1
    % if multiple files are selected
    for n = 1:size(files2get, 2)
        % converts name of each file from cell structure to char
        FullFileName = char(files2get(n));
        % Analyzes current C3D file
        [variance_ratios, residuals, ave_kinematics, stddev_kinematics, kinematics, TempSpat] = GDI_AnalyTrial(FullFileName);
        % builds structure of current trial's name
        c3d_files(n).name = FullFileName;
        % 9 kinematic curves
        c3d_files(n).kinematics = kinematics;
        % mean kinematic curves is a structure of 9 kinematic curves
        c3d_files(n).ave_kinematics = ave_kinematics;
        % standard deviation of kinematic curves
        c3d_files(n).stddev_kinematics = stddev_kinematics;
        % structure of each trial's temporal spatial variables
        c3d_files(n).tempspat = TempSpat;
        % structure of trial's variance ratios (12x1)
        c3d_files(n).variance_ratios = variance_ratios.trial;
        % ordered structure of trial's representative gait cycles for all left & right gait cycles
        c3d_files(n).residuals = residuals;
                
    end
else
    % if only one file is selected
    % FullFileName = char(cell2mat(files2get));
    FullFileName = char(files2get);
    % Analyzes current C3D file
    [variance_ratios, residuals, ave_kinematics, stddev_kinematics, kinematics, TempSpat] = GDI_AnalyTrial(FullFileName);
    % builds structure of current trial's name
    c3d_files.name = FullFileName;
    % 9 kinematic curves
    c3d_files.kinematics = kinematics;
    % mean kinematic curves is a structure of 9 kinematic curves
    c3d_files.ave_kinematics = ave_kinematics;
    % standard deviation of kinematic curves
    c3d_files.stddev_kinematics = stddev_kinematics;
    % structure of each trial's temporal spatial variables
    c3d_files.tempspat = TempSpat;
    % structure of trial's variance ratios (12x1)
    c3d_files.variance_ratios = variance_ratios.trial;
    % ordered structure of trial's representative gait cycles for all left & right gait cycles
    c3d_files.residuals = residuals;
end



%% Determines Representative Trial
% Determines the temporal spatial residuals
[sum_resids, c3d_files] = GAMS_TempSpatResids(c3d_files);

for t = 1:length(c3d_files)
    % sorts each trial's variance ratios
    var_ratios(:,t) = c3d_files(t).variance_ratios;
    
    % sorts each trial's tempspat residuals
    TS_resids(:,t) = c3d_files(t).residTS;
    
    % summed total of residuals and variance ratios
    trial_variance(:,t) =  var_ratios(:,t);
end
% ranks the representative trails, in ascending order, left = 1, right = 2
for i = 1:2
    if i == 1
        [variances, trial] = sort(trial_variance.left, 'ascend');
        Representative.trial_order.left = trial;
    else
        [variances, trial] = sort(trial_variance.right, 'ascend');
        Representative.trial_order.right = trial;
    end
end
% replicates matrix of trial names from structures format
rep_trials = (repmat({c3d_files.name}, 1,1));

% sorts names of trials by rank order
for t = 1:length(rep_trials)
    ordered_rep_trials(t) = rep_trials(trial(t));
end

% stores re-ordered trials in structure
% trial listed best to worst
Representative.trials = ordered_rep_trials;

% stores all c3d files' data in structure Representative
Representative.c3d_files = c3d_files;
