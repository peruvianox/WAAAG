function [KinematicsOptions, KinematicsPlotOptions, KineticsOptions, KineticsPlotOptions] = SelectOptions

% % -------------------------------------------------------------------------
% Selects Options from Selector_File
% -------------------------------------------------------------------------
% Author: Ricky Pimentel
% Date: December 21, 2016
%
% Description:	This function will read in the options from the
% selector_file_template and define the user inputs for kinemataics and
% kinetics data processing and plotting.
%
% Input:      Selector_File_Template.xlsx 
%
% Outputs: (all structural data)
%           KinematicsOptions           Options for Kinematics Processing
%           KinematicsPlotOptions     Options for Kinematics Plots
%           KineticsOptions                Options for Kinetics Processing
%           KineticsPlotOptions          Options for Kinetics Plots
%
% % DEBUG LOOP %%%%%%%%%
% close all
% clear
% clc
% %%% END DEBUG LOOP %%%%%

%% Load Selector_File_Template

[~,~,RAW] = xlsread('Selector_File.xlsx','Sheet1');

%% Find Regions for each output
% Define start of regions
Ind = strcmp(RAW(:,1),'Kinematics Data Processing Options') == 1;
KinemaStart = find(Ind);

Ind = strcmp(RAW(:,1),'Kinematics Plotting Options') == 1;
KinemaPlotStart = find(Ind);

Ind = strcmp(RAW(:,1),'Kinetics Data Processing Options') == 1;
KineticStart = find(Ind);

Ind = strcmp(RAW(:,1),'Kinetics Plotting Options') == 1;
KineticPlotStart = find(Ind);
% Define end of regions
KinemaEnd = KinemaPlotStart - 1;
KinemaPlotEnd = KineticStart -1;
KineticEnd = KineticPlotStart - 1;
Ind = strcmp(RAW(:,1),'END') == 1;
KineticPlotEnd = find(Ind) - 1; 

%% Define Kinematics Options
for i = 1: KinemaEnd - KinemaStart
    KinematicsOptions(i).Input = RAW(KinemaStart + i,1);
    KinematicsOptions(i).Choice = RAW(KinemaStart + i,2);
end

%% Define Kinematic Plot Options
for i = 1: KinemaPlotEnd - KinemaPlotStart
    KinematicsPlotOptions(i).Input = RAW(KinemaPlotStart + i,1);
    KinematicsPlotOptions(i).Choice = RAW(KinemaPlotStart + i,2);
end

%% Define Kinetics Options
for i = 1: KineticEnd - KineticStart
    KineticsOptions(i).Input = RAW(KineticStart + i,1);
    KineticsOptions(i).Choice = RAW(KineticStart + i,2);
end

%% Define Kinetic Plot Options
for i = 1: KineticPlotEnd - KineticPlotStart
    KineticsPlotOptions(i).Input = RAW(KineticPlotStart + i,1);
    KineticsPlotOptions(i).Choice = RAW(KineticPlotStart + i,2);
end

end


