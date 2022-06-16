function [GM] = GAMSmini_BatchGaitMeasures(files2get, age)
% -------------------------------------------------------------------------
% GAMS MINI PROGRAM - BATCH PROCESSING OF GAIT MEASURES
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: April 22, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This function batch processes all c3d files inputted
%               and returns the gait measures: overall,
%               spatial, and temporal for the GAMS_GUI.
%
% Input:        files2get        user selected C3D file(s) to be analyzed
%               Ages             consists of age for each c3d file (aka
%                                subject)
%
%
% Output:       GM               structure with each file's gait measures
%
%

%% DEBUG LOOP
% close all
% clear all
% clc
%
if exist('files2get','var') == 0
    [files2get, ~] = uigetfile('*.c3d*', 'Please select all C3D file(s) to be analyzed', 'MultiSelect', 'on');
end

%% Batch Processes All Files
if exist('age','var') == 0
    agechr = inputdlg('What is the age of the subject?');
    age = str2num(char(agechr));
end
    
if iscell(files2get)
    NumTrials = length(files2get);
    for z = 1:NumTrials
        FullFileName = char(files2get{z});
        GaitMeasures(z) = GAMSmini_GaitMeasuresCalc(FullFileName, age);
        GM(z).c3d_name = FullFileName;
    end
else
    NumTrials = 1;
    FullFileName = files2get;
end

%% Analyze C3D file(s)
if NumTrials > 1
    for z = 1:NumTrials
        % structure of each trial's gait measures
        %GM(z).OGP_vars = GaitMeasures(z).OGP_vars;
%         GM(z).Spatial_vars = GaitMeasures(z).Spatial_vars;
        GM(z).Temporal_vars = GaitMeasures(z).Temporal_vars;
%         GM(z).TempSpat_table_data = GaitMeasures(z).TempSpat_table_data;
%         GM(z).GaitParams = GaitMeasures(z).GaitParams;
        GM(z).TempSpat = GaitMeasures(z).TempSpat;
    end
    
else % only 1 trial
    GaitMeasures = GAMSmini_GaitMeasuresCalc(FullFileName, age);
    % builds structure of current trial's name
    GM.c3d_name = FullFileName;
    % structure of each trial's gait measures
    %GM.OGP_vars = GaitMeasures.OGP_vars;
%     GM.Spatial_vars = GaitMeasures.Spatial_vars;
    GM.Temporal_vars = GaitMeasures.Temporal_vars;
%     GM.TempSpat_table_data = GaitMeasures.TempSpat_table_data;
%     GM.GaitParams = GaitMeasures.GaitParams;
    GM.TempSpat = GaitMeasures.TempSpat;
    
end

end





