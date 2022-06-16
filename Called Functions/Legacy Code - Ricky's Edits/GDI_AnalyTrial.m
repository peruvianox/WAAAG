function [variance_ratios, residuals, ave_kinematics, stddev_kinematics, kinematics, TempSpat, GaitEvents] = GDI_AnalyTrial(FullFileName)
% % -------------------------------------------------------------------------
% ANALYZES TRIAL FOR GDI PROGRAM
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: July 14, 2009
% Updated by Colton Sauer July 2014
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Adjusted by Colton Sauer on July 7, 2014 to output left and right gait
% vectors separately.
%
% Description:	This function calculates the variance ratio for the 9
%               kinematic curves of the specified C3D file.
%               Kinematic curves used: Pelvis (all 3), Hip (all 3), Knee
%               (flexion/extension), Ankle (dorsi/plantar), FootProgression.
%
% Input:        FullFileName        user defined C3D file to be analyzed
%               
%
% Output:       variance_ratios         matrix of variance ratios for all 9
%                                       kinematic curves
%               residuals               residual (difference from mean) for
%                                       each gait cycle
%               ave kinematics          matrix of 9 mean kinematic curves for entire
%                                       trial, sampled at 2% of gait cycle
%               stddev_kinematics       matrix of standard deviation for
%                                       kinematics curves

% DEBUG LOOP
% close all
% clear all
% clc
% 
% [FullFileName, pathname] = uigetfile('*.c3d*', 'Please select a C3D file to be analyzed');

%%% END DEBUG LOOP %%%%%

%% Opens and Extract C3D file's data
[C3Dfile] = OpenC3D(FullFileName);

% Open and Read the user specified C3D file
% and returns: HeaderInfo, Point_xyz, PointLabels, and GaitEvents
[HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);


%% Organizes Kinematic Angles
% sorts kinematic angles into a structure
MkrSet_type = 'lowerbody';  alt_str = '';
[sortedPoints] = sortPoints(Point_xyz, PointLabels, MkrSet_type, alt_str);

% Currently, not including shank angles
Variables = sortedPoints(23:32);

%% Define Kinetics
Kinetics = sortedPoints(35:52);

%% Temporal Spatial Variables
% [TempSpat] = TemporalSpatial(HeaderInfo, Point_xyz, PointLabels, GaitEvents);
[TempSpat] = TemporalSpatial(HeaderInfo, Point_xyz, PointLabels, GaitEvents);

%% GAIT CYCLE
% Find mean ensemble of each gait cycle for the c3d file
FirstFrame = HeaderInfo(4).data;
LastFrame = HeaderInfo(5).data;

% Calculate times for left & right gait events and creates matrix of each
% side's gait cycles
[Events_times, LGaitCycles, RGaitCycles] = GaitEventTimes(GaitEvents);
N_GaitCycles = length(LGaitCycles(1,1,:)) + length(RGaitCycles(1,1,:));

%% Separate kinematic angles by gait cycle

% Sample Angles at 1% Increments of Gait Cycle (101pts)
new_samp_size = 101;

% Note: Left (odd counter) variables are listed before Right (even counter)
% variables
for n = 1:size(Variables,2)  % number of variables to process
    
    if (mod(n,2) == 1) 
        % LEFT SIDE
        % if current counter (n) is an odd number, then corresponds to left
        % side data      
        for p = 1:3  % for each plane (sag, coron, transv)
            % If variable exists then process, if not skip
            if (size(Variables(n).data,3) > 1)  % if data exhists in current variable
                % name of current variable
                Nvariable(n).VarName = Variables(n).name;
                
                % Parse current Variable (ie. kinematic, kinetic) into % gait cycle (0-100%)
                variable = Variables(n).data(p,:,:);  % Extract current variable's sagittal, coronal, & rotational angles
                GaitCycle_times = LGaitCycles;  % designate side of gait cycle times
                [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);  % calculate as % gait cycle
                Nvariable(n).VarPlane_GCs(p).left = var_GCs;  % parsed as % gait cycle for variable(n), for current plane (p)
               
                variable = Nvariable(n).VarPlane_GCs(p).left;  % parsed gait cycles for variable(n), for each plane p
                [mean_var, std_dev] = MeanEnsemble(variable);
                Nvariable(n).VarPlane_mean(p).left = mean_var;  % mean of variable(n) in plane 'p'
                Nvariable(n).VarPlane_stddev(p).left = std_dev;  % std dev of variable(n) in plane 'p'
                
                [resid] = Residual(mean_var, variable);  % calculate residual
                Nvariable(n).VarPlane_resid(p).left = resid;

                % Calculate gait cycle variance ratios
                [VR, VR_GC] = VarianceRatio(variable, mean_var);
                Nvariable(n).VarPlane_varratio(p).left = VR;  % variance ratio for this variable
                Nvariable(n).VarPlane_varratioGC(p).left = VR_GC;  % variance ratio for each gait cycle
            end
        end
        
    elseif (mod(n,2) == 0) 
        % RIGHT SIDE
        % if current counter (n) is an even number, then corresponds to
        % right side data 
        for p = 1:3  % for each plane (sag, coron, transv)
            % If variable exists then process, if not skip
            if (size(Variables(n).data,3) > 1)  % if data exhists in current variable

                % name of current variable
                Nvariable(n).VarName = Variables(n).name;

                % Parse current Variable (ie. kinematic, kinetic) into % gait cycle (0-100%)
                variable = Variables(n).data(p,:,:);
                GaitCycle_times = RGaitCycles;
                [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);
                Nvariable(n).VarPlane_GCs(p).right = var_GCs;
                
                %Find mean and standard deviation
                variable = Nvariable(n).VarPlane_GCs(p).right;  % parsed gait cycles for variable(n), for each plane p
                [mean_var, std_dev] = MeanEnsemble(variable);
                Nvariable(n).VarPlane_mean(p).right = mean_var;  % mean of variable(n) in plane 'p'
                Nvariable(n).VarPlane_stddev(p).right = std_dev;  % std dev of variable(n) in plane 'p'
                
                [resid] = Residual(mean_var, variable);  % calculate residual
                Nvariable(n).VarPlane_resid(p).right = resid;

                % Calculate gait cycle variance ratios
                [VR, VR_GC] = VarianceRatio(variable, mean_var);
                Nvariable(n).VarPlane_varratio(p).right = VR;  % variance ratio for this variable
                Nvariable(n).VarPlane_varratioGC(p).right = VR_GC;  % variance ratio for each gait cycle
            end
        end
    end
end
% Composite matrix of all variance ratios for the C3D file
L_variance_ratios_trial = [Nvariable(1).VarPlane_varratio(1).left; Nvariable(1).VarPlane_varratio(2).left; Nvariable(1).VarPlane_varratio(3).left;
                         Nvariable(3).VarPlane_varratio(1).left; Nvariable(3).VarPlane_varratio(2).left; Nvariable(3).VarPlane_varratio(3).left;
                         Nvariable(5).VarPlane_varratio(1).left; Nvariable(7).VarPlane_varratio(1).left; Nvariable(9).VarPlane_varratio(3).left];
R_variance_ratios_trial = [Nvariable(2).VarPlane_varratio(1).right; Nvariable(2).VarPlane_varratio(2).right; Nvariable(2).VarPlane_varratio(3).right;
                         Nvariable(4).VarPlane_varratio(1).right; Nvariable(4).VarPlane_varratio(2).right; Nvariable(4).VarPlane_varratio(3).right;
                         Nvariable(6).VarPlane_varratio(1).right; Nvariable(8).VarPlane_varratio(1).right; Nvariable(10).VarPlane_varratio(3).right];                     

% summed value of all curves' variance ratio, providing one variance ratio for trial
variance_ratios.trial.left = sum(L_variance_ratios_trial);
variance_ratios.trial.right = sum(R_variance_ratios_trial);

% Variance ratio for each gait cycle
L_variance_ratios_gaitycle = [Nvariable(1).VarPlane_varratioGC(1).left; Nvariable(1).VarPlane_varratioGC(2).left; Nvariable(1).VarPlane_varratioGC(3).left;
                         Nvariable(3).VarPlane_varratioGC(1).left; Nvariable(3).VarPlane_varratioGC(2).left; Nvariable(3).VarPlane_varratioGC(3).left;
                         Nvariable(5).VarPlane_varratioGC(1).left; Nvariable(7).VarPlane_varratioGC(1).left; Nvariable(9).VarPlane_varratioGC(3).left];
R_variance_ratios_gaitcycle = [Nvariable(2).VarPlane_varratioGC(1).right; Nvariable(2).VarPlane_varratioGC(2).right; Nvariable(2).VarPlane_varratioGC(3).right;
                         Nvariable(4).VarPlane_varratioGC(1).right; Nvariable(4).VarPlane_varratioGC(2).right; Nvariable(4).VarPlane_varratioGC(3).right;
                         Nvariable(6).VarPlane_varratioGC(1).right; Nvariable(8).VarPlane_varratioGC(1).right; Nvariable(10).VarPlane_varratioGC(3).right];


% Residuals for all gait cycles (column)
L_residuals_all = [Nvariable(1).VarPlane_resid(1).left; Nvariable(1).VarPlane_resid(2).left; Nvariable(1).VarPlane_resid(3).left;
                         Nvariable(3).VarPlane_resid(1).left; Nvariable(3).VarPlane_resid(2).left; Nvariable(3).VarPlane_resid(3).left;
                         Nvariable(5).VarPlane_resid(1).left; Nvariable(7).VarPlane_resid(1).left; Nvariable(9).VarPlane_resid(3).left];
R_residuals_all = [Nvariable(2).VarPlane_resid(1).right; Nvariable(2).VarPlane_resid(2).right; Nvariable(2).VarPlane_resid(3).right;
                         Nvariable(4).VarPlane_resid(1).right; Nvariable(4).VarPlane_resid(2).right; Nvariable(4).VarPlane_resid(3).right;
                         Nvariable(6).VarPlane_resid(1).right; Nvariable(8).VarPlane_resid(1).right; Nvariable(10).VarPlane_resid(3).right];
                     
% rank of all gait cycles using residuals
for t = 1:length(L_residuals_all(1,:))
    % sum of residuals for all gait cycles
    L_sum_residuals_all(:,t) = sum(L_residuals_all(:,t));
end
for t = 1:length(R_residuals_all(1,:))
    % sum of residuals for all gait cycles
    R_sum_residuals_all(:,t) = sum(R_residuals_all(:,t));
end

% sorts residuals
[L_resid_val, L_gaitcycle] = sort(L_sum_residuals_all, 'ascend');
[R_resid_val, R_gaitcycle] = sort(R_sum_residuals_all, 'ascend');
% provides the representative gait cycle, in ascending order
residuals.all_GCs.left = L_gaitcycle;
residuals.all_GCs.right = R_gaitcycle;

for s = 1:length(LGaitCycles(1,1,:))
    % for each left gait cycle, add string 'left' to matrix
    residuals.GC_sides(s).name = 'left';
end

for r = (length(LGaitCycles(1,1,:))+1):N_GaitCycles
    % add string 'right' for remaining gait cycles
    residuals.GC_sides(r).name = 'right';
end


%% Output Variables
% Composite structure of 9 kinematic curves for C3D file, sampled at 1%

kinematics(1).name = 'PelvisT';
kinematics(1).data.left = Nvariable(1).VarPlane_GCs(1).left;
kinematics(1).data.right = Nvariable(2).VarPlane_GCs(1).right; 
kinematics(2).name = 'PelvisO';
kinematics(2).data.left = Nvariable(1).VarPlane_GCs(2).left;
kinematics(2).data.right = Nvariable(2).VarPlane_GCs(2).right;
kinematics(3).name = 'PelviskiR';
kinematics(3).data.left = Nvariable(1).VarPlane_GCs(3).left;
kinematics(3).data.right = Nvariable(2).VarPlane_GCs(3).right;
kinematics(4).name = 'HipFE';
kinematics(4).data.left = Nvariable(3).VarPlane_GCs(1).left;
kinematics(4).data.right = Nvariable(4).VarPlane_GCs(1).right;
kinematics(5).name = 'HipAA';
kinematics(5).data.left = Nvariable(3).VarPlane_GCs(2).left;
kinematics(5).data.right = Nvariable(4).VarPlane_GCs(2).right;
kinematics(6).name = 'HipR';
kinematics(6).data.left = Nvariable(3).VarPlane_GCs(3).left;
kinematics(6).data.right = Nvariable(4).VarPlane_GCs(3).right;
kinematics(7).name = 'KneeFE';
kinematics(7).data.left = Nvariable(5).VarPlane_GCs(1).left;
kinematics(7).data.right = Nvariable(6).VarPlane_GCs(1).right;
kinematics(8).name = 'AnkleDP';
kinematics(8).data.left = Nvariable(7).VarPlane_GCs(1).left;
kinematics(8).data.right = Nvariable(8).VarPlane_GCs(1).right;
kinematics(9).name = 'FootP';
kinematics(9).data.left = Nvariable(9).VarPlane_GCs(3).left;
kinematics(9).data.right = Nvariable(10).VarPlane_GCs(3).right;

kinematics(10).name = 'KneeABAD';
kinematics(10).data.left = Nvariable(5).VarPlane_GCs(2).left;
kinematics(10).data.right = Nvariable(6).VarPlane_GCs(2).right;
kinematics(11).name = 'KneeR';
kinematics(11).data.left = Nvariable(5).VarPlane_GCs(3).left;
kinematics(11).data.right = Nvariable(6).VarPlane_GCs(3).right;
kinematics(12).name = 'FootR';
kinematics(12).data.left = Nvariable(7).VarPlane_GCs(3).left;
kinematics(12).data.right = Nvariable(8).VarPlane_GCs(3).right;

% Global Foot Angles other than Foot Prog 
kinematics(13).name = 'FootPitch';
kinematics(13).data.left = Nvariable(9).VarPlane_GCs(1).left;
kinematics(13).data.right = Nvariable(10).VarPlane_GCs(1).right;
kinematics(14).name = 'FootRoll';
kinematics(14).data.left = Nvariable(9).VarPlane_GCs(2).left;
kinematics(14).data.right = Nvariable(10).VarPlane_GCs(2).right;

% Composite structure of mean of 9 kinematic curves for C3D file, sampled at 1%
ave_kinematics(1).name = 'ave_PelvisT';
ave_kinematics(1).data.left = Nvariable(1).VarPlane_mean(1).left;
ave_kinematics(1).data.right = Nvariable(2).VarPlane_mean(1).right;
ave_kinematics(2).name = 'ave_PelvisO';
ave_kinematics(2).data.left = Nvariable(1).VarPlane_mean(2).left;
ave_kinematics(2).data.right = Nvariable(2).VarPlane_mean(2).right;
ave_kinematics(3).name = 'ave_PelvisR';
ave_kinematics(3).data.left = Nvariable(1).VarPlane_mean(3).left;
ave_kinematics(3).data.right = Nvariable(2).VarPlane_mean(3).right;
ave_kinematics(4).name = 'ave_HipFE';
ave_kinematics(4).data.left = Nvariable(3).VarPlane_mean(1).left;
ave_kinematics(4).data.right = Nvariable(4).VarPlane_mean(1).right;
ave_kinematics(5).name = 'ave_HipAA';
ave_kinematics(5).data.left = Nvariable(3).VarPlane_mean(2).left;
ave_kinematics(5).data.right = Nvariable(4).VarPlane_mean(2).right;
ave_kinematics(6).name = 'ave_HipR';
ave_kinematics(6).data.left = Nvariable(3).VarPlane_mean(3).left;
ave_kinematics(6).data.right = Nvariable(4).VarPlane_mean(3).right;
ave_kinematics(7).name = 'ave_KneeFE';
ave_kinematics(7).data.left = Nvariable(5).VarPlane_mean(1).left;
ave_kinematics(7).data.right = Nvariable(6).VarPlane_mean(1).right;
ave_kinematics(8).name = 'ave_AnkleDP';
ave_kinematics(8).data.left = Nvariable(7).VarPlane_mean(1).left; 
ave_kinematics(8).data.right = Nvariable(8).VarPlane_mean(1).right; 
ave_kinematics(9).name = 'ave_FootP';
ave_kinematics(9).data.left = Nvariable(9).VarPlane_mean(3).left;
ave_kinematics(9).data.right = Nvariable(10).VarPlane_mean(3).right;

ave_kinematics(10).name = 'KneeABAD';
ave_kinematics(10).data.left = Nvariable(5).VarPlane_mean(2).left;
ave_kinematics(10).data.right = Nvariable(6).VarPlane_mean(2).right;
ave_kinematics(11).name = 'KneeR';
ave_kinematics(11).data.left = Nvariable(5).VarPlane_mean(3).left;
ave_kinematics(11).data.right = Nvariable(6).VarPlane_mean(3).right;
ave_kinematics(12).name = 'FootR';
ave_kinematics(12).data.left = Nvariable(7).VarPlane_mean(3).left;
ave_kinematics(12).data.right = Nvariable(8).VarPlane_mean(3).right;

% Composite structure of standard deviation of 9 kinematic curves for C3D file, sampled at 1%
stddev_kinematics(1).name = 'stddev_PelvisT';
stddev_kinematics(1).data.left = Nvariable(1).VarPlane_stddev(1).left;
stddev_kinematics(1).data.right = Nvariable(2).VarPlane_stddev(1).right;
stddev_kinematics(2).name = 'stddev_PelvisO';
stddev_kinematics(2).data.left = Nvariable(1).VarPlane_stddev(2).left; 
stddev_kinematics(2).data.right = Nvariable(2).VarPlane_stddev(2).right; 
stddev_kinematics(3).name = 'stddev_PelvisR';
stddev_kinematics(3).data.left = Nvariable(1).VarPlane_stddev(3).left;
stddev_kinematics(3).data.right = Nvariable(2).VarPlane_stddev(3).right;
stddev_kinematics(4).name = 'stddev_HipFE';
stddev_kinematics(4).data.left = Nvariable(3).VarPlane_stddev(1).left;
stddev_kinematics(4).data.right = Nvariable(4).VarPlane_stddev(1).right;
stddev_kinematics(5).name = 'stddev_HipAA';
stddev_kinematics(5).data.left = Nvariable(3).VarPlane_stddev(2).left;
stddev_kinematics(5).data.right = Nvariable(4).VarPlane_stddev(2).right;
stddev_kinematics(6).name = 'stddev_HipR';
stddev_kinematics(6).data.left = Nvariable(3).VarPlane_stddev(3).left;
stddev_kinematics(6).data.right = Nvariable(4).VarPlane_stddev(3).right;
stddev_kinematics(7).name = 'stddev_KneeFE';
stddev_kinematics(7).data.left = Nvariable(5).VarPlane_stddev(1).left;
stddev_kinematics(7).data.right = Nvariable(6).VarPlane_stddev(1).right;
stddev_kinematics(8).name = 'stddev_AnkleDP';
stddev_kinematics(8).data.left = Nvariable(7).VarPlane_stddev(1).left;
stddev_kinematics(8).data.right = Nvariable(8).VarPlane_stddev(1).right; 
stddev_kinematics(9).name = 'stddev_FootP';
stddev_kinematics(9).data.left = Nvariable(9).VarPlane_stddev(3).left;
stddev_kinematics(9).data.right = Nvariable(10).VarPlane_stddev(3).right;

stddev_kinematics(10).name = 'KneeABAD';
stddev_kinematics(10).data.left = Nvariable(5).VarPlane_stddev(2).left;
stddev_kinematics(10).data.right = Nvariable(6).VarPlane_stddev(2).right;
stddev_kinematics(11).name = 'KneeR';
stddev_kinematics(11).data.left = Nvariable(5).VarPlane_stddev(3).left;
stddev_kinematics(11).data.right = Nvariable(6).VarPlane_stddev(3).right;
stddev_kinematics(12).name = 'FootR';
stddev_kinematics(12).data.left = Nvariable(7).VarPlane_stddev(3).left;
stddev_kinematics(12).data.right = Nvariable(8).VarPlane_stddev(3).right;

end
              