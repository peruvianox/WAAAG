function [Nvariable, ave_rotations, TempSpat] = GAMS_AnalyTrial(FullFileName)
% -------------------------------------------------------------------------
% ANALYZES TRIAL FOR GAMS PROGRAM
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: October 1, 2010
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This function calculates the variance ratio for all
%               kinematic curves of the specified C3D file. Updated to
%               calculate VRs more accurately (side independent).
%               
%
% Input:        FullFileName        user defined C3D file to be analyzed
%               
%
% Output:       variance_ratios         matrix of variance ratios for all 9
%                                       kinematic curves
%               residuals               residual (difference from mean) for each gait cycle
%               ave_kinematics          structure of all kinematic curves for
%                                       trial
%               TempSpat                structure of all temporal spatial
%                                       variables
%               


%% DEBUG LOOP
% close all
% clear all
% clc
% 
% side = 'bilat';
% 
% [FullFileName, pathname] = uigetfile('*.c3d', 'Please select a C3D file to be analyzed');

% %%%% END DEBUG LOOP %%%%%


%% Opens and Extract C3D file's data
[C3Dfile] = OpenC3D(FullFileName);

% Open and Read the user specified C3D file
% and returns: HeaderInfo, Point_xyz, PointLabels
[HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);


%% Organizes Kinematic Angles
% sorts kinematic angles and kinetics into a structure
[Kinematics, Kinetics] = sortPointLabels(Point_xyz, PointLabels);

% Lower body kinematics only
Kinematics_m = Kinematics(1:14);
Variables = Kinematics_m;


%% GAIT CYCLE
% Find mean ensemble of each gait cycle for the c3d file
FirstFrame = HeaderInfo(4).data;
LastFrame = HeaderInfo(5).data;

% Calculate times for left & right gait events and creates matrix of each
% side's gait cycles
[Events_times, LGaitCycles, RGaitCycles] = GaitEventTimes(GaitEvents);
LFS_times = Events_times(1).data;
LFO_times = Events_times(2).data;
RFS_times = Events_times(4).data;
RFO_times = Events_times(5).data;

%% Process (All but Powers)
% Sample Angles at 1% Increments of Gait Cycle (101pts)
new_samp_size = 101;

for n = 1:size(Variables,2)  % number of variables to process
    for p = 1:3  % for each plane (sag, coron, transv)
        % If variable exists then process, if not skip
        if (size(Variables(n).data,3) > 1)  % if data exhists in current variable

            % name of current (newly processed) variable
            Nvariable(n).VarName = Variables(n).name;
            
            % name of current (newly processed) variable
            Nvariable(n).VarName = Variables(n).name;

            % LEFT SIDE
            % Parse current Variable (ie. kinematic, kinetic) into % gait cycle (0-100%)
            variable = Variables(n).data(p,:,:);  % Extract current variable's sagittal, coronal, & rotational angles
            GaitCycle_times = LGaitCycles;  % designate side of gait cycle times
            [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);  % calculate as % gait cycle
            VarPlane_GC(p).left = var_GCs;  % parsed as % gait cycle for variable(n), for current plane (p)
            Nvariable(n).VarPlane_GCs = VarPlane_GC;

            % RIGHT SIDE
            % Parse current Variable (ie. kinematic, kinetic) into % gait cycle (0-100%)
            GaitCycle_times = RGaitCycles;
            [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);
            VarPlane_GC(p).right = var_GCs;  % parsed gait cycles for variable n
            Nvariable(n).VarPlane_GCs = VarPlane_GC;

            % LEFT SIDE
            % Calculate mean and standard deviation
            variable = Nvariable(n).VarPlane_GCs(p).left;  % parsed gait cycles for variable(n), for each plane p
            [mean_var, std_dev] = MeanEnsemble(variable);
            VarPlane_mean(p).left = mean_var;  % mean of variable(n) in plane 'p'
            Nvariable(n).VarPlane_mean = VarPlane_mean;  
            VarPlane_stddev(p).left = std_dev;  % std dev of variable(n) in plane 'p'
            Nvariable(n).VarPlane_stddev = VarPlane_stddev;  
            [resid] = Residual(mean_var, variable);  % calculate residual
            VarPlane_resid(p).left = resid;
            Nvariable(n).VarPlane_resid = VarPlane_resid;
            % Calculate gait cycle variance ratios
            [VR, VR_GC] = VarianceRatio(variable, mean_var);
            VarPlane_varratio(p).left = VR; 
            Nvariable(n).VarPlane_varratio = VarPlane_varratio;  % variance ratio for this variable
            VarPlane_varratioGC(p).left = VR_GC;
            Nvariable(n).VarPlane_varratioGC = VarPlane_varratioGC;  % variance ratio for each gait cycle

            % RIGHT SIDE
            % Calculate mean and standard deviation
            variable = Nvariable(n).VarPlane_GCs(p).right;  % parsed gait cycles for variable(n), for each plane p
            [mean_var, std_dev] = MeanEnsemble(variable);
            VarPlane_mean(p).right = mean_var;  % mean of variable(n) in plane 'p'
            Nvariable(n).VarPlane_mean = VarPlane_mean;  
            VarPlane_stddev(p).right = std_dev;  % std dev of variable(n) in plane 'p'
            Nvariable(n).VarPlane_stddev = VarPlane_stddev;  
            [resid] = Residual(mean_var, variable);  % calculate residual
            VarPlane_resid(p).right = resid;
            Nvariable(n).VarPlane_resid = VarPlane_resid;
            % Calculate gait cycle variance ratios
            [VR, VR_GC] = VarianceRatio(variable, mean_var);
            VarPlane_varratio(p).right = VR; 
            Nvariable(n).VarPlane_varratio = VarPlane_varratio;
            VarPlane_varratioGC(p).right = VR_GC;
            Nvariable(n).VarPlane_varratioGC = VarPlane_varratioGC;
        end
    end
end

                     
%% AVERAGE KINEMATIC ROTATION ANGLES
% Calculate averages of transverse plane kinematic angles

ave_rotations(1).name = 'ave_LPelvRot';
ave_rotations(1).data = mean(Nvariable(1).VarPlane_mean(3).left);
ave_rotations(2).name = 'ave_RPelvRot';
ave_rotations(2).data = mean(Nvariable(2).VarPlane_mean(3).right);
ave_rotations(3).name = 'ave_LHipRot';
ave_rotations(3).data = mean(Nvariable(3).VarPlane_mean(3).left);
ave_rotations(4).name = 'ave_RHipRot';
ave_rotations(4).data = mean(Nvariable(4).VarPlane_mean(3).right);
ave_rotations(5).name = 'ave_LKneeRot';
ave_rotations(5).data = mean(Nvariable(5).VarPlane_mean(3).left);
ave_rotations(6).name = 'ave_RKneeRot';
ave_rotations(6).data = mean(Nvariable(6).VarPlane_mean(3).right);
ave_rotations(7).name = 'ave_LAnkRot';
ave_rotations(7).data = mean(Nvariable(7).VarPlane_mean(3).left);
ave_rotations(8).name = 'ave_RAnkRot';
ave_rotations(8).data = mean(Nvariable(8).VarPlane_mean(3).right);
ave_rotations(9).name = 'ave_LFootProg';
ave_rotations(10).name = 'ave_RFootProg';
%Setup all Footprog data, note this is not time normalized
LFProgAll=squeeze(Variables(9).data(3,:,:));
RFProgAll=squeeze(Variables(10).data(3,:,:));
%readjust gait event times for start tof trial
LFS_timesAdj=LFS_times-HeaderInfo(4).data+1;
LFO_timesAdj=LFO_times-HeaderInfo(4).data+1;
RFS_timesAdj=RFS_times-HeaderInfo(4).data+1;
RFO_timesAdj=RFO_times-HeaderInfo(4).data+1;
%Average only stance phase
for n=1:length(LFS_timesAdj)
    m=find(LFO_timesAdj>LFS_timesAdj(n),1); % find the next index of foot off after the foot strike
        if isempty(m)==1 %if there is no foot off after the foot strike
            break
        end
    LFPAvgAll(n)=mean(LFProgAll(LFS_timesAdj(n):LFO_timesAdj(m)));    
end
for n=1:length(RFS_timesAdj)
    m=find(RFO_timesAdj>RFS_timesAdj(n),1); % find the next index of foot off after the foot strike
        if isempty(m)==1 %if there is no foot off after the foot strike
            break
        end
    RFPAvgAll(n)=mean(RFProgAll(RFS_timesAdj(n):RFO_timesAdj(m)));    
end
ave_rotations(9).data=mean(LFPAvgAll);
ave_rotations(10).data=mean(RFPAvgAll);

if (length(Kinematics_m(11).data) > 3)
    % For Gaia models only
    ave_rotations(11).name = 'ave_LShankRot';
    ave_rotations(11).data = mean(Nvariable(11).VarPlane_mean(3).left);
end

if (length(Kinematics_m(12).data) > 3)
    % For Gaia models only
    ave_rotations(12).name = 'ave_RShankRot';
    ave_rotations(12).data = mean(Nvariable(12).VarPlane_mean(3).right);
end


%% Temporal Spatial Variables for Trial
[TempSpat] = TemporalSpatial(HeaderInfo, Point_xyz, PointLabels, GaitEvents);


%% Variance Ratios & Residuals
% Composite matrix of all variance ratios and residuals for C3D file's lower body kinematics

% Variance ratio & residuals for each left gait cycle, lower body variables
for n = 1:length(LGaitCycles(1,1,:))
    LGC(n).varratio = [Nvariable(1).VarPlane_varratioGC(1).left(n); Nvariable(1).VarPlane_varratioGC(2).left(n); Nvariable(1).VarPlane_varratioGC(3).left(n);
                       Nvariable(3).VarPlane_varratioGC(1).left(n); Nvariable(3).VarPlane_varratioGC(2).left(n); Nvariable(3).VarPlane_varratioGC(3).left(n);
                       Nvariable(5).VarPlane_varratioGC(1).left(n); Nvariable(5).VarPlane_varratioGC(2).left(n); Nvariable(5).VarPlane_varratioGC(3).left(n);
                       Nvariable(7).VarPlane_varratioGC(1).left(n); Nvariable(9).VarPlane_varratioGC(3).left(n); Nvariable(7).VarPlane_varratioGC(3).left(n)];

    LGC(n).resid = [Nvariable(1).VarPlane_resid(1).left(n); Nvariable(1).VarPlane_resid(2).left(n); Nvariable(1).VarPlane_resid(3).left(n);
                    Nvariable(3).VarPlane_resid(1).left(n); Nvariable(3).VarPlane_resid(2).left(n); Nvariable(3).VarPlane_resid(3).left(n);
                    Nvariable(5).VarPlane_resid(1).left(n); Nvariable(5).VarPlane_resid(2).left(n); Nvariable(5).VarPlane_resid(3).left(n);
                    Nvariable(7).VarPlane_resid(1).left(n); Nvariable(9).VarPlane_resid(3).left(n); Nvariable(7).VarPlane_resid(3).left(n)];
end

% Variance ratio & residuals for each right gait cycle, lower body variables
for r = 1:length(RGaitCycles(1,1,:))
    RGC(r).varratio = [Nvariable(2).VarPlane_varratioGC(1).right(r); Nvariable(2).VarPlane_varratioGC(2).right(r); Nvariable(2).VarPlane_varratioGC(3).right(r);
                       Nvariable(4).VarPlane_varratioGC(1).right(r); Nvariable(4).VarPlane_varratioGC(2).right(r); Nvariable(4).VarPlane_varratioGC(3).right(r);
                       Nvariable(6).VarPlane_varratioGC(1).right(r); Nvariable(6).VarPlane_varratioGC(2).right(r); Nvariable(6).VarPlane_varratioGC(3).right(r);
                       Nvariable(8).VarPlane_varratioGC(1).right(r); Nvariable(10).VarPlane_varratioGC(3).right(r); Nvariable(8).VarPlane_varratioGC(3).right(r)];

    RGC(r).resid = [Nvariable(2).VarPlane_resid(1).right(r); Nvariable(2).VarPlane_resid(2).right(r); Nvariable(2).VarPlane_resid(3).right(r);
                    Nvariable(4).VarPlane_resid(1).right(r); Nvariable(4).VarPlane_resid(2).right(r); Nvariable(4).VarPlane_resid(3).right(r);
                    Nvariable(6).VarPlane_resid(1).right(r); Nvariable(6).VarPlane_resid(2).right(r); Nvariable(6).VarPlane_resid(3).right(r);
                    Nvariable(8).VarPlane_resid(1).right(r); Nvariable(10).VarPlane_resid(3).right(r); Nvariable(8).VarPlane_resid(3).right(r)];
end


% Sum all left gait cycle variance ratios & residuals
for j = 1:length(LGaitCycles(1,1,:))
    LGCs_VR_sums(j) = sum(LGC(j).varratio);
    LGCs_resid_sums(j) = sum(LGC(j).resid);
end
VR_trial_left = sum(LGCs_VR_sums);

% Sum all right gait cycle variance ratios & residuals
for w = 1:length(RGaitCycles(1,1,:))
    RGCs_VR_sums(w) = sum(RGC(w).varratio);
    RGCs_resid_sums(w) = sum(RGC(w).resid);
end
VR_trial_right = sum(RGCs_VR_sums);

% Sum of left and right side variance ratios
variance_ratio_trial = VR_trial_left + VR_trial_right;

% TRIAL'S TOTAL, BILATERAL VARIANCE RATIO
Nvariable(size(Variables,2)+1).variance_ratio_trial = variance_ratio_trial;

% Sort left side's residuals
[resid_val, gaitcycle] = sort(LGCs_resid_sums, 'ascend');
% provides the representative gait cycle, in ascending order
residuals.left = gaitcycle;

% Sort right side's residuals
[resid_val, gaitcycle] = sort(RGCs_resid_sums, 'ascend');
% provides the representative gait cycle, in ascending order
residuals.right = gaitcycle;

% TRIAL'S RESIDUALS FOR EACH SIDE'S GAIT CYCLES
Nvariable(size(Variables,2)+1).GCresiduals = residuals;
end



