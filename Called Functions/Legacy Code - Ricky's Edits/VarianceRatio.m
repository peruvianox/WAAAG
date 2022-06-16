function [VR, VR_GC] = VarianceRatio(variable, mean_var)

% -------------------------------------------------------------------------
% VARIANCE RATIO (aka F Distribution)
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: July 20, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This is function calculates the variance ratio 
%               of 'variable'.
%
% Input:        variable                 variable argument input to be parsed
%               mean_var                 mean of input variable
%
% OUTPUTS:  
%               VR                       variance ratio of each gait cycle
%                                        for all temporal points
%
% AKNOWLEDGEMENTS:
%               Variance ratio calculated by method described in "Repeatability of Phasic
%               Muscle Activity: Performance of Surface and Intramuscular Wire Electrodes
%               in Gait Analysis" by M.P. Kadaba, et al.
%               NOTE: all variables from this paper are provided in ().

%% DEBUG LOOP
% close all
% clear all
% clc

% GaitCycle_num = 0;
% side = 'Bilateral';
% 
% [FullFileName, pathname] = uigetfile('*.c3d*', 'Please select C3D file to be analyzed');
% 
% 
% % Opens and Extract C3D file's data
% [C3Dfile] = OpenC3D(FullFileName);
% 
% % Open and Read the user specified C3D file
% % and returns: HeaderInfo, Point_xyz, PointLabels
% [HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);
% 
% 
% % Temporal Spatial Variables
% [TempSpat] = TemporalSpatial(HeaderInfo, Point_xyz, PointLabels, GaitEvents);
% 
% 
% % Organizes Kinematic Angles
% % sorts kinematic angles into a structure
% [Kinematics] = sortKinematics(Point_xyz, PointLabels);
% 
% LPelvisAngles = Kinematics(1).data;
% RPelvisAngles = Kinematics(2).data;
% LHipAngles = Kinematics(3).data;
% RHipAngles = Kinematics(4).data;
% LKneeAngles = Kinematics(5).data;
% RKneeAngles = Kinematics(6).data;
% LShankAngles = Kinematics(11).data;
% RShankAngles = Kinematics(12).data;
% LAnkleAngles = Kinematics(7).data;
% RAnkleAngles = Kinematics(8).data;
% LFootProgressAngles = Kinematics(9).data;
% RFootProgressAngles = Kinematics(10).data;

%%%%% END DEBUG LOOP %%%%%

%% GAIT CYCLE
% find mean ensemble of each gait cycle for the c3d file
% find which cycle deviates the least for that entire trial(s)

% FirstFrame = HeaderInfo(4).data;
% LastFrame = HeaderInfo(5).data;
% 
% % Calculate times for left & right gait events and creates matrix of each
% % side's gait cycles
% [LFS_times, LFO_times, RFS_times, RFO_times, LGaitCycles, RGaitCycles] = GaitEventTimes(GaitEvents);
% new_samp_size = 101;


% % Only left
% variable = LPelvisAngles(1,:,:);
% GaitCycle_times = LGaitCycles;
% [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);
% LPelvTilt_GC = var_GCs;
% 
% variable = LPelvTilt_GC;
% [mean_var, std_dev] = MeanEnsemble(variable); % mean ensemble function


% % Only right
% variable = RPelvisAngles(1,:,:);
% GaitCycle_times = RGaitCycles;
% [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);
% RPelvTilt_GC = var_GCs;
% 
% variable = RPelvTilt_GC;
% [mean_var, std_dev] = MeanEnsemble(variable); % mean ensemble function


% % Bilateral
% LPelvTilt_GC(:,:,:) = LPelvTilt_GC((1:length(LPelvTilt_GC)),:);
% RPelvTilt_GC(:,:,:) = RPelvTilt_GC((1:length(RPelvTilt_GC)),:);
% PelvTilt_GC(:,:,:) = [LPelvTilt_GC((1:length(LPelvTilt_GC)),:) RPelvTilt_GC((1:length(RPelvTilt_GC)),:)];
% 
% variable = PelvTilt_GC;
% [mean_var, std_dev] = MeanEnsemble(variable); % mean ensemble function
% mean_PelvTilt = mean_var; % mean for all gait cycles
% stddev_PelvTilt = std_dev; % std dev of each gait cycle



%% Variance Ratio

% n is the number of gait cycles
n = length(variable(1,:));

% m is the total number of temporal points (100 points, each a % of gait
% cycle)
m = length(variable(:,1))-1;

% (Eij) value of the jth data point of 'variable' at time epoch t1

for t = 1:m
    % (Ei_bar) averaged values of 'variable' at time epoch t, averaged over all gait cycles
    Ei_bar(t,1) = mean_var(t,1);
end

% (E) grand mean average of Ei_bar
% grand mean average: average for all of the cases in all of the groups on the
% dependent variable for all temporal points
Ebar = (1/m)*sum(Ei_bar);

% end    

Eij = variable;


% numerator and denominator at each time epoch for all gait cycles
for t = 1:m
    for numGC = 1:n
        numer(t,numGC) = ((Eij(t,numGC) - Ei_bar(t,1))^2) / (m*(n-1));
        denom(t,numGC) = ((Eij(t,numGC) - Ebar)^2) / ((m*n)-1);
    end
end

for numGC = 1:n
    % numerator & denominator summed for all temporal points
    % columns represents value for each gait cycle
    numer_inner_sum(1,numGC) = sum(numer(:,numGC));
    denom_inner_sum(1,numGC) = sum(denom(:,numGC));
    
    % variance ration for each gait cycle
    VR_GC(1,numGC) = numer_inner_sum(1,numGC)/denom_inner_sum(1,numGC);
end


% outer sum of numerator and denominator for all gait cycles
numerator = sum(numer_inner_sum);
denominator = sum(denom_inner_sum);

% Variance Ratio (for all gait cycles and all time epochs)
VR = numerator/denominator;


