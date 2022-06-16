function [C3D,HeaderInfo, Point_xyz, PointLabels, GaitEvents, sortedPoints] = analyze_C3D(FullFileName, side)
% -------------------------------------------------------------------------
% ANALYZES C3D FILE
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: June 21, 2011
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	
%
% Input:        FullFileName        user defined C3D file to be analyzed% 
%               side                side to be analyzed
%               
%
% Output:       C3D                   structure of processed data                                    
%               C3D.name              name of each variable
%                                          processed, Kinetics or Kinematics                        
%               C3D.mean_GC           mean value at each % of GC                      
%               C3D.stddev_GC         standard deviation at each % of GC
%               C3D.meanVar           mean value over the entire GC
%


%% DEBUG LOOP
% close all
% clear all
% clc
% 
% [FullFileName, pathname] = uigetfile('*.c3d', 'Please select a C3D file to be analyzed');
% 
% side = 'bilateral';
% 
%%%%% END DEBUG LOOP %%%%%

%% Opens and Extract C3D file's data
[C3Dfile] = OpenC3D(FullFileName);

% Open and Read the user specified C3D file
% and returns: HeaderInfo, Point_xyz, PointLabels
[HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);

%% Organizes Kinematic Angles
% Define marker set
MkrSet_type = 'lowerbody';
alt_str = '';

% sorts kinematic angles and kinetics into a structure
[sortedPoints] = sortPoints(Point_xyz, PointLabels, MkrSet_type, alt_str);

% Modify Kinetics units
Kinetics_m(1).data = sortedPoints(34).data;
Kinetics_m(2).data = sortedPoints(35).data;
Kinetics_m(3).data = sortedPoints(36).data;
Kinetics_m(4).data = sortedPoints(37).data;
Kinetics_m(5).data = sortedPoints(38).data;
Kinetics_m(6).data = sortedPoints(39).data;
Kinetics_m(7).data = sortedPoints(40).data/1000;
Kinetics_m(8).data = sortedPoints(41).data/1000;
Kinetics_m(9).data = sortedPoints(42).data/1000;
Kinetics_m(10).data = sortedPoints(43).data/1000;
Kinetics_m(11).data = sortedPoints(44).data/1000;
Kinetics_m(12).data = sortedPoints(45).data/1000;
Kinetics_m(13).data = sortedPoints(46).data;
Kinetics_m(14).data = sortedPoints(47).data;
Kinetics_m(15).data = sortedPoints(48).data;
Kinetics_m(16).data = sortedPoints(49).data;
Kinetics_m(17).data = sortedPoints(50).data;
Kinetics_m(18).data = sortedPoints(51).data;

Kinetics_m(1).name = sortedPoints(34).name;
Kinetics_m(2).name = sortedPoints(35).name;
Kinetics_m(3).name = sortedPoints(36).name;
Kinetics_m(4).name = sortedPoints(37).name;
Kinetics_m(5).name = sortedPoints(38).name;
Kinetics_m(6).name = sortedPoints(39).name;
Kinetics_m(7).name = sortedPoints(40).name;
Kinetics_m(8).name = sortedPoints(41).name;
Kinetics_m(9).name = sortedPoints(42).name;
Kinetics_m(10).name = sortedPoints(43).name;
Kinetics_m(11).name = sortedPoints(44).name;
Kinetics_m(12).name = sortedPoints(45).name;
Kinetics_m(13).name = sortedPoints(46).name;
Kinetics_m(14).name = sortedPoints(47).name;
Kinetics_m(15).name = sortedPoints(48).name;
Kinetics_m(16).name = sortedPoints(49).name;
Kinetics_m(17).name = sortedPoints(50).name;
Kinetics_m(18).name = sortedPoints(51).name;

% Rebuild variables for processing
Variables = [sortedPoints(22:33) Kinetics_m];


%% Gait Cycle
% Find mean ensemble of each gait cycle for the c3d file
FirstFrame = HeaderInfo(4).data;
LastFrame = HeaderInfo(5).data;

% Calculate times for left & right gait events and creates matrix of each
% side's gait cycles
[Events_times, LGaitCycles, RGaitCycles] = GaitEventTimes(GaitEvents);


%% Process 
% Sample Angles at 1% Increments of Gait Cycle (101pts)
new_samp_size = 101;

% LEFT SIDE
for n = 1:2:size(Variables,2)  % number of variables to process
    for p = 1:3  % for each plane (sag, coron, transv)
        %If variable exists then process, if not skip
        if (size(Variables(n).data,3) > 1)  % if data exhists in current variable
            % if (strcmp(side, 'left') == 1)
            % name of current (newly processed) variable
            C3D(n).VarName = Variables(n).name;
            
            % Parse current Variable (ie. kinematic, kinetic) into % gait cycle (0-100%)
            variable = Variables(n).data(p,:,:);  % Extract current variable's sagittal, coronal, & rotational angles
            GaitCycle_times = LGaitCycles;
            [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);
            VarPlane_GC(p).left = var_GCs;  % parsed gait cycles for variable n
            C3D(n).VarPlane_GCs = VarPlane_GC;
            % Calculate mean and standard deviation
            variable = C3D(n).VarPlane_GCs(p).left;  % parsed gait cycles for variable(n), for each plane p
            [mean_var, std_dev] = MeanEnsemble(variable);
            mean_GC(p).left = mean_var;  % mean of variable(n) in plane 'p' for each % of gait cycle
            C3D(n).mean_GC = mean_GC;
            stddev_GC(p).left = std_dev;  % std dev of variable(n) in plane 'p' for each % of gait cycle
            C3D(n).stddev_GC = stddev_GC;
            average_of_mean(1:101,1) = mean(mean_var); % overall mean of variable
            C3D(n).meanVar(p).left = average_of_mean;
            % Calculate variance ratio
            [VR, VR_GC] = VarianceRatio(variable, mean_var);
            VarPlane_varratio(p).left = VR;
            C3D(n).VarPlane_varratio = VarPlane_varratio;
            VarPlane_varratioGC(p).left = VR_GC;
            C3D(n).VarPlane_varratioGC = VarPlane_varratioGC;
            % Calculate residual
            [resid] = Residual(mean_var, variable);
            C3D(n).VarPlane_resid(p).left = resid;
        end
    end
end

% RIGHT SIDE
for n = 2:2:size(Variables,2)  % number of variables to process
    for p = 1:3  % for each plane (sag, coron, transv)
        % If variable exists then process, if not skip
        if (size(Variables(n).data,3) > 1)  % if data exhists in current variable
            % name of current (newly processed) variable
            C3D(n).VarName = Variables(n).name;
            
            % Parse current Variable (ie. kinematic, kinetic) into % gait cycle (0-100%)
            variable = Variables(n).data(p,:,:);  % Extract current variable's sagittal, coronal, & rotational angles
            GaitCycle_times = RGaitCycles;
            [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);
            VarPlane_GC(p).right = var_GCs;  % parsed gait cycles for variable n
            C3D(n).VarPlane_GCs = VarPlane_GC;
            % Calculate mean and stddev
            variable = C3D(n).VarPlane_GCs(p).right;  % parsed gait cycles for variable(n), for each plane p
            [mean_var, std_dev] = MeanEnsemble(variable);
            mean_GC(p).right = mean_var;  % mean of variable(n) in plane 'p'
            C3D(n).mean_GC = mean_GC;
            stddev_GC(p).right = std_dev;  % std dev of variable(n) in plane 'p'
            C3D(n).stddev_GC = stddev_GC;
            average_of_mean(1:101,1) = mean(mean_var); % overall mean of variable
            C3D(n).meanVar(p).right = average_of_mean;
            % Calculate variance ratio
            [VR, VR_GC] = VarianceRatio(variable, mean_var);
            VarPlane_varratio(p).right = VR;
            C3D(n).VarPlane_varratio = VarPlane_varratio;
            VarPlane_varratioGC(p).right = VR_GC;
            C3D(n).VarPlane_varratioGC = VarPlane_varratioGC;
            % Calculate residual
            [resid] = Residual(mean_var, variable);
            C3D(n).VarPlane_resid(p).right = resid;
        end
    end
end
                
%             elseif (strcmp(side, 'bilateral') == 1)
%                 % name of current (newly processed) variable
%                 C3D(n).VarName = Variables(n).name;
% 
%                 % LEFT SIDE
%                 % Parse current Variable (ie. kinematic, kinetic) into % gait cycle (0-100%)
%                 variable = Variables(n).data(p,:,:);  % Extract current variable's sagittal, coronal, & rotational angles
%                 GaitCycle_times = LGaitCycles;  % designate side of gait cycle times
%                 [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);  % calculate as % gait cycle
%                 VarPlane_GC(p).left = var_GCs;  % parsed as % gait cycle for variable(n), for current plane (p)
%                 C3D(n).VarPlane_GCs = VarPlane_GC;
%                 % Calculate mean and stddev
%                 variable = C3D(n).VarPlane_GCs(p).left;  % parsed gait cycles for variable(n), for each plane p
%                 [mean_var, std_dev] = MeanEnsemble(variable);
%                 mean_GC(p).left = mean_var;  % mean of variable(n) in plane 'p'
%                 C3D(n).mean_GC = mean_GC;  
%                 stddev_GC(p).left = std_dev;  % std dev of variable(n) in plane 'p'
%                 C3D(n).stddev_GC = stddev_GC; 
%                 average_of_mean(1:101,1) = mean(mean_var); % overall mean of variable
%                 C3D(n).meanVar(p).left = average_of_mean;
%                 % Calculate variance ratio
%                 [VR, VR_GC] = VarianceRatio(variable, mean_var);
%                 VarPlane_varratio(p).left = VR; 
%                 C3D(n).VarPlane_varratio = VarPlane_varratio;
%                 VarPlane_varratioGC(p).left = VR_GC;
%                 C3D(n).VarPlane_varratioGC = VarPlane_varratioGC;
%                 % Calculate residual
%                 [resid] = Residual(mean_var, variable);
%                 C3D(n).VarPlane_resid(p).left = resid;
%                 
%                 
%                 % RIGHT SIDE
%                 % name of current (newly processed) variable
%                 C3D(n).VarName = Variables(n).name;
%                 
%                 % Parse current Variable (ie. kinematic, kinetic) into % gait cycle (0-100%)
%                 variable = Variables(n).data(p,:,:);  % Extract current variable's sagittal, coronal, & rotational angles
%                 GaitCycle_times = RGaitCycles;
%                 [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable);
%                 VarPlane_GC(p).right = var_GCs;  % parsed gait cycles for variable n
%                 C3D(n).VarPlane_GCs = VarPlane_GC;
%                 % Calculate mean and stddev
%                 variable = C3D(n).VarPlane_GCs(p).right;  % parsed gait cycles for variable(n), for each plane p
%                 [mean_var, std_dev] = MeanEnsemble(variable);
%                 mean_GC(p).right = mean_var;  % mean of variable(n) in plane 'p'
%                 C3D(n).mean_GC = mean_GC;  
%                 stddev_GC(p).right = std_dev;  % std dev of variable(n) in plane 'p'
%                 C3D(n).stddev_GC = stddev_GC;  
%                 average_of_mean(1:101,1) = mean(mean_var); % overall mean of variable
%                 C3D(n).meanVar(p).right = average_of_mean;
%                 % Calculate variance ratio
%                 [VR, VR_GC] = VarianceRatio(variable, mean_var);
%                 VarPlane_varratio(p).right = VR; 
%                 C3D(n).VarPlane_varratio = VarPlane_varratio;
%                 VarPlane_varratioGC(p).right = VR_GC;
%                 C3D(n).VarPlane_varratioGC = VarPlane_varratioGC;
%                 % Calculate residual
%                 [resid] = Residual(mean_var, variable);
%                 C3D(n).VarPlane_resid(p).right = resid;
%                 
% %             end
%             
%         end
%     end
% end


%% Temporal Spatial Variables for Trial
% Calculate temporal spatial measures for current c3d file
[TempSpat] = TemporalSpatial(HeaderInfo, Point_xyz, PointLabels, GaitEvents);
trial_level = size(C3D,2)+1;
C3D(trial_level).TempSpat = TempSpat;


%% Variance Ratios & Residuals
% Composite matrix of all variance ratios and residuals for C3D file's lower body kinematics
C3D(trial_level).VarName = 'Trial Level';

if (strcmp(side, 'left') == 1)
    % LEFT
    % Variance ratio & residuals for each left gait cycle, lower body variables
    for n = 1:length(LGaitCycles(1,1,:))
        LGC(n).varratio = [C3D(1).VarPlane_varratioGC(1).left(n); C3D(1).VarPlane_varratioGC(2).left(n); C3D(1).VarPlane_varratioGC(3).left(n);
                           C3D(3).VarPlane_varratioGC(1).left(n); C3D(3).VarPlane_varratioGC(2).left(n); C3D(3).VarPlane_varratioGC(3).left(n);
                           C3D(5).VarPlane_varratioGC(1).left(n); C3D(5).VarPlane_varratioGC(2).left(n); C3D(5).VarPlane_varratioGC(3).left(n);
                           C3D(7).VarPlane_varratioGC(1).left(n); C3D(9).VarPlane_varratioGC(3).left(n); C3D(7).VarPlane_varratioGC(3).left(n)];

        LGC(n).resid = [C3D(1).VarPlane_resid(1).left(n); C3D(1).VarPlane_resid(2).left(n); C3D(1).VarPlane_resid(3).left(n);
                        C3D(3).VarPlane_resid(1).left(n); C3D(3).VarPlane_resid(2).left(n); C3D(3).VarPlane_resid(3).left(n);
                        C3D(5).VarPlane_resid(1).left(n); C3D(5).VarPlane_resid(2).left(n); C3D(5).VarPlane_resid(3).left(n);
                        C3D(7).VarPlane_resid(1).left(n); C3D(9).VarPlane_resid(3).left(n); C3D(7).VarPlane_resid(3).left(n)];
    end
    
    % Sum all left gait cycle variance ratios & residuals
    for j = 1:length(LGaitCycles(1,1,:))
        LGCs_VR_sums(j) = sum(LGC(j).varratio);
        LGCs_resid_sums(j) = sum(LGC(j).resid);
    end
    VR_trial_left = sum(LGCs_VR_sums);
    
    % Sum of left and right side variance ratios
    variance_ratio_trial = VR_trial_left;

    % TRIAL'S TOTAL, LEFT VARIANCE RATIO
    C3D(trial_level).variance_ratio_trial = variance_ratio_trial;

    % Sort left side's residuals
    [resid_val, gaitcycle] = sort(LGCs_resid_sums, 'ascend');
    % provides the representative gait cycle, in ascending order
    residuals.left = gaitcycle;

    % TRIAL'S RESIDUALS FOR EACH SIDE'S GAIT CYCLES
    % listed in order of highest rank to lowest
    C3D(trial_level).GCranking = residuals;

elseif (strcmp(side, 'right') == 1)
    % RIGHT
    % Variance ratio & residuals for each right gait cycle, lower body variables
    for r = 1:length(RGaitCycles(1,1,:))
        RGC(r).varratio = [C3D(2).VarPlane_varratioGC(1).right(r); C3D(2).VarPlane_varratioGC(2).right(r); C3D(2).VarPlane_varratioGC(3).right(r);
                           C3D(4).VarPlane_varratioGC(1).right(r); C3D(4).VarPlane_varratioGC(2).right(r); C3D(4).VarPlane_varratioGC(3).right(r);
                           C3D(6).VarPlane_varratioGC(1).right(r); C3D(6).VarPlane_varratioGC(2).right(r); C3D(6).VarPlane_varratioGC(3).right(r);
                           C3D(8).VarPlane_varratioGC(1).right(r); C3D(10).VarPlane_varratioGC(3).right(r); C3D(8).VarPlane_varratioGC(3).right(r)];

        RGC(r).resid = [C3D(2).VarPlane_resid(1).right(r); C3D(2).VarPlane_resid(2).right(r); C3D(2).VarPlane_resid(3).right(r);
                        C3D(4).VarPlane_resid(1).right(r); C3D(4).VarPlane_resid(2).right(r); C3D(4).VarPlane_resid(3).right(r);
                        C3D(6).VarPlane_resid(1).right(r); C3D(6).VarPlane_resid(2).right(r); C3D(6).VarPlane_resid(3).right(r);
                        C3D(8).VarPlane_resid(1).right(r); C3D(10).VarPlane_resid(3).right(r); C3D(8).VarPlane_resid(3).right(r)];
    end
    
    % Sum all right gait cycle variance ratios & residuals
    for w = 1:length(RGaitCycles(1,1,:))
        RGCs_VR_sums(w) = sum(RGC(w).varratio);
        RGCs_resid_sums(w) = sum(RGC(w).resid);
    end
    VR_trial_right = sum(RGCs_VR_sums);

    % Sum of right side variance ratios
    variance_ratio_trial = VR_trial_right;

    % TRIAL'S TOTAL, RIGHT VARIANCE RATIO
    C3D(trial_level).variance_ratio_trial = variance_ratio_trial;

    % Sort right side's residuals
    [resid_val, gaitcycle] = sort(RGCs_resid_sums, 'ascend');
    % provides the representative gait cycle, in ascending order
    residuals.right = gaitcycle;

    % TRIAL'S RESIDUALS FOR EACH SIDE'S GAIT CYCLES
    % listed in order of highest rank to lowest
    C3D(trial_level).GCranking = residuals;
    
else
    % BILATERAL
    
    % Variance ratio & residuals for each left gait cycle, lower body variables
    for n = 1:length(LGaitCycles(1,1,:))
        LGC(n).varratio = [C3D(1).VarPlane_varratioGC(1).left(n); C3D(1).VarPlane_varratioGC(2).left(n); C3D(1).VarPlane_varratioGC(3).left(n);
                           C3D(3).VarPlane_varratioGC(1).left(n); C3D(3).VarPlane_varratioGC(2).left(n); C3D(3).VarPlane_varratioGC(3).left(n);
                           C3D(5).VarPlane_varratioGC(1).left(n); C3D(5).VarPlane_varratioGC(2).left(n); C3D(5).VarPlane_varratioGC(3).left(n);
                           C3D(7).VarPlane_varratioGC(1).left(n); C3D(9).VarPlane_varratioGC(3).left(n); C3D(7).VarPlane_varratioGC(3).left(n)];

        LGC(n).resid = [C3D(1).VarPlane_resid(1).left(n); C3D(1).VarPlane_resid(2).left(n); C3D(1).VarPlane_resid(3).left(n);
                        C3D(3).VarPlane_resid(1).left(n); C3D(3).VarPlane_resid(2).left(n); C3D(3).VarPlane_resid(3).left(n);
                        C3D(5).VarPlane_resid(1).left(n); C3D(5).VarPlane_resid(2).left(n); C3D(5).VarPlane_resid(3).left(n);
                        C3D(7).VarPlane_resid(1).left(n); C3D(9).VarPlane_resid(3).left(n); C3D(7).VarPlane_resid(3).left(n)];
    end

    % Variance ratio & residuals for each right gait cycle, lower body variables
    for r = 1:length(RGaitCycles(1,1,:))
        RGC(r).varratio = [C3D(2).VarPlane_varratioGC(1).right(r); C3D(2).VarPlane_varratioGC(2).right(r); C3D(2).VarPlane_varratioGC(3).right(r);
                           C3D(4).VarPlane_varratioGC(1).right(r); C3D(4).VarPlane_varratioGC(2).right(r); C3D(4).VarPlane_varratioGC(3).right(r);
                           C3D(6).VarPlane_varratioGC(1).right(r); C3D(6).VarPlane_varratioGC(2).right(r); C3D(6).VarPlane_varratioGC(3).right(r);
                           C3D(8).VarPlane_varratioGC(1).right(r); C3D(10).VarPlane_varratioGC(3).right(r); C3D(8).VarPlane_varratioGC(3).right(r)];

        RGC(r).resid = [C3D(2).VarPlane_resid(1).right(r); C3D(2).VarPlane_resid(2).right(r); C3D(2).VarPlane_resid(3).right(r);
                        C3D(4).VarPlane_resid(1).right(r); C3D(4).VarPlane_resid(2).right(r); C3D(4).VarPlane_resid(3).right(r);
                        C3D(6).VarPlane_resid(1).right(r); C3D(6).VarPlane_resid(2).right(r); C3D(6).VarPlane_resid(3).right(r);
                        C3D(8).VarPlane_resid(1).right(r); C3D(10).VarPlane_resid(3).right(r); C3D(8).VarPlane_resid(3).right(r)];
    end

    % Combine left and right GC variance ratios and residuals into one
    % large structure
    for nL = 1:length(LGaitCycles(1,1,:));
        AllGC_VRs(:,nL) = LGC(nL).varratio;
        AllGC_resids(:,nL) = LGC(nL).resid;
    end
    
    for nRGCs = 1:length(RGaitCycles(1,1,:))
        current_RGC_VR = RGC(nRGCs).varratio;
        current_RGC_resid = RGC(nRGCs).resid;
        nR = length(LGaitCycles(1,1,:)) + nRGCs;
        AllGC_VRs(:,nR) = current_RGC_VR;
        AllGC_resids(:,nR) = current_RGC_resid;
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
    
    % Sum the both left and right gait cycle's VR values
    AllGC_VR_sums = sum(AllGC_VRs);
    AllGC_resid_sums = sum(AllGC_resids);
    
    % TRIAL'S TOTAL, BILATERAL VARIANCE RATIO
    % Sum of left and right side variance ratios
    C3D(trial_level).variance_ratio_trial = AllGC_VR_sums;
    
    
    
    % Sort left side's residuals
    [resid_val, gaitcycle] = sort(LGCs_resid_sums, 'ascend');
    % provides the representative gait cycle, in ascending order
    residuals.left = gaitcycle;

    % Sort right side's residuals
    [resid_val, gaitcycle] = sort(RGCs_resid_sums, 'ascend');
    % provides the representative gait cycle, in ascending order
    residuals.right = gaitcycle;
    
    % Sort all gait cycles residuals
    [resid_val, gaitcycle] = sort(AllGC_resid_sums, 'ascend');
    % add a string (L or R) for which side the gait cycle belongs
    for b = 1:length(gaitcycle)
        current_num_GC = gaitcycle(1,b);
        
        if (current_num_GC <= length(LGaitCycles(1,1,:)))
            % if current gaitcycle value is equal to a LGC, add 'L'
            side_str = 'L';  GC = num2str(current_num_GC);
            id_gaitcycle(1,:,b) = strcat(side_str,GC);
        else
            % if current gaitcycle value is greater than number of left
            % GCs, add 'R' for right side and determine which RGC number
            side_str = 'R';  GC = num2str(current_num_GC - (length(LGaitCycles(1,1,:))));
            id_gaitcycle(1,:,b) = strcat(side_str,GC);
        end
    end
        
    residuals.allGCs = id_gaitcycle;
    
    % TRIAL'S RESIDUALS FOR EACH SIDE'S GAIT CYCLES
    % listed in order of highest rank to lowest
    C3D(trial_level).GCranking = residuals;
end

