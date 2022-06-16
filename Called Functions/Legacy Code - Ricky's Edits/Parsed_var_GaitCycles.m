function [var_GCs] = Parsed_var_GaitCycles(HeaderInfo, GaitCycle_times, new_samp_size, variable)

% -------------------------------------------------------------------------
% PARSED VARIABLE BY GAIT CYCLE
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: July 15, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This is function parses the input variable by the times for
%               each side's gait cycle(s).
%
% Input:        HeaderInfo       Header information of c3d file
%               LGaitCycles      Matrix of left gait cycle times
%               RGaitCycles      Matrix of right gait cycle times
%               new_samp_size    lenght variable (resampled size)
%               variable         input argument to be interpolated and
%                                resampled as 0-100% of gait cycle
%
% OUTPUTS:  
%               var_GCs          variable that is represented as percent of
%                                the gait cycle
%


%% DEBUG LOOP

% close all
% clear all
% clc
% 
% % Prompts User for Input
% [FullFileName, pathname] = uigetfile('*.*', 'Please select C3D file to be analyzed');
% 
% % Opens and Extract C3D file's data
% [C3Dfile] = OpenC3D(FullFileName);
% 
% % Open and Read the user specified C3D file
% % and returns: HeaderInfo, Point_xyz, PointLabels
% [HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);
% 
% % Organizes Kinematic Angles
% % sorts kinematic angles into a structure
% [Kinematics] = BSPT_sortKinematics(Point_xyz, PointLabels);
% 
% LPelvisAngles = Kinematics(1).data;
% RPelvisAngles = Kinematics(2).data;
% LThoraxAngles = Kinematics(13).data;
% RThoraxAngles = Kinematics(14).data;
% 
% % Calculate times for left & right gait events and creates matrix of each
% % side's gait cycles
% [LFS_times, LFO_times, RFS_times, RFO_times, LGaitCycles, RGaitCycles] = BSPT_GaitEventTimes(GaitEvents);
% 
% new_samp_size = 101;

% % Example data
% variable = LThoraxAngles(1,:,:);
% variable = LPelvisAngles(1,:,:);
% GaitCycle_times = LGaitCycles;
% variable = RPelvisAngles(1,:,:);
% GaitCycle_times = RGaitCycles;




%% Organize Input Variable by Gait Cycle(s)

% Gets first frame of trial
FirstFrame = HeaderInfo(4).data;

% GAIT CYCLE(S)
N_GaitCycles = length(GaitCycle_times(1,1,:));

if N_GaitCycles ~= 0
    % Find frames for each gait cycle
    for numGC = 1:N_GaitCycles
        % start & end times for current gait cycle
        start_GC(numGC) = GaitCycle_times(1,1,numGC);
        end_GC(numGC) = GaitCycle_times(1,2,numGC);
        
    end  
    
    % Counting from the first frame of data, finds the gait cycle
    % start time and adjustes so counting by total number of points in varargin
    for t = 1:length(GaitCycle_times(1,1,:))
        if FirstFrame == GaitCycle_times(1,:,1)
            adj_GCtimes(:,:,t) = GaitCycle_times(:,:,t);
        else
            adj_GCtimes(:,:,t) = GaitCycle_times(1,:,t)-FirstFrame;
        end
            
    end
    
    % Parses out the input variable into each gait cycle
    if isempty(GaitCycle_times) ~= 1
        % Parsed_var_GCs must be a structure array so can handle varying
        % data lengths for each gait cycle

        for numGC = 1:N_GaitCycles
            % number of gait cycles  
            Parsed_var_GCs(numGC).gaitcycle = numGC;
            
            % data for each gait cycle (of varying length)
            Parsed_var_GCs(numGC).data = variable(:,(adj_GCtimes(1,1,numGC):adj_GCtimes(1,2,numGC)) );
        end            

    end

end


%% Variable by Percent Gait Cycle
for numGC = 1:N_GaitCycles
    % Interpolates a curve to fit variable for each gait cycle
    Parsed_var_GC = Parsed_var_GCs(numGC).data;

    Parsed_var_x = 1:length(Parsed_var_GC); % data points to interpolate from
    Parsed_var_y = Parsed_var_GC(1,:); % data points to interpolate from
    Parsed_var_xinterp = 1:0.1:length(Parsed_var_GC(1,:)); % new sample rate (interpolated curve's x points)
    Parsed_var_yinterp = interp1(Parsed_var_x, Parsed_var_y, Parsed_var_xinterp, 'cubic'); % interpolated curve's y points
    
    % Resamples data from the interpolated curve so each point of the
    % gait cycle represents a percent (0-100%) of the gait cycle   
    var_GCs(:,numGC) = interpft(Parsed_var_GC, new_samp_size);

end


    

