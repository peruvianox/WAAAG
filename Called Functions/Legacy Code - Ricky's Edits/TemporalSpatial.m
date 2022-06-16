function [TempSpat] = TemporalSpatial(HeaderInfo, Point_xyz, PointLabels, GaitEvents)
% -------------------------------------------------------------------------
% TEMPORAL SPATIAL CALCULATIONS
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: Setpember 16, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Indexing fixed 23Dec2015
%
% Description:	This function calculates temporal spatial variabls for the
%               C3D file.
%
% Input:        GaitEvents          gait events from c3d file
%               
%
% Output:       TempSpat            structure of temporal spatial variables
%
%         

% %% DEBUG LOOP
% close all
% clear all
% clc
% 
% % Prompts User for Input
% [FullFileName, pathname] = uigetfile('*.c3d*', 'Please select C3D file to be analyzed');
% 
% [C3Dfile] = OpenC3D(FullFileName);
% 
% % Open and Read the user specified C3D file
% % and returns: HeaderInfo, Point_xyz, PointLabels
% [HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);

%%%%% END DEBUG LOOP %%%%%

%% Data Extraction
% Markers
MkrSet_type = 'lowerbody';  alt_str = '';
[search_Pt_label] = defPoints_Labels(MkrSet_type);
[sortedPoints] = sortPoints(Point_xyz, PointLabels, MkrSet_type, alt_str);
LHEE=searchPoint(Point_xyz, PointLabels, 'LHEE'); % left heel marker
RHEE = searchPoint(Point_xyz, PointLabels, 'RHEE');  % right heel marker

% First & Last Frames
FirstFrame = HeaderInfo(4).data;
LastFrame = HeaderInfo(5).data;

% Calculates GaitEvents, Foot Strike & Foot Off Times
[Events_times, LGaitCycles, RGaitCycles] = GaitEventTimes(GaitEvents);
LFS_times = Events_times(1).data;
LFO_times = Events_times(2).data;
RFS_times = Events_times(4).data;
RFO_times = Events_times(5).data;

% Adjust Foot Strike times so in context of FirstFrame value
if FirstFrame < LFS_times(1)
    for n = 1:length(LFS_times)    
        adj_LFS_times(n) = LFS_times(n) - FirstFrame;
    end   
else
    for n = 1:length(LFS_times)    
        adj_LFS_times(n) = LFS_times(n);
    end
end
% Number of left gait cycles in trial
N_LGaitCycles = length(LFS_times)-1;
if N_LGaitCycles ~= 0
    for numGC = 1:N_LGaitCycles
        start_LGC(numGC) = adj_LFS_times(numGC);
        for numGC_end = 1:N_LGaitCycles
            end_LGC(numGC_end) = adj_LFS_times(numGC_end+1);
        end
        % adjusted left gait cycle times (foot strikes)
        adj_LGaitCycles_t(:,:,numGC) = [start_LGC(numGC), end_LGC(numGC)];  
    end
end
if FirstFrame < RFS_times(1)
    for m = 1:length(RFS_times)
        adj_RFS_times(m) = RFS_times(m) - FirstFrame;
    end
else
    for m = 1:length(RFS_times)
        adj_RFS_times(m) = RFS_times(m);
    end
end
% Number of right gait cycles in trial
N_RGaitCycles = length(RFS_times)-1;
if N_RGaitCycles ~= 0
    for numGC = 1:N_RGaitCycles
        start_RGC(numGC) = adj_RFS_times(numGC);
        for numGC_end = 1:N_RGaitCycles
            end_RGC(numGC_end) = adj_RFS_times(numGC_end+1);
        end
        % adjusted right gait cycle times (foot strikes)
        adj_RGaitCycles_t(:,:,numGC) = [start_RGC(numGC), end_RGC(numGC)];
    end
end
% Time from first foot strike to last foot strike
if adj_LGaitCycles_t(1,1,1) < adj_RGaitCycles_t(1,1,1)
    % if first foot strike occurs on left foot
    if adj_LGaitCycles_t(1,2,end) > adj_RGaitCycles_t(1,2,end)
        % if last foot strike occurs on left foot
        frame_time = (adj_LGaitCycles_t(1,2,end) - adj_LGaitCycles_t(1,1,1))*(1/120);  % seconds
        % distance traveled (m) by heel marker for duration of time period considered
        dist_travel = (LHEE(1,1,adj_LGaitCycles_t(1,2,end)) - LHEE(1,1,adj_LGaitCycles_t(1,1,1)))/1000;
    else
        % if last foot strike occurs on right foot
        frame_time = (adj_RGaitCycles_t(1,2,end) - adj_LGaitCycles_t(1,1,1))*(1/120);
        % distance traveled (m) by heel marker for duration of time period considered
        dist_travel = (RHEE(1,1,adj_RGaitCycles_t(1,2,end)) - LHEE(1,1,adj_LGaitCycles_t(1,1,1)))/1000;
    end
    % stride length (m) distance between 2 successive placements of same foot
    % first left foot strike to last left foot strike
    for t = 1:length(LGaitCycles(1,1,:))
        stride_length(:,:,t) = abs(LHEE(1,1,adj_LGaitCycles_t(1,2,t)) - LHEE(1,1,adj_LGaitCycles_t(1,1,t)))/1000;
    end
else
    % if first foot strike occurs on right foot
    if adj_LGaitCycles_t(1,2,end) > adj_RGaitCycles_t(1,2,end)
        % if last foot strike occurs on left foot
        % frame_time (sec) from first to last foot strike
        frame_time = (adj_LGaitCycles_t(1,2,end) - adj_RGaitCycles_t(1,1,1))*(1/120);
        % distance traveled (m) by heel marker for duration of time period considered
        dist_travel = (LHEE(1,1,adj_LGaitCycles_t(1,2,end)) - RHEE(1,1,adj_RGaitCycles_t(1,1,1)))/1000;
    else
        % if last foot strike occurs on right foot
        frame_time = (adj_RGaitCycles_t(1,2,end) - adj_RGaitCycles_t(1,1,1))*(1/120);   
        % distance traveled (m) by heel marker for duration of time period considered
        dist_travel = (RHEE(1,1,adj_RGaitCycles_t(1,2,end)) - RHEE(1,1,adj_RGaitCycles_t(1,1,1)))/1000;  
    end
    % stride length (m) first right foot strike to last right foot strike
    for t = 1:length(RGaitCycles(1,1,:))
        stride_length(:,:,t) = abs(RHEE(1,1,adj_RGaitCycles_t(1,2,t)) - RHEE(1,1,adj_RGaitCycles_t(1,1,t)))/1000;
    end
end


%% CADENCE (steps/min)
% Steps = number of steps in trial 
% where total foot strikes - 1 because 1st foot strike indicates beginning of measurement
N_steps = (length(LFS_times)+length(RFS_times))-1;

% Subject's cadence (steps/minute)
Cadence = N_steps*60/frame_time;

%% WALKING SPEED (m/min)
% distance traveled by heel marker throughout entire time period considered
% (first foot strike to last foot strike, aka frame_time)
WalkSpeed = abs(dist_travel)/frame_time;  % meters/second
Walking_Speed = WalkSpeed*60;  % meters/minute

%% STRIDE LENGTH (meters)
% distance between two successive placements of the same foot (meters)
% subject's stride length is average for all strides for time period considered
Ave_Stride_Length = mean(stride_length);

%% STRIDE TIME (seconds)
% duration of time for each side's average stride length

%%%%% LEFT %%%%%
for n = 1:length(adj_LGaitCycles_t(1,1,:))
    % times for each left gait cycle (seconds)
    Lframe_times(n) = (adj_LGaitCycles_t(1,2,n) - adj_LGaitCycles_t(1,1,n))*(1/120);
end
% Average Left Stride Time (sec)
Ave_LStride_Time = sum(Lframe_times)/length(Lframe_times);
%%%%% RIGHT %%%%%
for n = 1:length(adj_RGaitCycles_t(1,1,:))
    % times for each right gait cycle (seconds)
    Rframe_times(n) = (adj_RGaitCycles_t(1,2,n) - adj_RGaitCycles_t(1,1,n))*(1/120);
end
% Average Right Stride Time (sec)
Ave_RStride_Time = sum(Rframe_times)/length(Rframe_times);

%% STEP TIME (seconds), LENGTH & WIDTH (centimeters)

%%%%% LEFT %%%%%
% duration of time from right foot strike to next left foot strike
if adj_LFS_times(1) < adj_RFS_times(1)
    % Step Time
    % if left foot strike occurs first,
    % ignore first left foot strike and count from first right foot strike
    new_LFS_times = adj_LFS_times(2:end);
    for n = 1:length(new_LFS_times)
        % times for each left step (seconds)
        Lstep_times(n) = abs(adj_RFS_times(n) - new_LFS_times(n))*(1/120);
    end
    % if uneven number of left and right foot strike times, adjust times
    % ????? Not sure if needed, given new_LFS_times when would this occur????
    if length(adj_RFS_times) < length(new_LFS_times)
        % if more left foot strike times than right, remove last left
        new_LFS_times = new_LFS_times(1:end-1);
    end
    % Step Length
    if length(adj_RFS_times) > length(new_LFS_times)
        % if there are more right foot strikes than (new adjusted) left
        % AGAIN NOT SURE IF NEEDED
        for t = 1:length(new_LFS_times)
            Lstep_lengths(t) = abs(RHEE(1,1,adj_RFS_times(t)) - LHEE(1,1,new_LFS_times(t)))/1000;
            Lstep_widths(t) = abs(RHEE(2,1,adj_RFS_times(t)) - LHEE(2,1,new_LFS_times(t)))/1000;
        end
    else
        % if there are more left foot strikes than right
        for t = 1:length(adj_RFS_times)
            Lstep_lengths(t) = abs(RHEE(1,1,adj_RFS_times(t)) - LHEE(1,1,new_LFS_times(t)))/1000;
            Lstep_widths(t) = abs(RHEE(2,1,adj_RFS_times(t)) - LHEE(2,1,new_LFS_times(t)))/1000;  
        end
    end
        
else
    % Step Time
    % if right foot strike occurs first, count from first right foot strike
    if length(adj_RFS_times) > length(adj_LFS_times)
        % if more right foot strikes than left, count by number of lefts
        nLsteps = length(adj_LFS_times);
    else
        nLsteps = length(adj_RFS_times);
    end
    for n = 1:nLsteps
        % times for each left step (seconds)
        Lstep_times(n) = abs(adj_RFS_times(n) - adj_LFS_times(n))*(1/120);
    end
    % Step Length
    if length(adj_RFS_times) > length(adj_LFS_times)
        % if there are more right foot strikes than (new adjusted) left
        for t = 1:length(Lstep_times)
            Lstep_lengths(t) = abs(RHEE(1,1,adj_RFS_times(t)) - LHEE(1,1,adj_LFS_times(t)))/1000;
            Lstep_widths(t) = abs(RHEE(2,1,adj_RFS_times(t)) - LHEE(2,1,adj_LFS_times(t)))/1000; 
        end
    else
        % if there are more left foot strikes than right
        for t = 1:length(adj_RFS_times)
            Lstep_lengths(t) = abs(RHEE(1,1,adj_RFS_times(t)) - LHEE(1,1,adj_LFS_times(t)))/1000;
            Lstep_widths(t) = abs(RHEE(2,1,adj_RFS_times(t)) - LHEE(2,1,adj_LFS_times(t)))/1000; 
        end
    end
end

Ave_LStep_Time = sum(Lstep_times)/length(Lstep_times);
Ave_LStep_Length = abs((sum(Lstep_lengths)/length(Lstep_lengths))*100);
Std_LStep_Length = std(Lstep_lengths)*100; 

%%%%% RIGHT %%%%%
% duration of time from left foot strike to next right foot strike
if adj_RFS_times(1) < adj_LFS_times(1)
    % Step Time
    % if right foot strike occurs first
    % ignore first right foot strike and count from first left foot strike
    new_RFS_times = adj_RFS_times(2:end);
    % if uneven number of left and right foot strike times, adjust times
    if length(adj_LFS_times) < length(new_RFS_times)
        % if more right foot strike times than left, remove last right
        new_RFS_times = new_RFS_times(1:end-1);
    end
    for n = 1:length(new_RFS_times)
        % calculate times for each right step (seconds)
        Rstep_times(n) = abs(adj_LFS_times(n) - new_RFS_times(n))*(1/120);
    end   
    % Step Length and Width
    if length(adj_LFS_times) > length(new_RFS_times)
        % if there are more left foot strikes than (new adjusted) right
        for t = 1:length(new_RFS_times)
            Rstep_lengths(t) = abs(LHEE(1,1,adj_LFS_times(t)) - RHEE(1,1,new_RFS_times(t)))/1000;
            Rstep_widths(t) = abs(LHEE(2,1,adj_LFS_times(t)) - RHEE(2,1,new_RFS_times(t)))/1000;
        end
    else
        % if there are more right foot strikes than left
        for t = 1:length(adj_LFS_times)
            Rstep_lengths(t) = abs(LHEE(1,1,adj_LFS_times(t)) - RHEE(1,1,new_RFS_times(t)))/1000;
            Rstep_widths(t) = abs(LHEE(2,1,adj_LFS_times(t)) - RHEE(2,1,new_RFS_times(t)))/1000; 
        end
    end
        
    
else
    % Step Time
    % if left foot strike occurs first, count from first left foot strike
    % for number of right foot strikes
    
    if length(adj_LFS_times) > length(adj_RFS_times)
        % if more left foot strikes than right, count by number of rights
        nRsteps = length(adj_RFS_times);
    else
        nRsteps = length(adj_LFS_times);
    end
    for n = 1:nRsteps
        % calculate times for each right step (seconds)
        Rstep_times(n) = abs(adj_LFS_times(n) - adj_RFS_times(n))*(1/120);
    end
    % Step Length and Width
    if length(adj_LFS_times) > length(adj_RFS_times)
        % if there are more left foot strikes than (adjusted) right
        for t = 1:length(adj_RFS_times)
            Rstep_lengths(t) = abs(LHEE(1,1,adj_LFS_times(t)) - RHEE(1,1,adj_RFS_times(t)))/1000;
            Rstep_widths(t) = abs(LHEE(2,1,adj_LFS_times(t)) - RHEE(2,1,adj_RFS_times(t)))/1000; 
        end
    else
        % if there are more right foot strikes than left
        for t = 1:length(adj_LFS_times)
            Rstep_lengths(t) = abs(LHEE(1,1,adj_LFS_times(t)) - RHEE(1,1,adj_RFS_times(t)))/1000;
            Rstep_widths(t) = abs(LHEE(2,1,adj_LFS_times(t)) - RHEE(2,1,adj_RFS_times(t)))/1000; 
        end
    end
end


Ave_RStep_Time = sum(Rstep_times)/length(Rstep_times);
Ave_RStep_Length = abs((sum(Rstep_lengths)/length(Rstep_lengths))*100);
Std_RStep_Length = std(Rstep_lengths)*100; 

Ave_Step_Width = abs(((sum(Rstep_widths)+sum(Lstep_widths)))/(length(Rstep_widths)+length(Lstep_widths))*100);
Std_Step_Width = (std(Rstep_widths) + std(Lstep_widths))*100/2; 


%% FOOT OFF (% gait cycle)
% Find the average percent of gait cycle when foot off occurs

%%%%% LEFT %%%%%
for t = 1:length(LGaitCycles(1,1,:))
    % Find only foot off times that are in a valid gait cycle
    for n = 1:length(LFO_times)
        if (LGaitCycles(1,1,t) < LFO_times(1,n)) && (LFO_times(1,n) < LGaitCycles(1,2,t))
            LFOtimes(1,n) = LFO_times(1,n);
        end
    end
    % Calculate percent of gait cycle when foot off occurs for each gait cycle
    for m = 1:length(LFOtimes)
        percent_Lfoot_offs(:,t) = (LFOtimes(1,m)-LGaitCycles(1,1,t)) / (LGaitCycles(1,2,t)-LGaitCycles(1,1,t));
    end
end

% Average percentage of left foot off during gait cycle
Ave_Percent_LFootOff = (mean(percent_Lfoot_offs)*100);
Std_Percent_LFootOff = std(percent_Lfoot_offs)*100;

%%%%% RIGHT %%%%%
for t = 1:length(RGaitCycles(1,1,:))
    % Find only foot off times that are in a valid gait cycle
    for n = 1:length(RFO_times)
        if (RGaitCycles(1,1,t) < RFO_times(1,n)) && (RFO_times(1,n) < RGaitCycles(1,2,t))
            RFOtimes(1,n) = RFO_times(1,n);
        end
    end
    % Calculate percent of gait cycle when foot off occurs for each gait cycle
    for m = 1:length(RFOtimes)
        percent_Rfoot_offs(:,t) = (RFOtimes(1,m)-RGaitCycles(1,1,t)) / (RGaitCycles(1,2,t)-RGaitCycles(1,1,t));
    end
end

% Average percentage of right foot off during gait cycle
Ave_Percent_RFootOff = (mean(percent_Rfoot_offs)*100);
Std_Percent_RFootOff = std(percent_Rfoot_offs)*100;


%% OPPOSITE FOOT OFF (% gait cycle)
% Find the average percent of gait cycle when opposite foot off occurs

%%%%% LEFT %%%%%
if RFO_times(1) < LFS_times(1)
    % if right foot off occurs before first left foot contact
    % ignore first right foot off
    new_RFO_times = RFO_times(2:end);
else
    % else keep first right foot off time
    new_RFO_times = RFO_times;
end
if new_RFO_times(end) > LFS_times(end)
    % if last right foot off times occurs after last left foot strike
    % remove last right foot off
    new_RFO_times = new_RFO_times(1:end-1);
end
if length(new_RFO_times) <= 2
    percent_L_OppFootOff = ((new_RFO_times(1)-LFS_times(1))/(LFS_times(2)-LFS_times(1)));
elseif length(new_RFO_times) == length(LFS_times)
    for k = 1:length(new_RFO_times)-1
        percent_L_OppFootOff(k) = ((new_RFO_times(k)-LFS_times(k))/(LFS_times(k+1)-LFS_times(k)));
    end
else
    for k = 1:length(new_RFO_times)
        percent_L_OppFootOff(k) = ((new_RFO_times(k)-LFS_times(k))/(LFS_times(k+1)-LFS_times(k)));
    end
end
% Average percentage of right foot off during left gait cycle
Ave_Percent_L_OppFootOff = (mean(percent_L_OppFootOff)*100);

%%%%% RIGHT %%%%%
if LFO_times(1) < RFS_times(1)
    % if left foot off occurs before first right foot contact
    % ignore first left foot off
    new_LFO_times = LFO_times(2:end);
else
    % else keep first left foot off time
    new_LFO_times = LFO_times;
end
if new_LFO_times(end) > RFS_times(end)
    % if last left foot off occurs after last right foot strike
    % remove last left foot off
    new_LFO_times = new_LFO_times(1:end-1);
end
if length(new_LFO_times) <= 2
    percent_R_OppFootOff = ((new_LFO_times(1)-RFS_times(1))/(RFS_times(2)-RFS_times(1)));
elseif length(new_LFO_times) == length(RFS_times)
    for k = 1:length(new_LFO_times)-1
        percent_R_OppFootOff(k) = ((new_LFO_times(k)-RFS_times(k))/(RFS_times(k+1)-RFS_times(k)));
    end
else
    for k = 1:length(new_LFO_times)
        percent_R_OppFootOff(k) = ((new_LFO_times(k)-RFS_times(k))/(RFS_times(k+1)-RFS_times(k)));
    end
end

% Average percentage of left foot off during right gait cycle
Ave_Percent_R_OppFootOff = (mean(percent_R_OppFootOff)*100);

%% OPPOSITE FOOT CONTACT (% gait cycle)
% Find the average percent of gait cycle when opposite foot contact occurs

%%%%% LEFT %%%%%
if RFS_times(1) < LFS_times(1)
    % if first right foot contact occurs before first left foot contact
    % ignore first right foot contact
    new_RFS_times = RFS_times(2:end);
else
    % else keep first right foot contact time
    new_RFS_times = RFS_times;
end
if new_RFS_times(end) > LFS_times(end)
    % if last right foot strike time occurs after last left foot strike
    % remove last right foot strike
    new_RFS_times = new_RFS_times(1:end-1);
end
if length(new_RFS_times) == length(LFS_times)
    for k = 1:length(new_RFS_times)-1
        percent_L_OppFootOn(k) = ((new_RFS_times(k)-LFS_times(k))/(LFS_times(k+1)-LFS_times(k)));
    end
else
    for k = 1:length(new_RFS_times)
        percent_L_OppFootOn(k) = ((new_RFS_times(k)-LFS_times(k))/(LFS_times(k+1)-LFS_times(k)));
    end
end
% Average percentage of right foot off during left gait cycle
Ave_Percent_L_OppFootOn = (mean(percent_L_OppFootOn)*100);
%%%%% RIGHT %%%%%
if LFS_times(1) < RFS_times(1)
    % if first left foot contact occurs before first right foot contact
    % ignore first left foot contact
    new_LFS_times = LFS_times(2:end);
else
    % else keep first left foot strike time
    new_LFS_times = LFS_times;
end
if new_LFS_times(end) > RFS_times(end)
    % if last left foot strike time occurs after last right foot strike
    % remove last left foot strike
    new_LFS_times = new_LFS_times(1:end-1);
end
if length(new_LFS_times) == length(RFS_times)
    for k = 1:length(new_LFS_times)-1
        percent_R_OppFootOn(k) = ((new_LFS_times(k)-RFS_times(k))/(RFS_times(k+1)-RFS_times(k))); % frames
    end
else
    for k = 1:length(new_LFS_times)
        percent_R_OppFootOn(k) = ((new_LFS_times(k)-RFS_times(k))/(RFS_times(k+1)-RFS_times(k))); % frames
    end
end
% Average percentage of right foot off during left gait cycle
Ave_Percent_R_OppFootOn = (mean(percent_R_OppFootOn)*100);

%% SINGLE LIMB SUPPORT (seconds)
% Single limb support is the time between opposite foot off and opposite foot on
% duration of Gait Cycles
for g = 1:length(adj_LGaitCycles_t(1,1,:))
    time_LGaitCycle(g) = (adj_LGaitCycles_t(1,2,g) - adj_LGaitCycles_t(1,1,g))*(1/120);  % seconds
end
for g = 1:length(adj_RGaitCycles_t(1,1,:))
    time_RGaitCycle(g) = (adj_RGaitCycles_t(1,2,g) - adj_RGaitCycles_t(1,1,g))*(1/120);  % seconds
end
% Adjust LFO_times and RFO_times
for n = 1:length(LFO_times)
    adj_LFO_times(n) = LFO_times(n) - FirstFrame;
end
for n = 1:length(RFO_times)
    adj_RFO_times(n) = RFO_times(n) - FirstFrame;
end
%%%%% Left %%%%%
% if first right foot on occurs before first right foot off
% ignore first right foot on
if adj_RFO_times(1) > adj_RFS_times(1)
    SS_RFS_times = adj_RFS_times(2:end);
else
    SS_RFS_times = adj_RFS_times;
end

% if last right foot off occurs after last right foot on
% ignore last right foot off
if adj_RFO_times(end) > adj_RFS_times(end)
    SS_RFO_times = adj_RFO_times(1:end-1);
else
    SS_RFO_times = adj_RFO_times;
end
for s = 1:length(SS_RFO_times)
    LeftSS(s) = abs(SS_RFO_times(s)-SS_RFS_times(s));  % frames
end
LSS_sec = mean(LeftSS*(1/120));  % seconds
LSingleSupport = 100 - ((abs(time_LGaitCycle-LSS_sec)/time_LGaitCycle)*100);  % as percent gait cycle
Std_LSingleSupport = std(100 - ((abs(mean(time_LGaitCycle) - (LeftSS * (1/120))) / (mean(time_LGaitCycle))) * 100));

%%%%% Right %%%%%
% if first left foot on occurs before first left foot off
% ignore first left foot on
if adj_LFO_times(1) > adj_LFS_times(1)
    SS_LFS_times = adj_LFS_times(2:end);
else
    SS_LFS_times = adj_LFS_times;
end

% if last left foot off occurs after last left foot on
% ignore last left foot off
if adj_LFO_times(end) > adj_LFS_times(end)
    SS_LFO_times = adj_LFO_times(1:end-1);
else
    SS_LFO_times = adj_LFO_times;
end
for s = 1:length(SS_LFO_times)
    RightSS(s) = abs(SS_LFO_times(s)-SS_LFS_times(s));  % frames
end
RSS_sec = mean(RightSS*(1/120));  % seconds
RSingleSupport = 100 - ((abs(time_RGaitCycle-RSS_sec)/time_RGaitCycle)*100);  % as percent gait cycle
Std_RSingleSupport = std(100 - ((abs(mean(time_RGaitCycle) - (RightSS * (1/120))) / (mean(time_RGaitCycle))) * 100));

%% DOUBLE LIMB SUPPORT (seconds)

%%%%% LEFT INITIAL DOUBLE SUPPORT %%%%%
if adj_RFO_times(1) < adj_LFS_times(1)
    % if first right foot off occurs before first left foot strike
    % ignore first right foot off
    DS_RFO_times = adj_RFO_times(2:end);
else
    DS_RFO_times = adj_RFO_times;
end
if DS_RFO_times(end) < adj_LFS_times(end)
    % if last left foot strike occurs after last right foot off
    % ignore last left foot strike
    DS_LFS_times = adj_LFS_times(1:end-1);
else
    DS_LFS_times = adj_LFS_times;
end
if length(DS_LFS_times) > length(DS_RFO_times)
    if DS_RFO_times(1,2) > DS_LFS_times(1,1)
        % remove first DS_RFS_time
        DS_LFS_times = DS_LFS_times(2:end);        
    end
end
for d = 1:length(DS_LFS_times)
    LeftDS(d) = abs(DS_LFS_times(d)-DS_RFO_times(d));  % frames
end
Init_LDS = mean(LeftDS*(1/120));  % seconds
Init_LDoubleSupport = (Init_LDS/mean(time_LGaitCycle)*100);  % as percent gait cycle
Std_Init_LDoubleSupport = std(LeftDS*(1/120)*mean(time_LGaitCycle))*100;

%%%%% RIGHT INITIAL DOUBLE SUPPORT %%%%%
if adj_LFO_times(1) < adj_LFS_times(1)
    % if first left foot off occurs before first right foot strike
    % ignore first left foot off
    DS_LFO_times = adj_LFO_times(2:end);
else
    DS_LFO_times = adj_LFO_times;
end
if DS_LFO_times(end) < adj_RFS_times(end)
    % if last right foot strike occurs after last left foot off
    % ignore last right foot strike
    DS_RFS_times = adj_RFS_times(1:end-1);
else
    DS_RFS_times = adj_RFS_times;
end
if length(DS_RFS_times) > length(DS_LFO_times)
    if DS_LFO_times(1,2) > DS_RFS_times(1,1)
        % remove first DS_RFS_time
        DS_RFS_times = DS_RFS_times(2:end);
    end
end
for d = 1:length(DS_LFO_times)
    RightDS(d) = abs(DS_RFS_times(d)-DS_LFO_times(d));  % frames
end
Init_RDS = mean(RightDS*(1/120));  % seconds
Init_RDoubleSupport = (Init_RDS/mean(time_RGaitCycle)*100);  % as percent gait cycle
Std_Init_RDoubleSupport = std(RightDS*(1/120)*mean(time_RGaitCycle))*100;

%%%%% LEFT FINAL %%%%%
% From opposite (right) foot strike to ipislateral (left) foot off
for f = 1:length(percent_Lfoot_offs)
    Final_LDS(f) = percent_Lfoot_offs(f) - percent_L_OppFootOn(f);
end
Final_LDoubleSupport = mean(Final_LDS)*100;  % percent gait cycle
Std_Final_LDoubleSupport = std(Final_LDS)*100; 

%%%%% RIGHT FINAL %%%%%
% From opposite (left) foot strike to ipislateral (right) foot off
for f = 1:length(percent_Rfoot_offs)
    Final_RDS(f) = percent_Rfoot_offs(f) - percent_R_OppFootOn(f);
end
Final_RDoubleSupport = mean(Final_RDS)*100;  % percent gait cycle
Std_Final_RDoubleSupport = std(Final_RDS)*100; 

%% Organizes Variables into Structure
% Overall
TempSpat(1).name = 'Cadence';
TempSpat(1).data = Cadence;
TempSpat(2).name = 'Walking_Speed';
TempSpat(2).data = Walking_Speed;
TempSpat(3).name = 'Ave_Stride_Length';
TempSpat(3).data = Ave_Stride_Length;

% Left
TempSpat(4).name = 'Ave_LStride_Time';
TempSpat(4).data = Ave_LStride_Time;
TempSpat(5).name = 'Ave_LStep_Time';
TempSpat(5).data = Ave_LStep_Time;
TempSpat(6).name = 'Ave_LStep_Length';
TempSpat(6).data = Ave_LStep_Length;
TempSpat(7).name = 'Ave_L_OppFootOff';
TempSpat(7).data = Ave_Percent_L_OppFootOff;
TempSpat(8).name = 'Ave_L_OppFootOn';
TempSpat(8).data = Ave_Percent_L_OppFootOn;
TempSpat(9).name = 'Ave_Percent_LFootOff';
TempSpat(9).data = Ave_Percent_LFootOff;
TempSpat(10).name = 'Init_LDoubleSupport';
TempSpat(10).data = Init_LDoubleSupport;
TempSpat(11).name = 'LSingleSupport';
TempSpat(11).data = LSingleSupport;
TempSpat(12).name = 'Final_LDoubleSupport';
TempSpat(12).data = Final_LDoubleSupport;

% Right
TempSpat(13).name = 'Ave_RStride_Time';
TempSpat(13).data = Ave_RStride_Time;
TempSpat(14).name = 'Ave_RStep_Time';
TempSpat(14).data = Ave_RStep_Time;
TempSpat(15).name = 'Ave_RStep_Length';
TempSpat(15).data = Ave_RStep_Length;
TempSpat(16).name = 'Ave_R_OppFootOff';
TempSpat(16).data = Ave_Percent_R_OppFootOff;
TempSpat(17).name = 'Ave_R_OppFootOn';
TempSpat(17).data = Ave_Percent_R_OppFootOn;
TempSpat(18).name = 'Ave_Percent_RFootOff';
TempSpat(18).data = Ave_Percent_RFootOff;
TempSpat(19).name = 'Init_RDoubleSupport';
TempSpat(19).data = Init_RDoubleSupport;
TempSpat(20).name = 'RSingleSupport';
TempSpat(20).data = RSingleSupport;
TempSpat(21).name = 'Final_RDoubleSupport';
TempSpat(21).data = Final_RDoubleSupport;

TempSpat(22).name = 'Ave_Step_Width';
TempSpat(22).data = Ave_Step_Width;

% STDs
TempSpat(23).name = 'STD_Step_Width';
TempSpat(23).data = Std_Step_Width;
% Left
TempSpat(24).name = 'STD_LStep_Length';
TempSpat(24).data = Std_LStep_Length;
TempSpat(25).name ='STD_Percent_LFootOff';
TempSpat(25).data = Std_Percent_LFootOff;
TempSpat(26).name ='STD_Init_LDoubleSupport';
TempSpat(26).data = Std_Init_LDoubleSupport;
TempSpat(27).name = 'STD_LSingleSupport';
TempSpat(27).data = Std_LSingleSupport;
TempSpat(28).name = 'STD_Final_LDoubleSupport';
TempSpat(28).data = Std_Final_LDoubleSupport;
% Right
TempSpat(29).name = 'STD_RStep_Length';
TempSpat(29).data = Std_RStep_Length;
TempSpat(30).name ='STD_Percent_RFootOff';
TempSpat(30).data = Std_Percent_RFootOff;
TempSpat(31).name ='STD_Init_RDoubleSupport';
TempSpat(31).data = Std_Init_RDoubleSupport;
TempSpat(32).name = 'STD_RSingleSupport';
TempSpat(32).data = Std_RSingleSupport;
TempSpat(33).name = 'STD_Final_RDoubleSupport';
TempSpat(33).data = Std_Final_RDoubleSupport;

% number of cycles on each side
TempSpat(34).name = 'N_L_Gait_Cycles';
TempSpat(34).data = N_LGaitCycles;
TempSpat(35).name = 'N_L_Gait_Cycles';
TempSpat(35).data = N_LGaitCycles;

end


%% Normal Age Matched Temporal Spatial Values
% norm_tempspat = [1.5yr    2yr     2.5yr   3yr     4yr     5yr     6yr      7yr     8yr     9yr     10yr    11yr    12yr    13yr    14yr    15yr   16yr  Adult]
% norm_tempspat =   [180      160.8	158.64	153.6	154.8	152.4	150     150     140.16	132.96	133.92	126.96	123     117     115.2	114     110.4	113.52; %Cadence
%                   0.68      0.743	0.758	0.777	0.773	0.79	0.82	0.835	0.858	0.9022	0.896	0.945	0.9756	1.027	1.035	1.036	1.029	1.057; %StrideTime
%                   17.5      16.7	15.5	15.5	14.15	13.45	13.4	12.35	12.6	11.9	11      11      13      13      13      14      14      13.3; %OppFootOff
%                   49.53     50.1	50.1	50.3	50      50      50      50      50      50.49	50.49	50.49	50.423	50.4	50.4	50.4	50.4	50.2; %OffFootOn
%                   0.34      0.367	0.38	0.39	0.386	0.39	0.41	0.415	0.43	0.4511	0.45	0.4724	0.4878	0.51	0.5172	0.5128	0.5145	0.5285; %StepTime
%                   32.3  	33.4	34.5	34.8	36      36.8	36.5	37.7	37.9	37.9	37.9	37.5	37.2	37.2	37.2	37.2	37.2	38; %SingleSupport
%                   35.2      33.6	31.2	30.2	27.6	26.9	27.2	25.2	25.12	25.1	25.1	25      25.1	25.1	25.1	25.1	25.1	24; %DoubleSupport
%                   67.6      66.8	65.7	65.5	63.8	63.6	63.6	62.4	62.4	61      61      62      65.48	65.5	65.5	65.5	65.5	65.5; %FootOff
%                   495       555     618     670     779     848     900     1050	1087	1150	1170	1170	1240	1350	1295	1310	1350	1320; %StrideLength
%                   247.54	277     309     335     390     424     450     525     543     530     585     585     622     675     647 	655     675     660; %StepLength
%                   44.58    44.52	49.02	51.48	60.36	64.62	67.5	78.78	76.08	70.44	78.36	74.28	76.5	78.96	74.52	74.7	74.52	74.94]; %WalkSpeed (m/min)
%                   %74.3      74.2	81.7	85.8	100.6	107.7	112.5	131.3	126.8	117.4	130.6	123.8	127.5	131.6	124.2	124.5	124.2	124.9]; %WalkSpeed



              
% %% Bar Graph
% 
% % Swing Time
% LSwingTime = round(100 - Ave_Percent_LFootOff);  % percent gait cycle
% RSwingTime = round(100 - Ave_Percent_RFootOff);  % percent gait cycle
% 
% % Temporal variables in bar graph
% Tvars = [round([TempSpat(10).data 12 TempSpat(19).data]);  % Initial Double Support
%           round([TempSpat(11).data 38 TempSpat(20).data]);  % Single Support
%           round([TempSpat(8).data 50 TempSpat(17).data]);  % Opposite Foot Off
%           round([TempSpat(12).data 12 TempSpat(21).data]);  % Final Double Support
%           round([LSwingTime 38 RSwingTime])];  % Swing time
% 
% % figure(1)
% % % Plot data in vertical bar graph
% b = bar(Tvars);
% % 
% hold on
% % xlabel('Gait Cycle Phases')
% % ylabel('Percent Gait Cycle (0-100%)')
% % title('Temporal Gait Cycle Variables')
% % 
% % % Set bar graph y axis limits
% % set(gca, 'ylim', [0 100])
% % 
% % 
% % % Label bar groups
% % set(gca, 'XTickLabel', {'IDS', 'SS', 'OFO', 'FDS', 'SW'})
% % 
% % % Bar colors
% ch = get(b,'children');
% set(ch{1},'FaceColor', [0.8 0 0])  % sets left values to red
% set(ch{2},'FaceColor', [0.7 0.7 0.7])  % sets normal values to grey
% set(ch{3},'FaceColor', [0 0.6 0.2])  % sets right values to green
% 
% 
% % figure(2)
% % % % horizontal bar graph
% % TSvarsh = [round([TempSpat(19).data TempSpat(20).data TempSpat(17).data TempSpat(21).data RSwingTime]); 
% %            12 38 50 12 38;
% %            round([TempSpat(10).data TempSpat(11).data TempSpat(8).data TempSpat(12).data TempSpat(12).data])];
% %            
% % 
% % % Plot data in horizontal bar graph
% % b = barh(TSvars', 'stack');
% % 
% % 
% % 
% % hold on
% % xlabel('Percent Gait Cycle (0-100%)')
% % ylabel('Side')
% % title('Temporal Spatial Gait Cycle Variables')
% % grid on
% % 
% % % Set bar graph y axis limits
% % % set(gca, 'xlim', [0 100])
% % 
% % % Label bar groups
% % set(gca, 'YTickLabel', {'Left', 'Normal', 'Right'})
% % 
% % % Bar colors
% % hch = get(b,'children');
% % set(hch{1},'FaceColor', [0.8 0 0])  % sets left values to red
% % set(hch{2},'FaceColor', [0.7 0.7 0.7])  % sets normal values to grey
% % set(hch{3},'FaceColor', [0 0.6 0.2])  % sets right values to green
% % set(hch{3},'FaceColor', [0 0 0.8])  % sets right values to green
% % set(hch{3},'FaceColor', [0 0.6 0.2])  % sets right values to green
% 
% 
% %% Tempora Spatial Variables Bar Graph
% 
% % Spatial variables in bar graph
% % TSvars = [TempSpat(4).data norm_tempspat(2,18) TempSpat(13).data;  % Stride Time (sec)
% %          TempSpat(5).data norm_tempspat(5,18) TempSpat(14).data;  % Step Time (sec)
% %          TempSpat(6).data (norm_tempspat(10,18)/10) TempSpat(15).data];  % Step Length (cm)
% 
% 
% % Expressed as percentage of normal
% nStpL = norm_tempspat(10,18)/10; %cm
% 
% LStrT = (TempSpat(4).data/norm_tempspat(2,18)) *100;
% LStpT = (TempSpat(5).data/norm_tempspat(5,18)) *100;
% LStpL = ((nStpL - abs(TempSpat(6).data-nStpL))/nStpL)*100;
% 
% RStrT = (norm_tempspat(2,18) - abs(TempSpat(13).data-norm_tempspat(2,18)))*100;
% RStpT = (norm_tempspat(5,18) - abs(TempSpat(14).data-norm_tempspat(5,18)))*100;
% RStpL = ((nStpL - abs(TempSpat(15).data-nStpL))/nStpL)*100;
% 
% TSvars = [LStrT 100 RStrT;  % Stride Time (sec)
%           LStpT 100 RStpT;  % Step Time (sec)
%           LStpL 100 RStpL];  % Step Length (cm)
% 
% hold on
% % figure(2)
% % Plot data in vertical bar graph
% b2 = bar(TSvars);
% 
% hold on
% xlabel('Parameter')
% ylabel('Percentage of Normal')
% title('Temporal Spatial Variables')
% 
% % Set bar graph y axis limits
% set(gca, 'ylim', [0 120])
% 
% % add the text to the plot
% % text(TSvars,TSvars, sL);
% 
% % Label bar groups
% set(gca, 'XTickLabel', {'Stride Time', 'Step Time', 'Step Length'})
% 
% % Bar colors
% ch2 = get(b2,'children');
% set(ch2{1},'FaceColor', [0.8 0 0])  % sets left values to red
% set(ch2{2},'FaceColor', [0.7 0.7 0.7])  % sets normal values to grey
% set(ch2{3},'FaceColor', [0 0.6 0.2])  % sets right values to green
% 
% [Lmax, Lind] = sort(TSvars(:,1));
% [Nmax, Nind] = sort(TSvars(:,2));
% [Rmax, Rind] = sort(TSvars(:,3));
% plot(Lind,Lmax, 'r.')
% plot(Nind,Nmax, 'k.')
% plot(Rind,Rmax, 'g.')
% hold off
% 
% Ltxt = text(Lind,Lmax,{round(TSvars(:,1))});
% set(Ltxt, 'Rotation', 90)
% Ntxt = text(Nind,Nmax,{round(TSvars(:,2))});
% set(Ntxt, 'Rotation', 90)
% Rtxt = text(Rind,Rmax,{round(TSvars(:,3))});
% set(Rtxt, 'Rotation', 90)
% 
% % Text outside of plot axes
% % str(1) = {'Plot of the function:'};
% % str(2) = {' y = A{\ite}^{-\alpha{\itt}}'};
% % str(3) = {'With the values:'};
% % str(3) = {' A = 0.25'};
% % str(4) = {' \alpha = .005'};
% % str(5) = {' t = 0:900'};
% % set(gcf,'CurrentAxes')
% % text(1,0,str,'FontSize',8)
% 
