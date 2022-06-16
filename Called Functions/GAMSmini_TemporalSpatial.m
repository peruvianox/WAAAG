function [TempSpat] = GAMSmini_TemporalSpatial(HeaderInfo, Point_xyz, PointLabels, GaitEvents)
% -------------------------------------------------------------------------
% TEMPORAL SPATIAL CALCULATIONS
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: Setpember 16, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
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

%% DEBUG LOOP
close all
clear all
clc

% Prompts User for Input
[FullFileName, pathname] = uigetfile('*.c3d*', 'Please select C3D file to be analyzed');

[C3Dfile] = OpenC3D(FullFileName);

% Open and Read the user specified C3D file
% and returns: HeaderInfo, Point_xyz, PointLabels
[HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);

%% Data Extraction
% Extract left heel marker
pt_label = 'LHEE';
[pt_xyz] = searchPoint(Point_xyz, PointLabels, pt_label);
LHEE = pt_xyz;
% Extract right heel marker
pt_label = 'RHEE';
[pt_xyz] = searchPoint(Point_xyz, PointLabels, pt_label);
RHEE = pt_xyz;

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
for n = 1:length(LFS_times)
    adj_LFS_times(n) = LFS_times(n) - FirstFrame;
end

% Number of left gait cycles in trial
N_LGaitCycles = length(LFS_times)-1;
if N_LGaitCycles ~= 0
    for numGC = 1:N_LGaitCycles
        start_LGC(numGC) = adj_LFS_times(numGC);

        for numGC_end = 1:N_LGaitCycles;
            end_LGC(numGC_end) = adj_LFS_times(numGC_end+1);
        end
        % adjusted left gait cycle times (foot strikes)
        adj_LGaitCycles_t(:,:,numGC) = [start_LGC(numGC), end_LGC(numGC)];

    end
end

for n = 1:length(RFS_times)
    adj_RFS_times(n) = RFS_times(n) - FirstFrame;
end

% Number of right gait cycles in trial
N_RGaitCycles = length(RFS_times)-1;
if N_RGaitCycles ~= 0
    for numGC = 1:N_RGaitCycles
        start_RGC(numGC) = adj_RFS_times(numGC);

        for numGC_end = 1:N_RGaitCycles;
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
% where total foot strikes - 1 because 1st foot strike indicates begining of measurement
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
Ave_Stride_Length = abs(sum(stride_length)/length(stride_length));


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


%% STEP TIME (seconds) & LENGTH (centimeters)

%%%%% LEFT %%%%%
% duration of time from right foot strike to next left foot strike
if adj_LFS_times(1) < adj_RFS_times(1)
    % Step Time
    % if left foot strike occurs first
    % ignore first left foot strike and count from first right foot strike
    new_LFS_times = adj_LFS_times(2:end);

    for n = 1:length(new_LFS_times)
        % times for each left step (seconds)
        Lstep_times(n) = abs(adj_RFS_times(n) - new_LFS_times(n))*(1/120);
    end
    
    % if uneven number of left and right foot strike times, adjust times
    if length(adj_RFS_times) < length(new_LFS_times)
        % if more left foot strike times than right, remove last left
        new_LFS_times = new_LFS_times(1:end-1);
    end
    
    
    % Step Length
    if length(adj_RFS_times) > length(new_LFS_times)
        % if there are more right foot strikes than (new adjusted) left
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
    for n = 1:length(adj_LFS_times)
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

        
    % Step Length
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
    for n = 1:length(adj_RFS_times)
        % calculate times for each right step (seconds)
        Rstep_times(n) = abs(adj_LFS_times(n) - adj_RFS_times(n))*(1/120);
    end
    
    
    % Step Length
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
Ave_Percent_LFootOff = (mean(percent_Lfoot_offs(:,t))*100);
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
Ave_Percent_RFootOff = (mean(percent_Rfoot_offs(:,t))*100);
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
elseif length(new_LFO_times) == length(LFS_times)
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

%%%%% LEFT INITIAL %%%%%
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

for d = 1:length(DS_RFO_times)
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

TempSpat(34).name = 'Num_Left_Cycles';
TempSpat(34).data = length(LGaitCycles);
TempSpat(35).name = 'Num_Right_Cycles';
TempSpat(35).data = length(RGaitCycles);


end
