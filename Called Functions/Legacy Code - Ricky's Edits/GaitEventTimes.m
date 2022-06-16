function [Events_times, LGaitCycles, RGaitCycles] = GaitEventTimes(GaitEvents)
% -------------------------------------------------------------------------
% GAIT EVENT TIMES
% -------------------------------------------------------------------------
% AUTHOR: Kate Worster
% DATE: April 13, 2010 (newest version)
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% DESCRIPTION:	This is function determines the times at which all
%               designated gait events (in c3d file) occur and organizes them by side
%               (left/right) into gait cycles.
%
% INPUT:        GaitEvents        Cell array of gait event information from
%                                 c3d file (events, side, times)
% 
%
% OUTPUT:  
%               LFS_times         Left foot strike times for respective side
%               LFO_times         Left foot off times for respective side
%               RFS_times         Right foot strike times for respective side
%               RFO_times         Right foot off times for respective side
%               LGaitCycles_t     Matrix of left gait cycle times
%               RGaitCycles_t     Matrix of right gait cycle times
%

% % DEBUG LOOP
% close all
% clear all
% clc
% 
% side = 'left';
% 
% [FullFileName, pathname] = uigetfile('*.c3d*', 'Please select a C3D file to be analyzed');
% 
% % Opens and Extract C3D file's data
% [C3Dfile] = OpenC3D(FullFileName);
% 
% % Open and Read the user specified C3D file
% % and returns: HeaderInfo, Point_xyz, PointLabels
% [HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);
% 
% %%%% END DEBUG LOOP %%%%%

%% Organize Events 

%%%%% FOOT STRIKE TIMES %%%%%
% Determines number of Left and Right Foot Strikes
% by searching for specific string in GaitEvents from C3D file
N_FS1 = strcmp(GaitEvents(:,:,1), 'Foot Strike');
N_Left = strcmp(GaitEvents(:,:,2), 'Left');
N_Right = strcmp(GaitEvents(:,:,2), 'Right');

% Calculates the number of gait cycles per side
LFootStrikes = 0;
RFootStrikes = 0;
m = 1;

for n = 1:length(GaitEvents(1,:,1))
    if (N_FS1(n) == 1) && (N_Left(n) == 1)
         LFootStrikes = LFootStrikes + 1;  
         LFootStrike_times(1,m) = GaitEvents(:,n,3);
         m = m+1;
    end

    if (N_FS1(n) == 1) && (N_Right(n) == 1)
         RFootStrikes = RFootStrikes + 1;
         RFootStrike_times(1,m) = GaitEvents(:,n,3);
         m = m+1;
    end
end

if exist('LFootStrike_times') == 1
    % Removes empty cells, converts from cell to matrix
    LFootStrike_times(cellfun(@isempty,LFootStrike_times)) = [];
    % stores the frames for when L foot strike occurs
    LFS_times = double(cell2mat(LFootStrike_times))+1;     
else
    LFS_times = 0;
end

if exist('RFootStrike_times') == 1
    % Removes empty cells, converts from cell to matrix
    RFootStrike_times(cellfun(@isempty,RFootStrike_times)) = [];
    % stores the frames for when R foot strike occurs
    RFS_times = double(cell2mat(RFootStrike_times))+1;
else
    RFS_times = 0;
end

%%%%% FOOT OFF & EVENT TIMES %%%%%
% Determines number of Left and Right Foot Off
% by searching for specific string in GaitEvents from C3D file
N_FO1 = strcmp(GaitEvents(:,:,1), 'Foot Off');
N_Event = strcmp(GaitEvents(:,:,1), 'Event');
N_Left_str = strcmp(GaitEvents(:,:,2), 'Left');
N_Right_str = strcmp(GaitEvents(:,:,2), 'Right');

% Calculates the number of gait cycles per side
LFootOffs = 0;
RFootOffs = 0;
LEvents = 0;
REvents = 0;
k = 1;

for n2 = 1:length(GaitEvents(1,:,1))
    % Looks for left 'Foot Off' events
    if (N_FO1(n2) == 1) && (N_Left_str(n2) == 1)
         LFootOffs = LFootOffs + 1;  
         LFootOff_times(1,k) = GaitEvents(1,n2,3);
         k = k+1;
    end
    
    % Looks for left general 'Event' events
    if (N_Event(n2) == 1) && (N_Left_str(n2) == 1)
        LEvents = LEvents + 1;  
        LEvents_times(1,k) = GaitEvents(1,n2,3);
        k = k+1;
    end
    
     % Looks for right 'Foot Off' events
    if (N_FO1(n2) == 1) && (N_Right_str(n2) == 1)
         RFootOffs = RFootOffs + 1;
         RFootOff_times(1,k) = GaitEvents(1,n2,3);
         k = k+1;
    end
    
    % Looks for right general 'Event' events
    if (N_Event(n2) == 1) && (N_Right_str(n2) == 1)
        REvents = REvents + 1;  
        REvents_times(1,k) = GaitEvents(1,n2,3);
        k = k+1;
    end
end

if exist('LFootOff_times') == 1
    % Removes empty cells, converts from cell to matrix
    LFootOff_times(cellfun(@isempty,LFootOff_times)) = [];
    % stores the frames for when L foot off occurs
    LFO_times = double(cell2mat(LFootOff_times))+1;   
else
    LFO_times = 0;
end

if exist('RFootOff_times') == 1
    % Removes empty cells, converts from cell to matrix
    RFootOff_times(cellfun(@isempty,RFootOff_times)) = [];
    % stores the frames for when R foot off occurs
    RFO_times = double(cell2mat(RFootOff_times))+1;
else
    RFO_times = 0;
end

if exist('LEvents_times') == 1
    % Removes empty cells, converts from cell to matrix
    LEvents_times(cellfun(@isempty,LEvents_times)) = [];
    % stores the frames for when L foot off occurs
    LEvent_times = double(cell2mat(LEvents_times))+1;   
else
    LEvent_times = 0;
end

if exist('REvents_times') == 1
    % Removes empty cells, converts from cell to matrix
    REvents_times(cellfun(@isempty,REvents_times)) = [];
    % stores the frames for when L foot off occurs
    REvent_times = double(cell2mat(REvents_times))+1;   
else
    REvent_times = 0;
end

% Stores results in structure
Events_times(1).name = 'LFS_times';
Events_times(1).data = LFS_times;
Events_times(2).name = 'LFO_times';
Events_times(2).data = LFO_times;
Events_times(3).name = 'LEvent_times';
Events_times(3).data = LEvent_times;
Events_times(4).name = 'RFS_times';
Events_times(4).data = RFS_times;
Events_times(5).name = 'RFO_times';
Events_times(5).data = RFO_times;
Events_times(6).name = 'REvent_times';
Events_times(6).data = REvent_times;

%% Gait Cycles

% Number of left gait cycles in trial
N_LGaitCycles = length(LFS_times)-1;


if N_LGaitCycles ~= 0
    for numGC = 1:N_LGaitCycles
        start_LGC(numGC) = LFS_times(numGC);

        for numGC_end = 1:N_LGaitCycles;
            end_LGC(numGC_end) = LFS_times(numGC_end+1);

        end

        LGaitCycles(:,:,numGC) = [start_LGC(numGC), end_LGC(numGC)];

    end

else
    % displays if no gait cycles are present for selected side
    disp('No left gait cycles present')
    LGaitCycles = 0;
end


% Number of right gait cycles in trial
N_RGaitCycles = length(RFS_times)-1;

if N_RGaitCycles ~= 0
    for numGC = 1:N_RGaitCycles
        start_RGC(numGC) = RFS_times(numGC);

        for numGC_end = 1:N_RGaitCycles;
            end_RGC(numGC_end) = RFS_times(numGC_end+1);

        end

        RGaitCycles(:,:,numGC) = [start_RGC(numGC), end_RGC(numGC)];

    end

else
    % displays if no gait cycles are present for selected side
    disp('No right gait cycles present')
    RGaitCycles = 0;
end
