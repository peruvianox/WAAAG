 function GDI = GDI_HemiInitial(file2get,Hemi)
% -------------------------------------------------------------------------
% GDI Calculator: Decides what side to select representative gait vector
% from and returns left, right, and average GDI score for that patient.
% -------------------------------------------------------------------------
%
% Modifying Author: Kayla Burnim
% Start Date: 09Jan2015
% Modified GDI_Prompt for use in GAMS
%   Author: Colton Sauer
%   Date: July 9, 2014
%
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	Determines what side to select representative gait cycle
%               from in order to calculate GDI and calculates it.
% 
% Input:        Hemi Side
%               
%
% Output:       GDI accounting for hemiplegia. GDI(1) is Left GDI, GDI(2)
%               is Right GDI, and GDI(3) is the average of the left and right sides.

%% DEBUG LOOP
% close all
% clear all
% clc
% %Prompt user to choose what subject to extract kinematic data for.
% [file2get, pathname] = uigetfile('*.c3d', 'Please select SUBJECT C3D file to be analyzed');
% Hemi='No';

%%%%% END DEBUG LOOP %%%%%%%%%%
warning off; 
% load control 
load('GDI_Calc_Fast_Data.mat'); 

%% Opens and Extract C3D file's data
[C3Dfile] = OpenC3D(file2get);

% Open and Read the user specified C3D file
% and returns: HeaderInfo, Point_xyz, PointLabels
[~, ~, ~, GaitEvents] = ReadC3D(C3Dfile);

%Find when first left/right gait cycle starts
[~, LGaitCycles, RGaitCycles] = GaitEventTimes(GaitEvents);
LCycleStart = LGaitCycles(1,1,1);
RCycleStart = RGaitCycles(1,1,1);

%Find representative trial for each side
[~, residuals, ~, ~, ~, ~] = GDI_AnalyTrial(file2get);
L_Rep = residuals.all_GCs.left(1);
R_Rep = residuals.all_GCs.right(1);

%% User Input and Calculating GDI
%Prompt to determine if patient is hemiplegic and which side is affected to
%find representative trial.
if strcmp(Hemi, 'Right') == 1
    % If the first cycle occurs on the left foot
    if LCycleStart < RCycleStart
        %If the representative cycle is the first cycle on the right
        %foot, or is any cycle before the last recorded cycle for the
        %right foot.
        if R_Rep == 1 || R_Rep < length(RGaitCycles)
            L_Rep = R_Rep + 1;
            %If the representative cycle occurs on the last cycle for the
            %right foot.
        else
            %If there are the same amount of L and R gait cycles
            if length(LGaitCycles) == length(RGaitCycles)
                L_Rep = R_Rep;
                %If there are more left gait cycles than right.
            else
                L_Rep = R_Rep + 1;
            end
        end
        %If the first cycle occurs on the right foot
    else
        %If the representative cycle is the first cycle on the right
        %foot, or is any cycle before the last recorded cycle for the
        %right foot.
        if R_Rep == 1 || R_Rep < length(RGaitCycles)
            L_Rep = R_Rep;
            %If the representative cycle occurs on the last cycle for the
            %right foot.
        else
            %If there are the same amount of L and R gait cycles
            if length(LGaitCycles) == length(RGaitCycles)
                L_Rep = R_Rep;
                %If there are more right gait cycles than left.
            else
                L_Rep = R_Rep - 1;
            end
        end
    end
    GDI = GDI_Calculator(L_Rep, R_Rep, file2get);
elseif strcmp(Hemi, 'Left') == 1
    %If the first cycle occurs on the left foot
    if LCycleStart < RCycleStart
        %If the representative cycle is the first cycle on the left
        %foot, or if any cycle before the last recorded cycle for the
        %left foot.
        if L_Rep == 1 || L_Rep < length(LGaitCycles)
            R_Rep = L_Rep;
            %If the the reperesentative cycle occurs on the last cycle of
            %the left foot.
        else
            %If there are the same amount of L and R gait cycles
            if length(LGaitCycles) == length(RGaitCycles)
                R_Rep = L_Rep;
                %If there are more left gait cycles than right
            else
                R_Rep = L_Rep - 1;
            end
        end
        %If the first cycle occurs on the right foot
    else
        %If the representative cycle is the first cycle on the left
        %foot or any cycle before the last recorded cycle for the left
        %foot.
        if L_Rep == 1 || L_Rep < length(LGaitCycles)
            R_Rep = L_Rep + 1;
            %If the representative cycle occurs on the last cycle of the
            %right foot.
        else
            %If the left and right sides have the same amount of cycles
            if length(LGaitCycles) == length(RGaitCycles)
                R_Rep = L_Rep;
                %If there are more right gait cycles than left
            else
                R_Rep = L_Rep + 1;
            end
        end
    end
    GDI = GDI_Calculator(L_Rep, R_Rep, file2get);
elseif strcmp(Hemi, 'No') == 1
    %Find which side has the lowest GDI
    GDI = GDI_Calculator_Fast(L_Rep, R_Rep, file2get, 0);
    %If left GDI is less than the right use left representative
    if GDI(1) < GDI(2)
        clear GDI
        if LCycleStart < RCycleStart
            %If the representative cycle is the first cycle on the left
            %foot, or if any cycle before the last recorded cycle for the
            %left foot.
            if L_Rep == 1 || L_Rep < length(LGaitCycles)
                R_Rep = L_Rep;
                %If the the reperesentative cycle occurs on the last cycle of
                %the left foot.
            else
                %If there are the same amount of L and R gait cycles
                if length(LGaitCycles) == length(RGaitCycles)
                    R_Rep = L_Rep;
                    %If there are more left gait cycles than right
                else
                    R_Rep = L_Rep - 1;
                end
            end
            %If the first cycle occurs on the right foot
        else
            %If the representative cycle is the first cycle on the right
            %foot or any cycle before the last recorded cycle for the left
            %foot.
            if L_Rep == 1 || L_Rep < length(LGaitCycles)
                R_Rep = L_Rep + 1;
                %If the representative cycle occurs on the last cycle of the
                %right foot.
            else
                %If the left and right sides have the same amount of cycles
                if length(LGaitCycles) == length(RGaitCycles)
                    R_Rep = L_Rep;
                    %If there are more right gait cycles than left
                else
                    R_Rep = L_Rep + 1;
                end
            end
        end
        GDI = GDI_Calculator_Fast(L_Rep, R_Rep, file2get,0);
        %If right GDI is lower than the left
    else
        % If the first cycle occurs on the left foot
        if LCycleStart < RCycleStart
            %If the representative cycle is the first cycle on the right
            %foot, or is any cycle before the last recorded cycle for the
            %right foot.
            if R_Rep == 1 || R_Rep < length(RGaitCycles)
                L_Rep = R_Rep + 1;
                %If the representative cycle occurs on the last cycle for the
                %right foot.
            else
                %If there are the same amount of L and R gait cycles
                if length(LGaitCycles) == length(RGaitCycles)
                    L_Rep = R_Rep;
                    %If there are more left gait cycles than right.
                else
                    L_Rep = R_Rep + 1;
                end
            end
            %If the first cycle occurs on the right foot
        else
            %If the representative cycle is the first cycle on the right
            %foot, or is any cycle before the last recorded cycle for the
            %right foot.
            if R_Rep == 1 || R_Rep < length(RGaitCycles)
                L_Rep = R_Rep;
                %If the representative cycle occurs on the last cycle for the
                %right foot.
            else
                %If there are the same amount of L and R gait cycles
                if length(LGaitCycles) == length(RGaitCycles)
                    L_Rep = R_Rep;
                    %If there are more right gait cycles than left.
                else
                    L_Rep = R_Rep - 1;
                end
            end
        end
        GDI = GDI_Calculator_Fast(L_Rep, R_Rep, file2get,0);
    end
end