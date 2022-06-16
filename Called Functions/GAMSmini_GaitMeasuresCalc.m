function [GaitMeasures] = GAMSmini_GaitMeasuresCalc(FullFileName, age)
% -------------------------------------------------------------------------
% GAMS MINI PROGRAM - BATCH PROCESSING OF C3D FILES
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: April 22, 2010
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This function batch processes all the c3d files and
%               determines most representative trial using each trial's kinematics
%               and temporal spatial variables.
% 
% Input:        files2get        user selected C3D file(s) to be analyzed
%               
%
% Output:       Representative  structure with representative trial's
%               kinematics, temporal spatial data and trial name
%              

%% DEBUG LOOP
% close all
% clear all
% clc
% 
% [FullFileName, pathname] = uigetfile('*.c3d*', 'Please select all C3D file(s) to be analyzed',...
%                                        'MultiSelect', 'off');
%                                    
% age = 10;


%% Opens and Extract C3D file's data
[C3Dfile] = OpenC3D(FullFileName);

% Open and Read the user specified C3D file
% and returns: HeaderInfo, Point_xyz, PointLabels
[HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);
[TempSpat] = TemporalSpatial(HeaderInfo, Point_xyz, PointLabels, GaitEvents);
                       
%% Overall Gait Measures    

% Stance Percent of Gait Cycle
Ave_Percent_LFootOff = TempSpat(9).data;
LStanceGC = round(Ave_Percent_LFootOff);
Ave_Percent_RFootOff = TempSpat(18).data;
RStanceGC = round(Ave_Percent_RFootOff);

% Swing Percent of Gait Cycle
LSwingGC = round(100 - Ave_Percent_LFootOff);  % percent gait cycle
RSwingGC = round(100 - Ave_Percent_RFootOff);  % percent gait cycle

% Cadence and Walking Speed
Cadence = round(TempSpat(1).data); 
Walk_Speed = round(TempSpat(2).data);


% Overall Gait Performance variables matrix
% as a percentage of normal/reference
% Cadence_n = norm_tempspat(1,age); % age matched cadence
% Cadence_p = (Cadence/Cadence_n)*100;
% 
% Walk_Speed_n = norm_tempspat(11,age); % age matched walking speed
% Walk_Speed_p = (Walk_Speed/Walk_Speed_n)*100;


% 
% OGP_vars = [Cadence_p 0;
%             0         Walk_Speed_p];
% 
% % store in GaitMeasures structure
% GaitMeasures.OGP_vars = OGP_vars;
    
    
%% SPATIAL MEASURES
% Stride length (m)
%nStrdL = norm_tempspat(9,age)/1000;  % Normal Stride Length
Ave_Stride_Length = (TempSpat(3).data);  % Trial's Average Stride Length

% Expressed as percentage of normal/reference
% if Ave_Stride_Length > nStrdL
%     Stride_Length = ((abs(Ave_Stride_Length-nStrdL)/nStrdL)*100)+100;
% else
%     Stride_Length = 100-((abs(Ave_Stride_Length-nStrdL)/nStrdL)*100);
% end
% Stride_Length = (Ave_Stride_Length/nStrdL)*100;

% Step Length (cm)
%nStpL = norm_tempspat(10,age)/10;  % Normal Step Length
Left_StpL = TempSpat(6).data;
% Expressed as percentage of normal/reference
% if Left_StpL > nStpL
%     LStep_Length = ((abs(Left_StpL-nStpL)/nStpL)*100)+100;
% else
%     LStep_Length = 100-((abs(Left_StpL-nStpL)/nStpL)*100);
% end
% LStep_Length = (Left_StpL/nStpL)*100;

Right_StpL = TempSpat(15).data;
% Expressed as percentage of normal/reference
% if Right_StpL > nStpL
%     RStep_Length = ((abs(Right_StpL-nStpL)/nStpL)*100)+100;
% else
%     RStep_Length = 100-((abs(Right_StpL-nStpL)/nStpL)*100);
% end
% RStep_Length = (Right_StpL/nStpL)*100;

% Spatial variables matrix (expressed as percent of normal)
% Spatial_vars = [Stride_Length  0             0;
%                 0              LStep_Length  RStep_Length];

% store in GaitMeasures structure
% GaitMeasures.Spatial_vars = Spatial_vars;
       

%% TEMPORAL MEASURES

% Temporal variables matrix, as a percent of gait cycle
Temporal_vars = [LStanceGC         62 RStanceGC;  % Stance Period
                 LSwingGC          38 RSwingGC;  % Swing Period
                 TempSpat(10).data 12 TempSpat(19).data;  % Initial Double Support
                 TempSpat(11).data 38 TempSpat(20).data;  % Single Limb Support
                 TempSpat(12).data 12 TempSpat(21).data];  % Final Double Support
    
% store in GaitMeasures structure
GaitMeasures.Temporal_vars = Temporal_vars;


%% TABLE DATA
StanceGC = round((round(LStanceGC) + round(RStanceGC))/2);
SwingGC = round((round(LSwingGC) + round(RSwingGC))/2);
IDS = round((round(TempSpat(10).data) + round(TempSpat(19).data))/2);
SLS = round((round(TempSpat(11).data) + round(TempSpat(20).data))/2);
FDS = round((round(TempSpat(12).data) + round(TempSpat(21).data))/2);

% values for table of gait measures
% TempSpat_table_data = {'steps/min'          'm'                         'm/min'              'cm'                        'cm'                       '% Gait Cycle'   '% Gait Cycle'    '% Gait Cycle'           '% Gait Cycle'           '% Gait Cycle';
%                        Cadence              TempSpat(3).data            Walk_Speed            TempSpat(6).data            TempSpat(15).data          StanceGC        SwingGC            IDS                      SLS                      FDS;
%                        norm_tempspat(1,age) (norm_tempspat(9,age)/1000) norm_tempspat(11,age) (norm_tempspat(10,age)/10)  (norm_tempspat(10,age)/10) 62              38                 12                       38                       12};
%     
% store in GaitMeasures structure
% GaitMeasures.TempSpat_table_data = TempSpat_table_data;


%% All Gait Parameters

% GaitParams = [Cadence            0                  0      0;
%               Walk_Speed         0                  0      0;
%               Ave_Stride_Length  0                  0      0;
%               0                  Left_StpL          nStpL  Right_StpL
%               0                  LStanceGC          62     RStanceGC;
%               0                  LSwingGC           38     RSwingGC;
%               0                  TempSpat(10).data  12     TempSpat(19).data;
%               0                  TempSpat(11).data  38     TempSpat(20).data;
%               0                  TempSpat(12).data  12     TempSpat(21).data];
% 
% GaitMeasures.GaitParams = GaitParams;

%% 
GaitMeasures.TempSpat = TempSpat; 

end





