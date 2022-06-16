function [Output] = GetControlKinematics (age)
% % -------------------------------------------------------------------------
% Obtains Control Kinematics 
% -------------------------------------------------------------------------
% Author: Ricky Pimentel
% Date: December 20, 2016
%
% Description:	This function will output a dataset of the mean and standard
%                       deviation of control kinematics for all lower body joints (unilateral) during normal
%                       walking. The N of subjects used for each control dataset vary by age
%                       group. The output data will be arranged in a 51x12 matrix, with 3 columns to
%                       each joint angle (Pelvis, Hip, Knee, Ankle). 
%
%                  * must have Control_Kinematics.mat file in order to function
%
% Input:       age       Age (to the nearest year) of desired kinetic controls 
%                               *if age is not input while called, it will
%                               be prompted for
%                               **general kinemaitics (unspecied ages)
%                               like those used for GDI, can be elected by
%                               inputting 0 for age
%
% Output:     Average and Standard deviations for joint angles
%                   of the desired age(s)

% General Kinetics?

% % DEBUG LOOP %%%%%%%%%
% close all
% clear
% clc
% %%% END DEBUG LOOP %%%%%



% General Kinematics can be chosen by inputting
% 0 as the age for:   R Baker's GDI and GPS controls -> http://dx.doi.org/10.1016/j.gaitpost.2008.05.001

% note - there are no Knee ABD/ADD, Knee Rot, or Ankle ABD/ADD for the general kinmatics dataset

%% Ask for age if no input
a = exist('age');
if a == 0
prompt = 'What is the age of the subject?';
age = inputdlg(prompt);
end

if iscell(age)
    age = cell2mat(age);
end

if isstr(age)
    age = str2double(age);
end

%% Load control kinetics data
load('Control_Kinematics.mat');

%% Find Data depending on age
if age == 4 % if they are 4 or 5
    CK = Control_Kinematics.Age4_5;
elseif age == 5
    CK = Control_Kinematics.Age4_5;
end

if age == 6 % if they are 6 or 7
    CK = Control_Kinematics.Age6_7;
elseif age == 7
    CK = Control_Kinematics.Age6_7;
end

if age == 8 % if they are 8 or 9
    CK = Control_Kinematics.Age8_9;
elseif age == 9
    CK = Control_Kinematics.Age8_9;
end

if age == 10 % if they are 10 or 11
    CK = Control_Kinematics.Age10_11;
elseif age == 11
    CK = Control_Kinematics.Age10_11;
end

if age == 12 % if they are 12 or 13
    CK = Control_Kinematics.Age12_13;
elseif age == 13
    CK = Control_Kinematics.Age12_13;
end

if age == 14 % if they are 14 or 15
    CK = Control_Kinematics.Age14_15;
elseif age == 15
    CK = Control_Kinematics.Age14_15;
end

if age == 16 % if they are 16 or 17
    CK = Control_Kinematics.Age16_17;
elseif age == 17
    CK = Control_Kinematics.Age16_17;
end

if age == 18 % if they are 18 or 19
    CK = Control_Kinematics.Age18_19;
elseif age ==19 
    CK = Control_Kinematics.Age18_19;
end

if age >= 20 && age <= 29 % if they are between 20 and 29
    CK = Control_Kinematics.Age20_29;
end

if age >= 30 && age <= 49 % if they are between 30 and 49
    CK = Control_Kinematics.Age30_49;
end

if age >= 50 && age <= 80 % if they are between 50 and 80
    CK = Control_Kinematics.Age50_80;
end

if age == 0 % if General Kinematics from R Baker's Data -> http://dx.doi.org/10.1016/j.gaitpost.2008.05.001
    CK = Control_Kinematics.General;
end

Output = CK;

end
    