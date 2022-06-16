function GDI= GDI_Calculator_Fast(L_Rep, R_Rep, files2get, KineData)
% -------------------------------------------------------------------------
%  Calculates GDI
% -------------------------------------------------------------------------
% Author: Colton Sauer
% Date: June 7th, 2014
% Edited by Ricky Pimentel to improve speed
%
% DESCRIPTION: Calculates GDI from Gillette control population.
%
% INPUT: Matrix of gait vectors for every control subject and feature
% components in excel spreadsheet, representative cycles for both sides.      
%
% OUTPUT:  GDI for a certain patient outlined by Schwartz' GDI article. Values are in the following order: [Left_GDI, Right_GDI, Average_GDI, GDI_Average_Vector] 
%   The Average GDI in the 3 position is calculated from average of the
%   left and right GDI's. GDI_Average_Vector is calculated from averaging
%   the left and right GAIT VECTORS for the subject and then GDI is
%   calculated from that single vector.  

% % DEBUG LOOP
% close all
% clear all
% clc
% %%%% END DEBUG LOOP %%%%%

warning('off');
%% Quickly load control GDI data
load('GDI_Calc_Fast_Data.mat');

%% If kinematic data already present, skip data loading
if KineData ~= 0
    L_PelvTilt = KineData(:,1);
    L_PelvOblq = KineData(:,2);
    L_PelvRot = KineData(:,3);
    L_HipFE = KineData(:,7);
    L_HipAbAd = KineData(:,8);
    L_HipRot = KineData(:,9);
    L_KneeFE = KineData(:,13);
    L_AnkleDP = KineData(:,19);
    L_FootProg = KineData(:,20);
    R_PelvTilt = KineData(:,4);
    R_PelvOblq = KineData(:,5);
    R_PelvRot = KineData(:,6);
    R_HipFE = KineData(:,10);
    R_HipAbAd = KineData(:,11);
    R_HipRot = KineData(:,12);
    R_KneeFE = KineData(:,16);
    R_AnkleDP = KineData(:,22);
    R_FootProg = KineData(:,23);
    
else  
    %% Kate's Code to get the Representative Trial of the Subject
    %Prompt user to choose what subject to extract kinematic data for.
    % [file2get, pathname] = uigetfile('*.c3d', 'Please select SUBJECT C3D file to be analyzed');
    % Batch process all input file(s)
    [Representative] = GDI_BatchProc_SubjTrials(files2get);
    
    for i = 1:2
        if i ==1
            % Representative Data
            % most representative trial (will be first value trial order)
            num_rep_trial_L = Representative.trial_order.left(1);
        else
            % Representative Data
            % most representative trial (will be first value trial order)
            num_rep_trial_R = Representative.trial_order.right(1);
        end
    end
    
    [~,N] = size(Representative.c3d_files(num_rep_trial_L).kinematics(1).data.left);
    while L_Rep > N
        L_Rep = L_Rep - 1;
    end
    
      [~,N] = size(Representative.c3d_files(num_rep_trial_L).kinematics(1).data.right);
    while R_Rep > N
        R_Rep = R_Rep - 1;
    end
    
    % Kinematics from representative trial's most representative gait cycle
    for i = 1:2
        if i ==1
            %LEFT
            L_PelvTilt = Representative.c3d_files(num_rep_trial_L).kinematics(1).data.left(:,L_Rep);
            L_PelvOblq = Representative.c3d_files(num_rep_trial_L).kinematics(2).data.left(:,L_Rep);
            L_PelvRot = Representative.c3d_files(num_rep_trial_L).kinematics(3).data.left(:,L_Rep);
            L_HipFE = Representative.c3d_files(num_rep_trial_L).kinematics(4).data.left(:,L_Rep);
            L_HipAbAd = Representative.c3d_files(num_rep_trial_L).kinematics(5).data.left(:,L_Rep);
            L_HipRot = Representative.c3d_files(num_rep_trial_L).kinematics(6).data.left(:,L_Rep);
            L_KneeFE = Representative.c3d_files(num_rep_trial_L).kinematics(7).data.left(:,L_Rep);
            L_AnkleDP = Representative.c3d_files(num_rep_trial_L).kinematics(8).data.left(:,L_Rep);
            L_FootProg = Representative.c3d_files(num_rep_trial_L).kinematics(9).data.left(:,L_Rep);
        else
            %RIGHT
            R_PelvTilt = Representative.c3d_files(num_rep_trial_R).kinematics(1).data.right(:,R_Rep);
            R_PelvOblq = Representative.c3d_files(num_rep_trial_R).kinematics(2).data.right(:,R_Rep);
            R_PelvRot = Representative.c3d_files(num_rep_trial_R).kinematics(3).data.right(:,R_Rep);
            R_HipFE = Representative.c3d_files(num_rep_trial_R).kinematics(4).data.right(:,R_Rep);
            R_HipAbAd = Representative.c3d_files(num_rep_trial_R).kinematics(5).data.right(:,R_Rep);
            R_HipRot = Representative.c3d_files(num_rep_trial_R).kinematics(6).data.right(:,R_Rep);
            R_KneeFE = Representative.c3d_files(num_rep_trial_R).kinematics(7).data.right(:,R_Rep);
            R_AnkleDP = Representative.c3d_files(num_rep_trial_R).kinematics(8).data.right(:,R_Rep);
            R_FootProg = Representative.c3d_files(num_rep_trial_R).kinematics(9).data.right(:,R_Rep);
        end
    end
end

%% If length of data is 101 points, decrease to 51
if length(L_PelvTilt) == 101
    x = 2:2:100;
    L_PelvTilt(x) = []; 
    L_PelvOblq(x) = [];
    L_PelvRot(x) = [];
    L_HipFE(x) = []; 
    L_HipAbAd(x) = [];
    L_HipRot(x) = [];
    L_KneeFE(x) = []; 
    L_AnkleDP(x) = [];
    L_FootProg(x) = [];
    R_PelvTilt(x) = []; 
    R_PelvOblq(x) = [];
    R_PelvRot(x) = [];
    R_HipFE(x) = []; 
    R_HipAbAd(x) = [];
    R_HipRot(x) = [];
    R_KneeFE(x) = []; 
    R_AnkleDP(x) = [];
    R_FootProg(x) = [];
end

%% Create subject's Left and Right gait vector (transpose) 459x1
subj_gait_vector_L = [L_PelvTilt' L_PelvOblq' L_PelvRot' L_HipFE' L_HipAbAd' L_HipRot' L_KneeFE' L_AnkleDP' L_FootProg']';
subj_gait_vector_R = [R_PelvTilt' R_PelvOblq' R_PelvRot' R_HipFE' R_HipAbAd' R_HipRot' R_KneeFE' R_AnkleDP' R_FootProg']';

%% Calculate Feature Components of Control Population (c)
%Preallocate space for control data structures
for i = 1:15
    Ctrl_data(i).singular_vector = zeros(length(all_ctrl_gait_vectors),1); 
end
for i = 1:num_ctrl_subjects
   Ctrl_data(i).gait_vector = zeros(length(all_ctrl_gait_vectors),1);
   Ctrl_data(i).feature_components  = zeros(15,1);
end

%Place gait and singular vectors into structures 
for i = 1:num_ctrl_subjects
    for j = 1:459
        gait_vector(j) = all_ctrl_gait_vectors(j,i);
    end
    Ctrl_data(i).gait_vector = gait_vector';
end

for i = 1:15
    Ctrl_data(i).singular_vector = gait_features(:,i);
end

%Calculate feature components and store in a structure
feature_components = zeros(1,15);
for i = 1:num_ctrl_subjects
    for j = 1:15
        feature_components(j) = dot(Ctrl_data(i).gait_vector,Ctrl_data(j).singular_vector);
    end
    Ctrl_data(i).feature_components = feature_components;
end

%Find average feature component for each feature direction
c_avg = zeros(1,15);

for j = 1:15
    numer_feature_components = 0;
    for i = 1:num_ctrl_subjects
        numer_feature_components = numer_feature_components + Ctrl_data(i).feature_components(j);
    end
    c_avg(j) = numer_feature_components/num_ctrl_subjects;
end

%% Calculate Feature Components of subject of interest
% Calculate feature components
c_subj_L = zeros(1,15);
c_subj_R = zeros(1,15);
% c_subj_avg = zeros(1,15);
for j = 1:2
    if j == 1
        for i = 1:15
            %LEFT
            c_subj_L(i) =  dot(subj_gait_vector_L,gait_features(:,i));
        end
    else
        for i = 1:15
            %RIGHT
            c_subj_R(i) =  dot(subj_gait_vector_R,gait_features(:,i));
        end
    end
end

%% Calculate GDI from Eucledian Distance
for i = 1:2
    if i == 1
        %LEFT
        %Find Eucledian distance between subject and control average feature
        %components
        euc_dist_L = norm(c_subj_L - c_avg);
        %Raw GDI for subject
        subj_GDI_raw_L = log(euc_dist_L);
    else
        %RIGHT
        %Find Eucledian distance between subject and control average feature
        %components
        euc_dist_R = norm(c_subj_R - c_avg);
        %Raw GDI for subject
        subj_GDI_raw_R = log(euc_dist_R);
    end
end

%Raw GDI for each control
for i = 1:num_ctrl_subjects
   Ctrl_data(i).raw_GDI = log(norm(Ctrl_data(i).feature_components - c_avg)); 
end

%Place control raw GDI into a vector
ctrl_GDI = zeros(1,num_ctrl_subjects);
for i = 1:num_ctrl_subjects
   ctrl_GDI(i) = Ctrl_data(i).raw_GDI; 
end

%Find z-score of the subject
for i = 1:2
    if i == 1
        %LEFT
        z_subj_L = (subj_GDI_raw_L - mean(ctrl_GDI))/std(ctrl_GDI);
    else
        %RIGHT
        z_subj_R = (subj_GDI_raw_R - mean(ctrl_GDI))/std(ctrl_GDI);
    end 
end

%Calculate GDI from z-score
for i = 1:2
    if i == 1
        %LEFT
        L_GDI = 100 - 10*z_subj_L;
    else
        %RIGHT
        R_GDI = 100 - 10*z_subj_R;
    end 
end

%Find average GDI from each side
avg_GDI = (L_GDI + R_GDI)/2;
GDI = [L_GDI R_GDI avg_GDI];

end



