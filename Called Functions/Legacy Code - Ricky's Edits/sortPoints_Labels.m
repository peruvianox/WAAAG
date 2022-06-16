function [search_Pt_label] = defPoints_Labels(MkrSet_type)

% -------------------------------------------------------------------------
% DEFINE POINT LABELS - List of Point Labels for various tasks
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: April 6, 2011
%
% DESCRIPTION:  Pairs the xyz coordinates of a point with point
%               labels.  
%
% INPUTS: 
%           MkrSet_type             type of marker set being used
%           
%
% OUTPUTS: 
%           sortedPoints            cells of each point's xyz coordinates
%                                   with corresponding point label

%% DEBUG LOOP
% close all
% clear all
% clc
% 
% % Prompts User for Input
% [FullFileName, pathname] = uigetfile('*.c3d', 'Please select C3D file to be analyzed');
% 
% [C3Dfile] = OpenC3D(FullFileName);
% 
% % Open and Read the user specified C3D file
% % and returns: HeaderInfo, Point_xyz, PointLabels
% [HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);
% 
% % alternative string
% alt_str = '';

%%%%% END DEBUG LOOP %%%%%


%% MARKER SETS
if (strcmp(MkrSet_type, 'lowerbody') == 1)
    % If user is wanting to search lower body marker sets
    
    %%%%% MARKERS %%%%%
    search_Pt_label(1).name = 'SACR';
    search_Pt_label(2).name = 'LASI';
    search_Pt_label(3).name = 'RASI';

    search_Pt_label(4).name = 'LFEP';  % PiG, OLGA
        search_Pt_label(4).alt_name1 = 'LHJC';  % Gaia, Gaia Helical
        search_Pt_label(4).alt_name2 = 'LFEz';  % old Gaia models
    search_Pt_label(5).name = 'LTHI';
    search_Pt_label(6).name = 'LKNE';
    search_Pt_label(7).name = 'LFEO';  % PiG, OLGA
        search_Pt_label(7).alt_name1 = 'LKJC';  % Gaia, Gaia Helical
    search_Pt_label(8).name = 'LTIB';
    search_Pt_label(9).name = 'LTIO';
        search_Pt_label(9).alt_name1 = 'LAJC';  % Gaia, Gaia Helical
        search_Pt_label(9).alt_name2 = 'LTIBO';  % old Gaia models
    search_Pt_label(10).name = 'LANK';
    search_Pt_label(11).name = 'LHEE';
    search_Pt_label(12).name = 'LTOE';

    search_Pt_label(13).name = 'RFEP';  % PiG, OLGA
        search_Pt_label(13).alt_name1 = 'RHJC';  % Gaia, Gaia Helical
        search_Pt_label(13).alt_name2 = 'RFEz';  % old Gaia models
    search_Pt_label(14).name = 'RTHI';
    search_Pt_label(15).name = 'RKNE';
    search_Pt_label(16).name = 'RFEO';  % PiG, OLGA
        search_Pt_label(16).alt_name1 = 'RKJC';  % Gaia, Gaia Helical
    search_Pt_label(17).name = 'RTIB';
    search_Pt_label(18).name = 'RTIO';  % PiG, OLGA
        search_Pt_label(18).alt_name1 = 'RAJC';  % Gaia, Gaia Helical
        search_Pt_label(18).alt_name1 = 'RTIBO';  % old Gaia models
    search_Pt_label(19).name = 'RANK';
    search_Pt_label(20).name = 'RHEE';
    search_Pt_label(21).name = 'RTOE';
    
    
    %%%%% KINEMATICS %%%%%
    search_Pt_label(22).name = 'LPelvisAngles';
    search_Pt_label(23).name = 'RPelvisAngles';
    search_Pt_label(24).name = 'LHipAngles';
    search_Pt_label(25).name = 'RHipAngles';
    search_Pt_label(26).name = 'LKneeAngles';
    search_Pt_label(27).name = 'RKneeAngles';
    search_Pt_label(28).name = 'LShankAngles';  % Gaia, GaiaHelical
    search_Pt_label(29).name = 'RShankAngles';  % Gaia, GaiaHelical
    search_Pt_label(30).name = 'LAnkleAngles';  
    search_Pt_label(31).name = 'RAnkleAngles'; 
    search_Pt_label(32).name = 'LFootProgressAngles';
        search_Pt_label(32).alt_name1 = 'LFootProgAngle';
        search_Pt_label(32).alt_name2 = 'LFootProgressAngl';
    search_Pt_label(33).name = 'RFootProgressAngles';
        search_Pt_label(33).alt_name1 = 'RFootProgAngle';
        search_Pt_label(33).alt_name2 = 'RFootProgressAngl';

        
    %%%%% KINETICS %%%%%
    search_Pt_label(34).name = 'LHipForce';
    search_Pt_label(35).name = 'RHipForce';
    search_Pt_label(36).name = 'LKneeForce';
    search_Pt_label(37).name = 'RKneeForce';
    search_Pt_label(38).name = 'LAnkleForce';
    search_Pt_label(39).name = 'RAnkleForce';
    search_Pt_label(40).name = 'LHipMoment';  
    search_Pt_label(41).name = 'RHipMoment';  
    search_Pt_label(42).name = 'LKneeMoment';  
    search_Pt_label(43).name = 'RKneeMoment'; 
    search_Pt_label(44).name = 'LAnkleMoment';
    search_Pt_label(45).name = 'RAnkleMoment';
    search_Pt_label(46).name = 'LHipPower';
    search_Pt_label(47).name = 'RHipPower';
    search_Pt_label(48).name = 'LKneePower';
    search_Pt_label(49).name = 'RKneePower';
    search_Pt_label(50).name = 'LAnklePower';
    search_Pt_label(51).name = 'RAnklePower';


elseif (strcmop(MkrSet_type, 'upperbody') == 1)
    % If user is wanting to search upper body marker sets
    
    %%%%% MARKERS %%%%%
    search_Pt_label(1).name = 'LFHD';
    search_Pt_label(2).name = 'RFHD';
    search_Pt_label(3).name = 'LBHD';
    search_Pt_label(4).name = 'RBHD';
    search_Pt_label(5).name = 'C7';
    search_Pt_label(6).name = 'T10';
    search_Pt_label(7).name = 'CLAV';
    search_Pt_label(8).name = 'STRN';
    search_Pt_label(9).name = 'RBAK';

    search_Pt_label(10).name = 'LSHO';
    search_Pt_label(11).name = 'LUPA';
    search_Pt_label(12).name = 'LUPB';  % Gaia
    search_Pt_label(13).name = 'LUPC';  % Gaia
    search_Pt_label(14).name = 'LELB';
    search_Pt_label(15).name = 'LFOR';  % Gaia
        search_Pt_label(15).alt_name1 = 'LFRA';  % AddRadius
    search_Pt_label(16).name = 'LRAD';  % Gaia
        search_Pt_label(16).alt_name1 = 'LWRA';  % AddRadius
    search_Pt_label(17).name = 'LULN';  % Gaia
        search_Pt_label(17).alt_name1 = 'LWRB';  % AddRadius
    search_Pt_label(18).name = 'LMC3';  % Gaia
        search_Pt_label(18).alt_name1 = 'LFIN';  % AddRadius
    search_Pt_label(19).name = 'LFIN1';  % Gaia, AddRadius
    search_Pt_label(19).name = 'LTHM';  % Gaia, AddRadius

    search_Pt_label(20).name = 'RSHO';
    search_Pt_label(21).name = 'RUPA';
    search_Pt_label(22).name = 'RUPB';  % Gaia
    search_Pt_label(23).name = 'RUPC';  % Gaia
    search_Pt_label(24).name = 'RELB';
    search_Pt_label(25).name = 'RFOR';  % Gaia
        search_Pt_label(25).alt_name1 = 'RFRA';  % AddRadius
    search_Pt_label(26).name = 'RRAD';  % Gaia
        search_Pt_label(26).alt_name1 = 'RWRA';  % AddRadius
    search_Pt_label(27).name = 'RULN';  % Gaia
        search_Pt_label(27).alt_name1 = 'RWRB';  % AddRadius;
    search_Pt_label(28).name = 'RMC3';  % Gaia
        search_Pt_label(28).alt_name1 = 'RFIN';  % AddRadius
    search_Pt_label(29).name = 'RFIN1';  % Gaia, AddRadius
    search_Pt_label(29).name = 'RTHM';  % Gaia, AddRadius
    
    search_Pt_label(30).name = 'GFS1';  % Gaia
    search_Pt_label(31).name = 'GFS2';  % Gaia
    
    
    %%%%% KINEMATICS %%%%%
    search_Pt_label(32).name = 'LHeadAngles';
    search_Pt_label(33).name = 'RHeadAngles';
    search_Pt_label(34).name = 'LNeckAngles';
    search_Pt_label(35).name = 'RNeckAngles';
    search_Pt_label(36).name = 'LThoraxAngles';
    search_Pt_label(37).name = 'RThoraxAngles';
    search_Pt_label(38).name = 'LSpineAngles';
    search_Pt_label(39).name = 'RSpineAngles';
    search_Pt_label(40).name = 'LClavAngles';
    search_Pt_label(41).name = 'RClavAngles';
    search_Pt_label(42).name = 'LShoulderAngles';
    search_Pt_label(43).name = 'RShoulderAngles';
    search_Pt_label(44).name = 'LElbowAngles';
    search_Pt_label(45).name = 'RElbowAngles';
    search_Pt_label(46).name = 'LForearmAngles';
    search_Pt_label(47).name = 'RForearmAngles';
    search_Pt_label(48).name = 'LWristAngles';
    search_Pt_label(49).name = 'RWristAngles';


elseif (strcmp(MkrSet_type, 'fullbody') == 1)
    %%%%% MARKERS %%%%%
    
    
    %%%%% KNEMATICS %%%%%
    
    
    %%%%% KINETICS %%%%%
    
    
elseif (strcmp(MkrSet_type, 'Helios') == 1)
    % For Helios Model
    search_Pt_label(1).name = 'SACR';
    search_Pt_label(2).name = 'LASI';
    search_Pt_label(3).name = 'RASI';

    search_Pt_label(4).name = 'LFEP';  % PiG, OLGA
        search_Pt_label(4).alt_name1 = 'LHJC';  % Gaia, Gaia Helical
        search_Pt_label(4).alt_name2 = 'LFEz';  % old Gaia models
    search_Pt_label(5).name = 'LGT';
    search_Pt_label(6).name = 'LTHI';
    search_Pt_label(7).name = 'LITB';
    search_Pt_label(8).name = 'LKNE';
    search_Pt_label(9).name = 'LFEO';  % PiG, OLGA
        search_Pt_label(9).alt_name1 = 'LKJC';  % Gaia, Gaia Helical
    search_Pt_label(10).name = 'LTT';
    search_Pt_label(11).name = 'LTC';
    search_Pt_label(12).name = 'LTIB';
    search_Pt_label(13).name = 'LTIO';
        search_Pt_label(13).alt_name1 = 'LAJC';  % Gaia, Gaia Helical
        search_Pt_label(13).alt_name1 = 'LTIBO';  % old Gaia models
    search_Pt_label(14).name = 'LANK';
    search_Pt_label(15).name = 'LHEE';
    search_Pt_label(16).name = 'LTOE';

    search_Pt_label(17).name = 'RFEP';  % PiG, OLGA
        search_Pt_label(17).alt_name1 = 'RHJC';  % Gaia, Gaia Helical
        search_Pt_label(17).alt_name2 = 'RFEz';  % old Gaia models
    search_Pt_label(18).name = 'RGT';
    search_Pt_label(19).name = 'RTHI';
    search_Pt_label(20).name = 'RITB';
    search_Pt_label(21).name = 'RKNE';
    search_Pt_label(22).name = 'RFEO';  % PiG, OLGA
        search_Pt_label(22).alt_name1 = 'RKJC';  % Gaia, Gaia Helical
    search_Pt_label(23).name = 'RTT';
    search_Pt_label(24).name = 'RTC';
    search_Pt_label(25).name = 'RTIB';
    search_Pt_label(26).name = 'RTIO';  % PiG, OLGA
        search_Pt_label(26).alt_name1 = 'RAJC';  % Gaia, Gaia Helical
        search_Pt_label(26).alt_name1 = 'RTIBO';  % old Gaia models
    search_Pt_label(27).name = 'RANK';
    search_Pt_label(28).name = 'RHEE';
    search_Pt_label(29).name = 'RTOE';


        
end
