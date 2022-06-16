function [search_Pt_label] = defPoints_Labels(MkrSet_type)
% -------------------------------------------------------------------------
% DEFINE POINT LABELS - List of Point Labels for various tasks
% -------------------------------------------------------------------------
%
% AUTHOR: Kate Worster
% DATE: April 6, 2011
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
%
% NOTES: All left point labels should be an odd number and all right point
%        labels should be an even number.  Other functions rely on this
%        convention since they use the corresponding gait cycle (left or
%        right) to parse the marker data appropriately.  If a point is not
%        side dependent (e.g. sacrum marker), then default to both left and
%        right.
%


% % DEBUG LOOP
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
% 
% %%%% END DEBUG LOOP %%%%%


%% MARKER SETS
if (strcmp(MkrSet_type, 'lowerbody') == 1)
    % If user is wanting to search lower body marker sets
    
    %%%%% MARKERS %%%%%
    search_Pt_label(1).name = 'SACR';
    search_Pt_label(2).name = 'SACR';
    search_Pt_label(3).name = 'LASI';
    search_Pt_label(4).name = 'RASI';

    search_Pt_label(5).name = 'LFEP';  % PiG, OLGA
        search_Pt_label(5).alt_name1 = 'LHJC';  % Gaia, Gaia Helical
        search_Pt_label(5).alt_name2 = 'LFEz';  % old Gaia models
    search_Pt_label(6).name = 'RFEP';  % PiG, OLGA
        search_Pt_label(6).alt_name1 = 'RHJC';  % Gaia, Gaia Helical
        search_Pt_label(6).alt_name2 = 'RFEz';  % old Gaia models
    search_Pt_label(7).name = 'LTHI';
    search_Pt_label(8).name = 'RTHI';
    search_Pt_label(9).name = 'LKNE';
    search_Pt_label(10).name = 'RKNE';
    search_Pt_label(11).name = 'LFEO';  % PiG, OLGA
        search_Pt_label(11).alt_name1 = 'LKJC';  % Gaia, Gaia Helical
    search_Pt_label(12).name = 'RFEO';  % PiG, OLGA
        search_Pt_label(12).alt_name1 = 'RKJC';  % Gaia, Gaia Helical
    search_Pt_label(13).name = 'LTIB';
    search_Pt_label(14).name = 'RTIB';
    search_Pt_label(15).name = 'LTIO';
        search_Pt_label(15).alt_name1 = 'LAJC';  % Gaia, Gaia Helical
        search_Pt_label(15).alt_name2 = 'LTIBO';  % old Gaia models
    search_Pt_label(16).name = 'RTIO';  % PiG, OLGA
        search_Pt_label(16).alt_name1 = 'RAJC';  % Gaia, Gaia Helical
        search_Pt_label(16).alt_name2 = 'RTIBO';  % old Gaia models
    search_Pt_label(17).name = 'LANK';
    search_Pt_label(18).name = 'RANK';
    search_Pt_label(19).name = 'LHEE';
    search_Pt_label(20).name = 'RHEE';
    search_Pt_label(21).name = 'LTOE';
    search_Pt_label(22).name = 'RTOE';
    
    
    %%%%% KINEMATICS %%%%%
    search_Pt_label(23).name = 'LPelvisAngles';
    search_Pt_label(24).name = 'RPelvisAngles';
    search_Pt_label(25).name = 'LHipAngles';
    search_Pt_label(26).name = 'RHipAngles';
    search_Pt_label(27).name = 'LKneeAngles';
    search_Pt_label(28).name = 'RKneeAngles';
    search_Pt_label(29).name = 'LAnkleAngles';  
    search_Pt_label(30).name = 'RAnkleAngles'; 
    search_Pt_label(31).name = 'LFootProgressAngles';
        search_Pt_label(31).alt_name1 = 'LFootProgAngle';
        search_Pt_label(31).alt_name2 = 'LFootProgressAngl';
    search_Pt_label(32).name = 'RFootProgressAngles';
        search_Pt_label(32).alt_name1 = 'RFootProgAngle';
        search_Pt_label(32).alt_name2 = 'RFootProgressAngl';
%     search_Pt_label(33).name = 'LShankAngles';  % Gaia, GaiaHelical
%     search_Pt_label(34).name = 'RShankAngles';  % Gaia, GaiaHelical

        
    %%%%% KINETICS %%%%%
    search_Pt_label(35).name = 'LHipForce';
    search_Pt_label(36).name = 'RHipForce';
    search_Pt_label(37).name = 'LKneeForce';
    search_Pt_label(38).name = 'RKneeForce';
    search_Pt_label(39).name = 'LAnkleForce';
    search_Pt_label(40).name = 'RAnkleForce';
    search_Pt_label(41).name = 'LHipMoment';  
    search_Pt_label(42).name = 'RHipMoment';  
    search_Pt_label(43).name = 'LKneeMoment';  
    search_Pt_label(44).name = 'RKneeMoment'; 
    search_Pt_label(45).name = 'LAnkleMoment';
    search_Pt_label(46).name = 'RAnkleMoment';
    search_Pt_label(47).name = 'LHipPower';
    search_Pt_label(48).name = 'RHipPower';
    search_Pt_label(49).name = 'LKneePower';
    search_Pt_label(50).name = 'RKneePower';
    search_Pt_label(51).name = 'LAnklePower';
    search_Pt_label(52).name = 'RAnklePower';


elseif (strcmp(MkrSet_type, 'upperbody') == 1)
    % If user is wanting to search upper body marker sets
    
    %%%%% MARKERS %%%%%
    search_Pt_label(1).name = 'LFHD';
    search_Pt_label(2).name = 'RFHD';
    search_Pt_label(3).name = 'LBHD';
    search_Pt_label(4).name = 'RBHD';
    search_Pt_label(5).name = 'C7';
    search_Pt_label(6).name = 'C7';
    search_Pt_label(7).name = 'T10';
    search_Pt_label(8).name = 'T10';
    search_Pt_label(9).name = 'CLAV';
    search_Pt_label(10).name = 'CLAV';
    search_Pt_label(11).name = 'STRN';
    search_Pt_label(12).name = 'STRN';
    search_Pt_label(13).name = 'LBAK';
    search_Pt_label(14).name = 'RBAK';

    search_Pt_label(15).name = 'LSHO';
    search_Pt_label(16).name = 'RSHO';    
    search_Pt_label(17).name = 'LUPA';
    search_Pt_label(18).name = 'RUPA';
    search_Pt_label(19).name = 'LUPB';  % Gaia
    search_Pt_label(20).name = 'RUPB';  % Gaia
    search_Pt_label(21).name = 'LUPC';  % Gaia
    search_Pt_label(22).name = 'RUPC';  % Gaia
    search_Pt_label(23).name = 'LELB';
    search_Pt_label(24).name = 'RELB';
    search_Pt_label(25).name = 'LFOR';  % Gaia
        search_Pt_label(25).alt_name1 = 'LFRA';  % AddRadius
    search_Pt_label(26).name = 'RFOR';  % Gaia
        search_Pt_label(26).alt_name1 = 'RFRA';  % AddRadius
    search_Pt_label(27).name = 'LRAD';  % Gaia
        search_Pt_label(27).alt_name1 = 'LWRA';  % AddRadius
    search_Pt_label(28).name = 'RRAD';  % Gaia
        search_Pt_label(28).alt_name1 = 'RWRA';  % AddRadius
    search_Pt_label(29).name = 'LULN';  % Gaia
        search_Pt_label(29).alt_name1 = 'LWRB';  % AddRadius
    search_Pt_label(30).name = 'RULN';  % Gaia
        search_Pt_label(30).alt_name1 = 'RWRB';  % AddRadius;
    search_Pt_label(31).name = 'LMC3';  % Gaia
        search_Pt_label(31).alt_name1 = 'LFIN';  % AddRadius
    search_Pt_label(32).name = 'RMC3';  % Gaia
        search_Pt_label(32).alt_name1 = 'RFIN';  % AddRadius
    search_Pt_label(33).name = 'LFIN1';  % Gaia, AddRadius
    search_Pt_label(34).name = 'RFIN1';  % Gaia, AddRadius
    search_Pt_label(33).name = 'LTHM';  % Gaia, AddRadius
    search_Pt_label(36).name = 'RTHM';  % Gaia, AddRadius
    
    % Object Markers
    search_Pt_label(37).name = 'OBJ1';  % Gaia
    search_Pt_label(38).name = 'OBJ2';  % Gaia
    
    
    %%%%% KINEMATICS %%%%%
    search_Pt_label(39).name = 'LHeadAngles';
    search_Pt_label(40).name = 'RHeadAngles';
    search_Pt_label(41).name = 'LNeckAngles';
    search_Pt_label(42).name = 'RNeckAngles';
    search_Pt_label(43).name = 'LThoraxAngles';
    search_Pt_label(44).name = 'RThoraxAngles';
    search_Pt_label(45).name = 'LSpineAngles';
    search_Pt_label(46).name = 'RSpineAngles';
    search_Pt_label(47).name = 'LClavAngles';
    search_Pt_label(48).name = 'RClavAngles';
    search_Pt_label(49).name = 'LShoulderAngles';
    search_Pt_label(50).name = 'RShoulderAngles';
    search_Pt_label(51).name = 'LElbowAngles';
    search_Pt_label(52).name = 'RElbowAngles';
    search_Pt_label(53).name = 'LForearmAngles';
    search_Pt_label(54).name = 'RForearmAngles';
    search_Pt_label(55).name = 'LWristAngles';
    search_Pt_label(56).name = 'RWristAngles';


elseif (strcmp(MkrSet_type, 'fullbody') == 1)
    %%%%% MARKERS %%%%%
    
    % Upper body markers
    search_Pt_label(1).name = 'LFHD';
    search_Pt_label(2).name = 'RFHD';
    search_Pt_label(3).name = 'LBHD';
    search_Pt_label(4).name = 'RBHD';
    search_Pt_label(5).name = 'C7';
    search_Pt_label(6).name = 'C7';
    search_Pt_label(7).name = 'T10';
    search_Pt_label(8).name = 'T10';
    search_Pt_label(9).name = 'CLAV';
    search_Pt_label(10).name = 'CLAV';
    search_Pt_label(11).name = 'STRN';
    search_Pt_label(12).name = 'STRN';
    search_Pt_label(13).name = 'LBAK';
    search_Pt_label(14).name = 'RBAK';

    search_Pt_label(15).name = 'LSHO';
    search_Pt_label(16).name = 'RSHO';    
    search_Pt_label(17).name = 'LUPA';
    search_Pt_label(18).name = 'RUPA';
    search_Pt_label(19).name = 'LUPB';  % Gaia
    search_Pt_label(20).name = 'RUPB';  % Gaia
    search_Pt_label(21).name = 'LUPC';  % Gaia
    search_Pt_label(22).name = 'RUPC';  % Gaia
    search_Pt_label(23).name = 'LELB';
    search_Pt_label(24).name = 'RELB';
    search_Pt_label(25).name = 'LFOR';  % Gaia
        search_Pt_label(25).alt_name1 = 'LFRA';  % AddRadius
    search_Pt_label(26).name = 'RFOR';  % Gaia
        search_Pt_label(26).alt_name1 = 'RFRA';  % AddRadius
    search_Pt_label(27).name = 'LRAD';  % Gaia
        search_Pt_label(27).alt_name1 = 'LWRA';  % AddRadius
    search_Pt_label(28).name = 'RRAD';  % Gaia
        search_Pt_label(28).alt_name1 = 'RWRA';  % AddRadius
    search_Pt_label(29).name = 'LULN';  % Gaia
        search_Pt_label(29).alt_name1 = 'LWRB';  % AddRadius
    search_Pt_label(30).name = 'RULN';  % Gaia
        search_Pt_label(30).alt_name1 = 'RWRB';  % AddRadius;
    search_Pt_label(31).name = 'LMC3';  % Gaia
        search_Pt_label(31).alt_name1 = 'LFIN';  % AddRadius
    search_Pt_label(32).name = 'RMC3';  % Gaia
        search_Pt_label(32).alt_name1 = 'RFIN';  % AddRadius
    search_Pt_label(33).name = 'LFIN1';  % Gaia, AddRadius
    search_Pt_label(34).name = 'RFIN1';  % Gaia, AddRadius
    search_Pt_label(33).name = 'LTHM';  % Gaia, AddRadius
    search_Pt_label(36).name = 'RTHM';  % Gaia, AddRadius
    
    % Lower body markers
    search_Pt_label(37).name = 'SACR';
    search_Pt_label(38).name = 'SACR';
    search_Pt_label(39).name = 'LASI';
    search_Pt_label(40).name = 'RASI';

    search_Pt_label(41).name = 'LFEP';  % PiG, OLGA
        search_Pt_label(41).alt_name1 = 'LHJC';  % Gaia, Gaia Helical
        search_Pt_label(41).alt_name2 = 'LFEz';  % old Gaia models
    search_Pt_label(42).name = 'RFEP';  % PiG, OLGA
        search_Pt_label(42).alt_name1 = 'RHJC';  % Gaia, Gaia Helical
        search_Pt_label(42).alt_name2 = 'RFEz';  % old Gaia models
    search_Pt_label(43).name = 'LTHI';
    search_Pt_label(44).name = 'RTHI';
    search_Pt_label(45).name = 'LKNE';
    search_Pt_label(46).name = 'RKNE';
    search_Pt_label(47).name = 'LFEO';  % PiG, OLGA
        search_Pt_label(47).alt_name1 = 'LKJC';  % Gaia, Gaia Helical
    search_Pt_label(48).name = 'RFEO';  % PiG, OLGA
        search_Pt_label(48).alt_name1 = 'RKJC';  % Gaia, Gaia Helical
    search_Pt_label(49).name = 'LTIB';
    search_Pt_label(50).name = 'RTIB';
    search_Pt_label(51).name = 'LTIO';
        search_Pt_label(51).alt_name1 = 'LAJC';  % Gaia, Gaia Helical
        search_Pt_label(51).alt_name2 = 'LTIBO';  % old Gaia models
    search_Pt_label(52).name = 'RTIO';  % PiG, OLGA
        search_Pt_label(52).alt_name1 = 'RAJC';  % Gaia, Gaia Helical
        search_Pt_label(52).alt_name1 = 'RTIBO';  % old Gaia models
    search_Pt_label(53).name = 'LANK';
    search_Pt_label(54).name = 'RANK';
    search_Pt_label(55).name = 'LHEE';
    search_Pt_label(56).name = 'RHEE';
    search_Pt_label(57).name = 'LTOE';
    search_Pt_label(58).name = 'RTOE';
    
    
    %%%%% KINEMATICS %%%%%
    search_Pt_label(59).name = 'LHeadAngles';
    search_Pt_label(60).name = 'RHeadAngles';
    search_Pt_label(61).name = 'LNeckAngles';
    search_Pt_label(62).name = 'RNeckAngles';
    search_Pt_label(63).name = 'LThoraxAngles';
    search_Pt_label(64).name = 'RThoraxAngles';
    search_Pt_label(65).name = 'LSpineAngles';
    search_Pt_label(66).name = 'RSpineAngles';
    search_Pt_label(67).name = 'LClavAngles';
    search_Pt_label(68).name = 'RClavAngles';
    search_Pt_label(69).name = 'LShoulderAngles';
    search_Pt_label(70).name = 'RShoulderAngles';
    search_Pt_label(71).name = 'LElbowAngles';
    search_Pt_label(72).name = 'RElbowAngles';
    search_Pt_label(73).name = 'LForearmAngles';
    search_Pt_label(74).name = 'RForearmAngles';
    search_Pt_label(75).name = 'LWristAngles';
    search_Pt_label(76).name = 'RWristAngles';
    
    search_Pt_label(77).name = 'LPelvisAngles';
    search_Pt_label(78).name = 'RPelvisAngles';
    search_Pt_label(79).name = 'LHipAngles';
    search_Pt_label(80).name = 'RHipAngles';
    search_Pt_label(81).name = 'LKneeAngles';
    search_Pt_label(82).name = 'RKneeAngles';
    search_Pt_label(83).name = 'LAnkleAngles';  
    search_Pt_label(84).name = 'RAnkleAngles'; 
    search_Pt_label(85).name = 'LFootProgressAngles';
        search_Pt_label(85).alt_name1 = 'LFootProgAngle';
        search_Pt_label(85).alt_name2 = 'LFootProgressAngl';
    search_Pt_label(86).name = 'RFootProgressAngles';
        search_Pt_label(86).alt_name1 = 'RFootProgAngle';
        search_Pt_label(86).alt_name2 = 'RFootProgressAngl';
    search_Pt_label(87).name = 'LShankAngles';  % Gaia, GaiaHelical
    search_Pt_label(88).name = 'RShankAngles';  % Gaia, GaiaHelical

        
    %%%%% KINETICS %%%%%
    search_Pt_label(89).name = 'LHipForce';
    search_Pt_label(90).name = 'RHipForce';
    search_Pt_label(91).name = 'LKneeForce';
    search_Pt_label(92).name = 'RKneeForce';
    search_Pt_label(93).name = 'LAnkleForce';
    search_Pt_label(94).name = 'RAnkleForce';
    search_Pt_label(95).name = 'LHipMoment';  
    search_Pt_label(96).name = 'RHipMoment';  
    search_Pt_label(97).name = 'LKneeMoment';  
    search_Pt_label(98).name = 'RKneeMoment'; 
    search_Pt_label(99).name = 'LAnkleMoment';
    search_Pt_label(100).name = 'RAnkleMoment';
    search_Pt_label(101).name = 'LHipPower';
    search_Pt_label(102).name = 'RHipPower';
    search_Pt_label(103).name = 'LKneePower';
    search_Pt_label(104).name = 'RKneePower';
    search_Pt_label(105).name = 'LAnklePower';
    search_Pt_label(106).name = 'RAnklePower';    
    
    
elseif (strcmp(MkrSet_type, 'Dionysos') == 1)
    % For Dionysos Program
    
    % Upper body markers
    search_Pt_label(1).name = 'LHEO';
    search_Pt_label(2).name = 'RHEO';
    search_Pt_label(3).name = 'LHEP';
        search_Pt_label(3).alt_name1 = 'LHEz';
    search_Pt_label(4).name = 'RHEP';
        search_Pt_label(4).alt_name1 = 'RHEz';
    search_Pt_label(5).name = 'LHEA';
        search_Pt_label(5).alt_name1 = 'LHEx';
    search_Pt_label(6).name = 'RHEA';
        search_Pt_label(6).alt_name1 = 'RHEx';
    search_Pt_label(7).name = 'LSJC';
    search_Pt_label(8).name = 'RSJC';
    search_Pt_label(9).name = 'LHUP';
        search_Pt_label(9).alt_name1 = 'LHUz';
    search_Pt_label(10).name = 'RHUP';
        search_Pt_label(10).alt_name1 = 'RHUz';
    search_Pt_label(11).name = 'LEJC';
    search_Pt_label(12).name = 'REJC';
    search_Pt_label(13).name = 'LWJC';
    search_Pt_label(14).name = 'RWJC';
    search_Pt_label(15).name = 'LFIN';
        search_Pt_label(15).alt_name1 = 'LMC3';
    search_Pt_label(16).name = 'RFIN';
        search_Pt_label(16).alt_name1 = 'RMC3';
    
    % Lower body markers
    search_Pt_label(17).name = 'PELO';
    search_Pt_label(18).name = 'PELO';
    search_Pt_label(19).name = 'LFEP';
    search_Pt_label(20).name = 'RFEP';
    search_Pt_label(21).name = 'LFEO';
    search_Pt_label(22).name = 'RFEO';
    search_Pt_label(23).name = 'LTIO';
    search_Pt_label(24).name = 'RTIO';
    search_Pt_label(25).name = 'LHEE';
    search_Pt_label(26).name = 'RHEE';
    search_Pt_label(27).name = 'LTOE';
    search_Pt_label(28).name = 'RTOE';    
        
        
elseif (strcmp(MkrSet_type, 'Helios') == 1)
    % For Helios Model
    search_Pt_label(1).name = 'SACR';
    search_Pt_label(2).name = 'SACR';
    search_Pt_label(3).name = 'LASI';
    search_Pt_label(4).name = 'RASI';

    search_Pt_label(5).name = 'LFEP';  % PiG, OLGA
        search_Pt_label(5).alt_name1 = 'LHJC';  % Gaia, Gaia Helical
        search_Pt_label(5).alt_name2 = 'LFEz';  % old Gaia models
    search_Pt_label(6).name = 'RFEP';  % PiG, OLGA
        search_Pt_label(6).alt_name1 = 'RHJC';  % Gaia, Gaia Helical
        search_Pt_label(6).alt_name2 = 'RFEz';  % old Gaia models
    search_Pt_label(7).name = 'LGT';
    search_Pt_label(8).name = 'RGT';
    search_Pt_label(9).name = 'LTHI';
    search_Pt_label(10).name = 'RTHI';
    search_Pt_label(11).name = 'LITB';
    search_Pt_label(12).name = 'RITB';
    search_Pt_label(13).name = 'LKNE';
    search_Pt_label(14).name = 'RKNE';
    search_Pt_label(15).name = 'LFEO';  % PiG, OLGA
        search_Pt_label(15).alt_name1 = 'LKJC';  % Gaia, Gaia Helical
    search_Pt_label(16).name = 'RFEO';  % PiG, OLGA
        search_Pt_label(16).alt_name1 = 'RKJC';  % Gaia, Gaia Helical
    search_Pt_label(17).name = 'LTT';
    search_Pt_label(18).name = 'RTT';
    search_Pt_label(19).name = 'LTC';
    search_Pt_label(20).name = 'RTC';
    search_Pt_label(21).name = 'LTIB';
    search_Pt_label(22).name = 'RTIB';
    search_Pt_label(23).name = 'LTIO';
        search_Pt_label(23).alt_name1 = 'LAJC';  % Gaia, Gaia Helical
        search_Pt_label(23).alt_name1 = 'LTIBO';  % old Gaia models
    search_Pt_label(24).name = 'RTIO';  % PiG, OLGA
        search_Pt_label(24).alt_name1 = 'RAJC';  % Gaia, Gaia Helical
        search_Pt_label(24).alt_name1 = 'RTIBO';  % old Gaia models
    search_Pt_label(25).name = 'LANK';
    search_Pt_label(26).name = 'RANK';
    search_Pt_label(27).name = 'LHEE';
    search_Pt_label(28).name = 'RHEE';
    search_Pt_label(29).name = 'LTOE';   
    search_Pt_label(30).name = 'RTOE';


        
end
