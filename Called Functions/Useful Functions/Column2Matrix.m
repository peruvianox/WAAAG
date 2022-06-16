function [DataOut] = Column2Matrix(DataIn)
% -------------------------------------------------------------------------
% Reorganize C3D kinematic data to 101x24 Matrix
% -------------------------------------------------------------------------
% Author: Ricky Pimentel
% Date: 12/1/2016
% For use at the Center for Gait and Movement Analysis Laboratory
% at the Children's Hospital in Denver, CO.
%
% Description: This function will translate kinematic data from Schwartz/Baker's columns
% of kinematic data for left and right sides into a matrix that is used by
% CGMA to analyze kinematics OR vice versa.
%
% Input:     DataIn                 Either a 918x1 or 51x24 Matrix
%
% Output:   DataOut               Whichever organization DataIn is not
% 

% The column is composed of 9 sets of 51-row kinematic data for each limb,
% or a 918x1matrix conatining kinematic data for each joint for GDI
% calcualtions.

% The matrix will be arranged in a 51x24 with 6 empty columns due to the 3
% missing curves from each side. 51*18 = 918 data points

% Variable definitions for reference
% P = pelvis, H = hip, K = knee, A = ankle
% X = sagittal plane, Y = frontal plane, Z = transverse plane
% L = Left, R = Right

if size(DataIn) == [918,1]
    
    LPX = DataIn(1:51);
    LPY = DataIn(52:102);
    LPZ = DataIn(103:153);
    LHX = DataIn(154:204);
    LHY = DataIn(205:255);
    LHZ = DataIn(256:306);
    LKX = DataIn(307:357);
    LKY = zeros(51,1);
    LKZ = zeros(51,1);
    LAX = DataIn(358:408);
    LAY = DataIn(409:459);
    LAZ = zeros(51,1);
    
    RPX = DataIn(460:510);
    RPY = DataIn(511:561);
    RPZ = DataIn(562:612);
    RHX = DataIn(613:663);
    RHY = DataIn(664:714);
    RHZ = DataIn(715:765);
    RKX = DataIn(766:816);
    RKY = zeros(51,1);
    RKZ = zeros(51,1);
    RAX = DataIn(817:867);
    RAY = DataIn(868:918);
    RAZ = zeros(51,1);
    
    DataOut = [LPX,LPY,LPZ,RPX,RPY,RPZ,LHX,LHY,LHZ,RHX,RHY,RHZ,...
        LKX,LKY,LKZ,RKX,RKY,RKZ,LAX,LAY,LAZ,RAX,RAY,RAZ];
    
elseif size(DataIn) == [51,24]
    
    
    LPX = DataIn(:,1);
    LPY = DataIn(:,2);
    LPZ = DataIn(:,3);
    
    RPX = DataIn(:,4);
    RPY = DataIn(:,5);
    RPZ = DataIn(:,6);
    
    LHX = DataIn(:,7);
    LHY = DataIn(:,8);
    LHZ = DataIn(:,9);
    
    RHX = DataIn(:,10);
    RHY = DataIn(:,11);
    RHZ = DataIn(:,12);
    
    LKX = DataIn(:,13);
    
    RKX = DataIn(:,16);
    
    LAX = DataIn(:,19);
    LAY = DataIn(:,20);
    
    RAX = DataIn(:,22);
    RAY = DataIn(:,23);
    
    
    DataOut = [LPX;LPY;LPZ;LHX;LHY;LHZ;LKX;LAX;LAY;...
        RPX;RPY;RPZ;RHX;RHY;RHZ;RKX;RAX;RAY];
    
elseif size(DataIn) == [101,24]
    x = 2:2:100;
    DataIn(x,:) = [];
    
    LPX = DataIn(:,1);
    LPY = DataIn(:,2);
    LPZ = DataIn(:,3);
    
    RPX = DataIn(:,4);
    RPY = DataIn(:,5);
    RPZ = DataIn(:,6);
    
    LHX = DataIn(:,7);
    LHY = DataIn(:,8);
    LHZ = DataIn(:,9);
    
    RHX = DataIn(:,10);
    RHY = DataIn(:,11);
    RHZ = DataIn(:,12);
    
    LKX = DataIn(:,13);
    
    RKX = DataIn(:,16);
    
    LAX = DataIn(:,19);
    LAY = DataIn(:,20);
    
    RAX = DataIn(:,22);
    RAY = DataIn(:,23);
    
    
    DataOut = [LPX;LPY;LPZ;LHX;LHY;LHZ;LKX;LAX;LAY;...
        RPX;RPY;RPZ;RHX;RHY;RHZ;RKX;RAX;RAY];
end

end

