function [Kinematics, Kinetics] = sortPointLabels(Point_xyz, PointLabels)

% -------------------------------------------------------------------------
% Sorts PointLabels Group
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: April 11, 2010

% Updated 07Apr2016
% Kayla Burnim
% Compatible with Score/Sara _M structure
%
% DESCRIPTION:  Pairs the lower body, thorax, & spine kinematic angles with the appropriate labels. For
%               use with the Anterior Knee Pain Program. For PiG model.
%               Lower body kinetics as well.
%
% INPUTS: 
%           Point_xyz               xyz coordinates of points from C3D file
%           PointLabels             labels of all points in C3D file
%
% OUTPUTS: 
%           Kinematics              cells of each Kinematic angle with corresponding label
%           Kinetics                cells of each Kinetic group with corresponding label
%


%% DEBUG Loop
% close all
% clear all
% clc
% 
% [FullFileName, pathname] = uigetfile('*.*', 'Please select a C3D file to be analyzed');
% 
% % Opens and Extract C3D file's data
% [C3Dfile] = OpenC3D(FullFileName);
% 
% % Open and Read the user specified C3D file
% % and returns: HeaderInfo, Point_xyz, PointLabels
% [HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile);


%% KINMATICS 
% Find and organize Point_xyz and PointLabels 

% LEFT THORAX ANGLES
searchPtLabels_LThorax = strncmp(PointLabels, 'LThoraxAngles',10);
if sum(searchPtLabels_LThorax == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LThorax(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LThorax = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LThorax(:,:,:) = Point_xyz(:,ncolumn_LThorax,:);
else

    LThorax = 0;
end

% RIGHT THORAX ANGLES
searchPtLabels_RThorax = strncmp(PointLabels, 'RThoraxAngles',10);
if sum(searchPtLabels_RThorax == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RThorax(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RThorax = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RThorax(:,:,:) = Point_xyz(:,ncolumn_RThorax,:);
else

    RThorax = 0;
end


% LEFT SPINE ANGLES
searchPtLabels_LSpine = strncmp(PointLabels, 'LSpineAngles',10);
if sum(searchPtLabels_LSpine == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LSpine(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LSpine = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LSpine(:,:,:) = Point_xyz(:,ncolumn_LSpine,:);
else

    LSpine = 0;
end

% RIGHT SPINE ANGLES
searchPtLabels_RSpine = strncmp(PointLabels, 'RSpineAngles',10);
if sum(searchPtLabels_RSpine == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RSpine(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RSpine = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RSpine(:,:,:) = Point_xyz(:,ncolumn_RSpine,:);
else

    RSpine = 0;
end



% LEFT PELVIS ANGLES
searchPtLabels_LPelv = strncmp(PointLabels, 'LPelvisAngles',10);
if sum(searchPtLabels_LPelv == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LPelv(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LPelv = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LPelv(:,:,:) = Point_xyz(:,ncolumn_LPelv,:);
else
    LPelv = 0;
end


% RIGHT PELVIS ANGLES
searchPtLabels_RPelv = strncmp(PointLabels, 'RPelvisAngles',10);
if sum(searchPtLabels_RPelv == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RPelv(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RPelv = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RPelv(:,:,:) = Point_xyz(:,ncolumn_RPelv,:);
else
    RPelv = 0;
end


% LEFT HIP ANGLES
searchPtLabels_LHip = strncmp(PointLabels, 'LHipAngles',10);
if sum(searchPtLabels_LHip == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LHip(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LHip = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LHip(:,:,:) = Point_xyz(:,ncolumn_LHip,:);
else
%     disp('Point Label "LHipAngles" does not exist');
    LHip = 0;
end

% RIGHT HIP ANGLES
searchPtLabels_RHip = strncmp(PointLabels, 'RHipAngles',10);
if sum(searchPtLabels_RHip == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RHip(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RHip = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RHip(:,:,:) = Point_xyz(:,ncolumn_RHip,:);
else
%     disp('Point Label "RHipAngles" does not exist');
    RHip = 0;
end

% LEFT KNEE ANGLES
searchPtLabels_LKnee = strncmp(PointLabels, 'LKneeAngles',10);
if sum(searchPtLabels_LKnee == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LKnee(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LKnee = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LKnee(:,:,:) = Point_xyz(:,ncolumn_LKnee,:);
else
%     disp('Point Label "LKneeAngles" does not exist');
    LKnee = 0;
end

% RIGHT KNEE ANGLES
searchPtLabels_RKnee = strncmp(PointLabels, 'RKneeAngles',10);
if sum(searchPtLabels_RKnee == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RKnee(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RKnee = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RKnee(:,:,:) = Point_xyz(:,ncolumn_RKnee,:);
else
%     disp('Point Label "RKneeAngles" does not exist');
    RKnee = 0;
end

% LEFT SHANK ANGLES
searchPtLabels_LShank = strncmp(PointLabels, 'LShankAngles',10);
if sum(searchPtLabels_LShank == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LShank(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LShank = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LShank(:,:,:) = Point_xyz(:,ncolumn_LShank,:);
else
%     disp('Point Label "LShankAngles" does not exist');
    LShank(1,:,:) = 0;
    LShank(2,:,:) = 0;
    LShank(3,:,:) = 0;
end

% RIGHT SHANK ANGLES
searchPtLabels_RShank = strncmp(PointLabels, 'RShankAngles',10);
if sum(searchPtLabels_RShank == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RShank(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RShank = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RShank(:,:,:) = Point_xyz(:,ncolumn_RShank,:);
else
%     disp('Point Label "RShankAngles" does not exist');
    RShank(1,:,:) = 0;
    RShank(2,:,:) = 0;
    RShank(3,:,:) = 0;
end

% LEFT ANKLE ANGLES
searchPtLabels_LAnk = strncmp(PointLabels, 'LAnkleAngles',10);
if sum(searchPtLabels_LAnk == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LAnk(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LAnk = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LAnk(:,:,:) = Point_xyz(:,ncolumn_LAnk,:);
else
%     disp('Point Label "LAnkleAngles" does not exist');
    LAnk = 0;
end


% RIGHT ANKLE ANGLES
searchPtLabels_RAnk = strncmp(PointLabels, 'RAnkleAngles',10);
if sum(searchPtLabels_RAnk == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RAnk(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RAnk = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RAnk(:,:,:) = Point_xyz(:,ncolumn_RAnk,:);
else
%     disp('Point Label "RAnkleAngles" does not exist');
    RAnk = 0;
end



% LEFT FOOT PROGRESSION ANGLES
searchPtLabels_LFootP = strncmp(PointLabels, 'LFootProgressAngles',10);
if sum(searchPtLabels_LFootP == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LFootP(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LFootP = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LFootP(:,:,:) = Point_xyz(:,ncolumn_LFootP,:);
    
elseif sum(searchPtLabels_LFootP == 0)
    % Try looking for older version of Gaia.mod's foot progression angles
    searchPtLabels_LFootP = strncmp(PointLabels, 'LFootProgressAngles',10);
    if sum(searchPtLabels_LFootP == 1)
        for c1 = 1:length(PointLabels) 
            m1 = 1;
            while m1 ~= 0        
                if searchPtLabels_LFootP(1,c1) == 1
                    % finds the column number where desired point is stored
                    ncolumn_LFootP = c1;
                    m1 = 0;
                else
                    m1 = 0;
                end
            end
        end
        LFootP(:,:,:) = Point_xyz(:,ncolumn_LFootP,:);
    end    
else
%     disp('Point Label "LFootProgessAngles" does not exist');
    LFootP = 0;
end


% RIGHT FOOT PROGRESSION ANGLES
searchPtLabels_RFootP = strncmp(PointLabels, 'RFootProgressAngles',10);
if sum(searchPtLabels_RFootP == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RFootP(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RFootP = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RFootP(:,:,:) = Point_xyz(:,ncolumn_RFootP,:);
    
elseif sum(searchPtLabels_RFootP == 0)
    % Try looking for older version of Gaia.mod's foot progression angles
    searchPtLabels_RFootP = strncmp(PointLabels, 'RFootProgAngle',10);
    if sum(searchPtLabels_RFootP == 1)
        for c1 = 1:length(PointLabels) 
            m1 = 1;
            while m1 ~= 0        
                if searchPtLabels_RFootP(1,c1) == 1
                    % finds the column number where desired point is stored
                    ncolumn_RFootP = c1;
                    m1 = 0;
                else
                    m1 = 0;
                end
            end
        end
        RFootP(:,:,:) = Point_xyz(:,ncolumn_RFootP,:);
    end
    
else
%     disp('Point Label "RFootProgessAngles" does not exist');
    RFootP = 0;
end


%% KINETICS

% LEFT HIP FORCE
searchPtLabels_LHipF = strncmp(PointLabels, 'LHipForce',10);
if sum(searchPtLabels_LHipF == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LHipF(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LHipF = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LHipForce(:,:,:) = Point_xyz(:,ncolumn_LHipF,:);
else

    LHipForce = 0;
end

% RIGHT HIP FORCE
searchPtLabels_RHipF = strncmp(PointLabels, 'RHipForce',10);
if sum(searchPtLabels_RHipF == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RHipF(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RHipF = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RHipForce(:,:,:) = Point_xyz(:,ncolumn_RHipF,:);
else

    RHipForce = 0;
end


% LEFT KNEE FORCE
searchPtLabels_LKneeF = strncmp(PointLabels, 'LKneeForce',10);
if sum(searchPtLabels_LKneeF == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LKneeF(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LKneeF = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LKneeForce(:,:,:) = Point_xyz(:,ncolumn_LKneeF,:);
else

    LKneeForce = 0;
end

% RIGHT KNEE FORCE
searchPtLabels_RKneeF = strncmp(PointLabels, 'RKneeForce',10);
if sum(searchPtLabels_RKneeF == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RKneeF(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RKneeF = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RKneeForce(:,:,:) = Point_xyz(:,ncolumn_RKneeF,:);
else

    RKneeForce = 0;
end


% LEFT ANKLE FORCE
searchPtLabels_LAnkleF = strncmp(PointLabels, 'LAnkleForce',10);
if sum(searchPtLabels_LAnkleF == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LAnkleF(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LAnkleF = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LAnkleForce(:,:,:) = Point_xyz(:,ncolumn_LAnkleF,:);
else

    LAnkleForce = 0;
end

% RIGHT ANKLE FORCE
searchPtLabels_RAnkleF = strncmp(PointLabels, 'RAnkleForce',10);
if sum(searchPtLabels_RAnkleF == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RAnkleF(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RAnkleF = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RAnkleForce(:,:,:) = Point_xyz(:,ncolumn_RAnkleF,:);
else

    RAnkleForce = 0;
end


% LEFT HIP MOMENT
searchPtLabels_LHipM = strncmp(PointLabels, 'LHipMoment',10);
if sum(searchPtLabels_LHipM == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LHipM(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LHipM = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LHipMoment(:,:,:) = Point_xyz(:,ncolumn_LHipM,:);
else

    LHipMoment = 0;
end

% RIGHT HIP MOMENT
searchPtLabels_RHipM = strncmp(PointLabels, 'RHipMoment',10);
if sum(searchPtLabels_RHipM == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RHipM(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RHipM = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RHipMoment(:,:,:) = Point_xyz(:,ncolumn_RHipM,:);
else

    RHipMoment = 0;
end


% LEFT KNEE MOMENT
searchPtLabels_LKneeM = strncmp(PointLabels, 'LKneeMoment',10);
if sum(searchPtLabels_LKneeM == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LKneeM(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LKneeM = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LKneeMoment(:,:,:) = Point_xyz(:,ncolumn_LKneeM,:);
else

    LKneeMoment = 0;
end


% RIGHT KNEE MOMENT
searchPtLabels_RKneeM = strncmp(PointLabels, 'RKneeMoment',10);
if sum(searchPtLabels_RKneeM == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RKneeM(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RKneeM = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RKneeMoment(:,:,:) = Point_xyz(:,ncolumn_RKneeM,:);
else

    RKneeMoment = 0;
end


% LEFT ANKLE MOMENT
searchPtLabels_LAnkleM = strncmp(PointLabels, 'LAnkleMoment',10);
if sum(searchPtLabels_LAnkleM == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LAnkleM(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LAnkleM = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LAnkleMoment(:,:,:) = Point_xyz(:,ncolumn_LAnkleM,:);
else

    LAnkleMoment = 0;
end

% RIGHT ANKLE MOMENT
searchPtLabels_RAnkleM = strncmp(PointLabels, 'RAnkleMoment',10);
if sum(searchPtLabels_RAnkleM == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RAnkleM(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RAnkleM = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RAnkleMoment(:,:,:) = Point_xyz(:,ncolumn_RAnkleM,:);
else

    RAnkleMoment = 0;
end


% LEFT HIP POWER
searchPtLabels_LHipP = strncmp(PointLabels, 'LHipPower',9);
if sum(searchPtLabels_LHipP == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LHipP(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LHipP = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LHipPower(:,:,:) = Point_xyz(:,ncolumn_LHipP,:);
else

    LHipPower = 0;
end

% RIGHT HIP POWER
searchPtLabels_RHipP = strncmp(PointLabels, 'RHipPower',9);
if sum(searchPtLabels_RHipP == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RHipP(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RHipP = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RHipPower(:,:,:) = Point_xyz(:,ncolumn_RHipP,:);
else

    RHipPower = 0;
end


% LEFT KNEE POWER
searchPtLabels_LKneeP = strncmp(PointLabels, 'LKneePower',10);
if sum(searchPtLabels_LKneeP == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LKneeP(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LKneeP = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LKneePower(:,:,:) = Point_xyz(:,ncolumn_LKneeP,:);
else

    LKneePower = 0;
end

% RIGHT KNEE POWER
searchPtLabels_RKneeP = strncmp(PointLabels, 'RKneePower',10);
if sum(searchPtLabels_RKneeP == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RKneeP(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RKneeP = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RKneePower(:,:,:) = Point_xyz(:,ncolumn_RKneeP,:);
else

    RKneePower = 0;
end


% LEFT ANKLE POWER
searchPtLabels_LAnkleP = strncmp(PointLabels, 'LAnklePower',10);
if sum(searchPtLabels_LAnkleP == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_LAnkleP(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_LAnkleP = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    LAnklePower(:,:,:) = Point_xyz(:,ncolumn_LAnkleP,:);
else

    LAnklePower = 0;
end

% RIGHT ANKLE POWER
searchPtLabels_RAnkleP = strncmp(PointLabels, 'RAnklePower',10);
if sum(searchPtLabels_RAnkleP == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_RAnkleP(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_RAnkleP = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    RAnklePower(:,:,:) = Point_xyz(:,ncolumn_RAnkleP,:);
else

    RAnklePower = 0;
end


%% Store sorted kinematic angles and kinetics with labels in structure array
Kinematics(1).name = 'LPelv';
Kinematics(1).data = LPelv;
Kinematics(2).name = 'RPelv';
Kinematics(2).data = RPelv;

Kinematics(3).name = 'LHip';
Kinematics(3).data = LHip;
Kinematics(4).name = 'RHip';
Kinematics(4).data = RHip;

Kinematics(5).name = 'LKnee';
Kinematics(5).data = LKnee;
Kinematics(6).name = 'RKnee';
Kinematics(6).data = RKnee;

Kinematics(7).name = 'LAnk';
Kinematics(7).data = LAnk;
Kinematics(8).name = 'RAnk';
Kinematics(8).data = RAnk;

Kinematics(9).name = 'LFootP';
Kinematics(9).data = LFootP;
Kinematics(10).name = 'RFootP';
Kinematics(10).data = RFootP;

Kinematics(11).name = 'LShank';
Kinematics(11).data = LShank;
Kinematics(12).name = 'RShank';
Kinematics(12).data = RShank;

Kinematics(13).name = 'LThorax';
Kinematics(13).data = LThorax;
Kinematics(14).name = 'RThorax';
Kinematics(14).data = RThorax;
Kinematics(15).name = 'LSpine';
Kinematics(15).data = LSpine;
Kinematics(16).name = 'RSpine';
Kinematics(16).data = RSpine;


Kinetics(1).name = 'LHipF';
Kinetics(1).data = LHipForce;
Kinetics(2).name = 'RHipF';
Kinetics(2).data = RHipForce;
Kinetics(3).name = 'LKneeF';
Kinetics(3).data = LKneeForce;
Kinetics(4).name = 'RKneeF';
Kinetics(4).data = RKneeForce;
Kinetics(5).name = 'LAnkleF';
Kinetics(5).data = LAnkleForce;
Kinetics(6).name = 'RAnkleF';
Kinetics(6).data = RAnkleForce;

Kinetics(7).name = 'LHipM';
Kinetics(7).data = LHipMoment;
Kinetics(8).name = 'RHipM';
Kinetics(8).data = RHipMoment;
Kinetics(9).name = 'LKneeM';
Kinetics(9).data = LKneeMoment;
Kinetics(10).name = 'RKneeM';
Kinetics(10).data = RKneeMoment;
Kinetics(11).name = 'LAnkleM';
Kinetics(11).data = LAnkleMoment;
Kinetics(12).name = 'RAnkleM';
Kinetics(12).data = RAnkleMoment;

Kinetics(13).name = 'LHipP';
Kinetics(13).data = LHipPower;
Kinetics(14).name = 'RHipP';
Kinetics(14).data = RHipPower;
Kinetics(15).name = 'LKneeP';
Kinetics(15).data = LKneePower;
Kinetics(16).name = 'RKneeP';
Kinetics(16).data = RKneePower;
Kinetics(17).name = 'LAnkleP';
Kinetics(17).data = LAnklePower;
Kinetics(18).name = 'RAnkleP';
Kinetics(18).data = RAnklePower;
