function [StrideMatrix, EnsAvgMatrix, STDMatrix] = C3D2Matrix (C3D, LStride,RStride)
% -------------------------------------------------------------------------
% Reorganize C3D kinematic data to 101x24 Matrix
% -------------------------------------------------------------------------
% Author: Ricky Pimentel
% Date: 12/1/2016
% For use at the Center for Gait and Movement Analysis Laboratory
% at the Children's Hospital in Denver, CO.
%
% Description: This function will reorganize C3D-structured kinematic data
% into a 101x 24 matrix with columns organized in XYZ triplets of LPelvis,
% RPelvis, LHip, RHip, LKnee, RKnee, LAnkle,RAnkle.
%
% Input:     C3D                        Structured Dataset of kinematic data
%                                               from analyze_C3D or similar 
%               LStride, RStride       strides to be analyzed (optional)
%                                              *if no strides are selected,
%                                              the representative strides
%                                              will be chosen 
%
% Output:      StrideMatrix       Data corresponding to the strides chosen
%                   EnsAvgMatrix     Ensemble average of all gait cycles in the trial      
%                   STDMatrix           Standard Deviations for the
%                                             ensemble average data
% 

%% If L and R strides not specified, use representative cycle
% LEFT
ise = evalin( 'base', 'exist(''LStride'',''var'') == 1' );
if ise == 1
    LStride = LStride;
else
    LStride = C3D(31).GCranking.left(1);
end
% RIGHT
ise = evalin( 'base', 'exist(''RStride'',''var'') == 1' );
if ise == 1
    RStride = RStride;
else
    RStride = C3D(31).GCranking.right(1);
end

%% Define Matrix of representative (or selected) strides

StrideMatrix(:,1) = C3D(1).VarPlane_GCs(1).left(:,LStride); % L Pelvis Angles
StrideMatrix(:,2) = C3D(1).VarPlane_GCs(2).left(:,LStride);
StrideMatrix(:,3) = C3D(1).VarPlane_GCs(3).left(:,LStride);
StrideMatrix(:,4) = C3D(2).VarPlane_GCs(1).right(:,RStride); % R Pelvis Angles
StrideMatrix(:,5) = C3D(2).VarPlane_GCs(2).right(:,RStride);
StrideMatrix(:,6) = C3D(2).VarPlane_GCs(3).right(:,RStride);
StrideMatrix(:,7) = C3D(3).VarPlane_GCs(1).left(:,LStride); % L Hip Angles
StrideMatrix(:,8) = C3D(3).VarPlane_GCs(2).left(:,LStride);
StrideMatrix(:,9) = C3D(3).VarPlane_GCs(3).left(:,LStride);
StrideMatrix(:,10) = C3D(4).VarPlane_GCs(1).right(:,RStride); % R Hip Anlges
StrideMatrix(:,11) = C3D(4).VarPlane_GCs(2).right(:,RStride);
StrideMatrix(:,12) = C3D(4).VarPlane_GCs(3).right(:,RStride);
StrideMatrix(:,13) = C3D(5).VarPlane_GCs(1).left(:,LStride); % L Knee Angles
StrideMatrix(:,14) = C3D(5).VarPlane_GCs(2).left(:,LStride); 
StrideMatrix(:,15) = C3D(5).VarPlane_GCs(3).left(:,LStride); 
StrideMatrix(:,16) = C3D(6).VarPlane_GCs(1).right(:,RStride); % R Knee Angles
StrideMatrix(:,17) = C3D(6).VarPlane_GCs(2).right(:,RStride);
StrideMatrix(:,18) = C3D(6).VarPlane_GCs(3).right(:,RStride);
StrideMatrix(:,19) = C3D(7).VarPlane_GCs(1).left(:,LStride); % L Ankle Angles
StrideMatrix(:,20) = C3D(9).VarPlane_GCs(2).left(:,LStride);  % L Prog Angle
StrideMatrix(:,21) = C3D(7).VarPlane_GCs(3).left(:,LStride); 
StrideMatrix(:,22) = C3D(8).VarPlane_GCs(1).right(:,RStride); % R Ankle Angles
StrideMatrix(:,23) = C3D(10).VarPlane_GCs(2).right(:,RStride); % R Prog Angle
StrideMatrix(:,24) = C3D(8).VarPlane_GCs(3).right(:,RStride);

%% Define Matrix of Ensemble Averages

EnsAvgMatrix(:,1) = C3D(1).mean_GC(1).left(:,1); % L Pelvis Angles
EnsAvgMatrix(:,2) = C3D(1).mean_GC(2).left(:,1);
EnsAvgMatrix(:,3) = C3D(1).mean_GC(3).left(:,1);
EnsAvgMatrix(:,4) = C3D(2).mean_GC(1).right(:,1); % R Pelvis Angles
EnsAvgMatrix(:,5) = C3D(2).mean_GC(2).right(:,1);
EnsAvgMatrix(:,6) = C3D(2).mean_GC(3).right(:,1);
EnsAvgMatrix(:,7) = C3D(3).mean_GC(1).left(:,1); % L Hip Angles
EnsAvgMatrix(:,8) = C3D(3).mean_GC(2).left(:,1);
EnsAvgMatrix(:,9) = C3D(3).mean_GC(3).left(:,1);
EnsAvgMatrix(:,10) = C3D(4).mean_GC(1).right(:,1); % R Hip Anlges
EnsAvgMatrix(:,11) = C3D(4).mean_GC(2).right(:,1);
EnsAvgMatrix(:,12) = C3D(4).mean_GC(3).right(:,1);
EnsAvgMatrix(:,13) = C3D(5).mean_GC(1).left(:,1); % L Knee Angles
EnsAvgMatrix(:,14) = C3D(5).mean_GC(2).left(:,1); 
EnsAvgMatrix(:,15) = C3D(5).mean_GC(3).left(:,1); 
EnsAvgMatrix(:,16) = C3D(6).mean_GC(1).right(:,1); % R Knee Angles
EnsAvgMatrix(:,17) = C3D(6).mean_GC(2).right(:,1);
EnsAvgMatrix(:,18) = C3D(6).mean_GC(3).right(:,1);
EnsAvgMatrix(:,19) = C3D(7).mean_GC(1).left(:,1); % L Ankle Angles
EnsAvgMatrix(:,20) = C3D(9).mean_GC(2).left(:,1);  % L Prog Angle
EnsAvgMatrix(:,21) = C3D(7).mean_GC(3).left(:,1); 
EnsAvgMatrix(:,22) = C3D(8).mean_GC(1).right(:,1); % R Ankle Angles
EnsAvgMatrix(:,23) = C3D(10).mean_GC(2).right(:,1); % R Prog Angle
EnsAvgMatrix(:,24) = C3D(8).mean_GC(3).right(:,1);


%% Define Matrix of Standard Deviations

STDMatrix(:,1) = C3D(1).stddev_GC(1).left(:,1); % L Pelvis Angles
STDMatrix(:,2) = C3D(1).stddev_GC(2).left(:,1);
STDMatrix(:,3) = C3D(1).stddev_GC(3).left(:,1);
STDMatrix(:,4) = C3D(2).stddev_GC(1).right(:,1); % R Pelvis Angles
STDMatrix(:,5) = C3D(2).stddev_GC(2).right(:,1);
STDMatrix(:,6) = C3D(2).stddev_GC(3).right(:,1);
STDMatrix(:,7) = C3D(3).stddev_GC(1).left(:,1); % L Hip Angles
STDMatrix(:,8) = C3D(3).stddev_GC(2).left(:,1);
STDMatrix(:,9) = C3D(3).stddev_GC(3).left(:,1);
STDMatrix(:,10) = C3D(4).stddev_GC(1).right(:,1); % R Hip Anlges
STDMatrix(:,11) = C3D(4).stddev_GC(2).right(:,1);
STDMatrix(:,12) = C3D(4).stddev_GC(3).right(:,1);
STDMatrix(:,13) = C3D(5).stddev_GC(1).left(:,1); % L Knee Angles
STDMatrix(:,14) = C3D(5).stddev_GC(2).left(:,1); 
STDMatrix(:,15) = C3D(5).stddev_GC(3).left(:,1); 
STDMatrix(:,16) = C3D(6).stddev_GC(1).right(:,1); % R Knee Angles
STDMatrix(:,17) = C3D(6).stddev_GC(2).right(:,1);
STDMatrix(:,18) = C3D(6).stddev_GC(3).right(:,1);
STDMatrix(:,19) = C3D(7).stddev_GC(1).left(:,1); % L Ankle Angles
STDMatrix(:,20) = C3D(9).stddev_GC(2).left(:,1);  % L Prog Angle
STDMatrix(:,21) = C3D(7).stddev_GC(3).left(:,1); 
STDMatrix(:,22) = C3D(8).stddev_GC(1).right(:,1); % R Ankle Angles
STDMatrix(:,23) = C3D(10).stddev_GC(2).right(:,1); % R Prog Angle
STDMatrix(:,24) = C3D(8).stddev_GC(3).right(:,1);


