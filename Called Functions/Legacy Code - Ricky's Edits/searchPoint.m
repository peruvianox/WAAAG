function [pt_xyz] = searchPoint(Point_xyz, PointLabels, pt_label)

% -------------------------------------------------------------------------
% SEARCHES FOR A POINT
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: June 4, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%

% DESCRIPTION:  Searches for the xyz coordinates of a point with provided point
%               label. 

% INPUTS: 
%           Point_xyz          xyz coordinates of points from C3D file
%           PointLabels        labels of all points in C3D file
%           pt_label           point to find

% OUTPUTS: 
%           pt_xyz             point's x,y,z global coordinate values
%
%
%% Find and organize Point_xyz and PointLabels 

% MARKER
searchPtLabels_pt = strcmp(PointLabels, pt_label);
if sum(searchPtLabels_pt == 1)
    for c1 = 1:length(PointLabels) 
        m1 = 1;
        while m1 ~= 0        
            if searchPtLabels_pt(1,c1) == 1
                % finds the column number where desired point is stored
                ncolumn_pt = c1;
                m1 = 0;
            else
                m1 = 0;
            end
        end
    end
    pt_xyz(:,:,:) = Point_xyz(:,ncolumn_pt,:);
    else
    pt_xyz = 0;
end
