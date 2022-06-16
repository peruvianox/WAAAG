function [] = CreateGCD(GCDfilename, typegcd, gcd_data)
% -------------------------------------------------------------------------
% CREATE GCD FILE
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: August 10, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This function creates a Polygon compatible .gcd file.
%               For use with GDI Control program.
% 
% Input:        GCDfilename     name of gcd file
%               typegcd         type of gcd file to be created
%               gcd_data        data for new gcd file
%               
%
% Output:       .gcd            creates a .gcd file of input data
%       
%


%% DEBUG LOOP
% GCDfilename = 'Ave_Rotations';
% typegcd = 'ave_rot';


%% Saves file as a .gcd
% file's extension
file_ext = '.gcd';

% Control file's name
file_name = [GCDfilename file_ext];

% Delete the contents of an existing file or create new file, and open it
% for reading and writing
fid = fopen(file_name, 'wt');


%% File Header & Data
if strcmp(typegcd, 'ave_rot') == 1
    % Must have this header to be compatible in Polygon
    fprintf(fid,'%s\n', '#!DST-0.GCD     Oxford  Metrics');
    fprintf(fid,'%s\n', '$REFERENCE');
    fprintf(fid,'%s\n', 'AVERAGE:1:4');
    
    % File data
    % 1 column, average of rotation for entire gait cycle
    fprintf(fid,'%s\n', '!LeftPelvicRotation');
    LPelvicRot = gcd_data(1).data;
    fprintf(fid,'%3.10f\n', LPelvicRot(:,:));

    fprintf(fid,'%s\n', '!LeftHipRotation');
    LHipRot = gcd_data(3).data;
    fprintf(fid,'%3.10f\n', LHipRot(:,:));

    fprintf(fid,'%s\n', '!LeftKneeRotation');
    LKneeRot = gcd_data(5).data;
    fprintf(fid,'%3.10f\n', LKneeRot(:,:));

    fprintf(fid,'%s\n', '!LeftFootRotation');
    LAnkRot = gcd_data(7).data;
    fprintf(fid,'%3.10f\n', LAnkRot(:,:));

    fprintf(fid,'%s\n', '!RightPelvicRotation');
    RPelvicRot = gcd_data(2).data;
    fprintf(fid,'%3.10f\n', RPelvicRot(:,:));

    fprintf(fid,'%s\n', '!RightHipRotation');
    RHipRot = gcd_data(4).data;
    fprintf(fid,'%3.10f\n', RHipRot(:,:));

    fprintf(fid,'%s\n', '!RightKneeRotation');
    LKneeRot = gcd_data(6).data;
    fprintf(fid,'%3.10f\n', LKneeRot(:,:));

    fprintf(fid,'%s\n', '!RightFootRotation');
    RAnkRot = gcd_data(8).data;
    fprintf(fid,'%3.10f\n', RAnkRot(:,:));

    fprintf(fid,'%s\n', '!LeftFootProgression');
    LFootProg = gcd_data(9).data;
    fprintf(fid,'%3.10f\n', LFootProg(:,:));
    
    fprintf(fid,'%s\n', '!RightFootProgression');
    RFootProg = gcd_data(10).data;
    fprintf(fid,'%3.10f\n', RFootProg(:,:));

    
elseif strcmp(typegcd, 'control') == 1
    % Must have this header to be compatible in Polygon
    fprintf(fid,'%s\n', '#!DST-0.GCD     Oxford  Metrics');
    fprintf(fid,'%s\n', '$PREFIXES');
    fprintf(fid,'%s\n', 'Left');
    fprintf(fid,'%s\n', 'Right');
    
    % File data
    % 1st column: Mean Values
    % 2nd column: Standard Deviation
    fprintf(fid,'%s\n', '!PelvicTilt 3');
    PelvicTilt = gcd_data(1).data;
    fprintf(fid,'%3.4f %3.4f\n', PelvicTilt(:,:));

    fprintf(fid,'%s\n', '!PelvicObliquity 3');
    PelvOblq = gcd_data(2).data;
    fprintf(fid,'%3.4f %3.4f\n', PelvOblq(:,:));

    fprintf(fid,'%s\n', '!PelvicRotation 3');
    PelvRot = gcd_data(3).data;
    fprintf(fid,'%3.4f %3.4f\n', PelvRot(:,:));

    fprintf(fid,'%s\n', '!HipFlexExt 3');
    HipFE = gcd_data(4).data;
    fprintf(fid,'%3.4f %3.4f\n', HipFE(:,:));

    fprintf(fid,'%s\n', '!HipAbAdduct 3');
    HipAbAd = gcd_data(5).data;
    fprintf(fid,'%3.4f %3.4f\n', HipAbAd(:,:));

    fprintf(fid,'%s\n', '!HipRotation 3');
    HipRot = gcd_data(6).data;
    fprintf(fid,'%3.4f %3.4f\n', HipRot(:,:));

    fprintf(fid,'%s\n', '!KneeFlexExt 3');
    KneeFE = gcd_data(7).data;
    fprintf(fid,'%3.4f %3.4f\n', KneeFE(:,:));

    fprintf(fid,'%s\n', '!DorsiPlanFlex 3');
    AnkleDP = gcd_data(8).data;
    fprintf(fid,'%3.4f %3.4f\n', AnkleDP(:,:));

    fprintf(fid,'%s\n', '!FootProgression 3');
    FootProg = gcd_data(9).data;
    fprintf(fid,'%3.4f %3.4f\n', FootProg(:,:));

end


fclose(fid);
