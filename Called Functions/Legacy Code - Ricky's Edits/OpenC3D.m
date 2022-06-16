function [C3Dfile] = OpenC3D(FullFileName)

% -------------------------------------------------------------------------
% Opens C3D file
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: Mar 17, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% DESCRIPTION:   When called, this function will open a gui for which the
%                user can navigate to and select (only) a C3D file.
%
% INPUTS: 
%           C3Dfile                 User specified C3D file to be analyzed
%
% OUTPUTS:  
%           C3Dfile                 User specified C3D file to be analyzed
%           FullFileName            C3D file's full file path


%% OPENS C3D File

file_level = findstr(FullFileName, '\');
if file_level > 0
    FileName = FullFileName(file_level(length(file_level))+1:length(FullFileName));
else
    FileName = FullFileName;
end


% C3Dfile = fopen(FILENAME, PERMISSION, MACHINEFORMAT)
C3Dfile = fopen(FullFileName,'r','n');

if C3Dfile == -1,
    % If fopen cannot open the file it returns C3Dfile == -1 and an error dialog
    % will appear and remain until user clears error dialog box.
    message = errordlg(['File: ',FileName,' could not be opened'],'application error');
    uiwait(message)
    return
end


NBlockFirstParamBlock = fread(C3Dfile, 1, 'int8');

% for all C3D files IDbyte = 80
IDbyte = fread(C3Dfile, 1, 'int8');

if IDbyte ~= 80
    message = errordlg(['File: ',FileName,' does not comply with the C3D format'],'application error');
    uiwait(message)
%     fclose(C3Dfile)
    return
end

%% Add Working Directory Retention?

