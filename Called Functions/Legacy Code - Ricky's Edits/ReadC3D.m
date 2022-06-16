function [HeaderInfo, Point_xyz, PointLabels, GaitEvents] = ReadC3D(C3Dfile)
% -------------------------------------------------------------------------
% Reads and Extracts Data from C3D file
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: Aug 19, 2008
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% DESCRIPTION:   When called, this function will obtain the output variables (listed
%                below) for the user specified C3D file (FullFileName).
%
% INPUTS: 
%           C3Dfile                 User specified C3D file to be analyzed

% OUTPUTS:  
%           HeaderInfo              Information from Header section of file
%           Point_xyz               3D coordinates for all markers in file
%           PointLabels             Labels for all 3D data points in file
%           GaitEvents              Contains the event type, side, and time
%                                   of occurance
%           AnalogData              Each channel of analog data
%
%
% ACKNOWLEDGMENTS:  The following codes for C3D file reading were used as
%                   references to create this function.
%
% Ver. 1.0 Creation (Alan Morris, Toronto, October 1998) [originally named "getc3d.m"]
% Ver. 2.0 Revision (Jaap Harlaar, Amsterdam, april 2002)
%
% Modified by May Liu, Dec 2004. Added the HeaderGroup, timeVector, and
% changed some of the type parameters (e.g., int8, float32, etc)
%
% Modified by Michael Schwartz, Dec 2004. Changed data reading; eliminated
% loops (read in blocks with skip) - dramatically faster (150x - 500x).
%
% Documentation about C3D files from C3D.org


% Point_xyz = [];
% VideoFrameRate = 0;
% AnalogData = [];
% AnalogFrameRate = 0;
% Event = [];
% ParameterGroup = [];
% CameraInfo = [];
% ResidualError = [];
% HeaderGroup = [];

%% Reads user specified C3D file
NBlockFirstParamBlock = fread(C3Dfile, 1, 'int8');

% for all C3D files IDbyte = 80
IDbyte = fread(C3Dfile, 1, 'int8');



%%              READS HEADER SECTION

% Determine the type of processor by jumping to this field
fseek(C3Dfile, 512*(NBlockFirstParamBlock-1)+3, 'bof');

% processor types: 1(INTEL-PC), 2(DEC-VAX), 3(MIPS-SUN/SGI)
processor_type = fread(C3Dfile, 1, 'int8')-83;

if processor_type == 2,
    fclose(C3Dfile);
    
    % DEC VAX D floating point and VAX ordering
    C3Dfile = fopen(FullFileName,'r','d');
end

% Resets the file position indicator to the beginning of the c3dfile
fseek(C3Dfile, 0, 'bof');

% fread(C3Dfile, size, precision)
% reads the c3dfile
% size reads the specified number of elements into a column vector (i.e. 1) 
% precision is the datatype specifier

% Word 1, Byte 1: first byte of C3Dfile points to the first block of the parameter section
N_ParamBlockStart = fread(C3Dfile, 1, 'int8');

% Word 1, Byte2: for all C3D files IDbyte = 80
IDbyte = fread(C3Dfile, 1, 'int8');

% Word 2: number of 3D trajectories in C3D file
N_3Dtrajects = fread(C3Dfile, 1, 'int16');

% Word 3:total number of analaog measurements per 3D frame
% (i.e. number of channels multiplied by the samples per channel)
N_AnalogMeasures_per3Dframe = fread(C3Dfile, 1, 'int16');

% Word 4: number of the first frame of 3D data
FirstFrame = fread(C3Dfile, 1, 'int16');

% Word 5: number of the last frame of 3D data
LastFrame = fread(C3Dfile, 1, 'int16');

% Word 6: maximum interpolation gap allowed in 3D frames
MaxInterpGap = fread(C3Dfile, 1, 'int16');

% Words 7-8: the 3D scale factor (floating-point) that converts signed integer 
% 3D data to reference system measurement units.  If this is negative then 
% the file is scaled in floating-point.
Scale = fread(C3Dfile, 1, 'float32');

% if the 3D scale factor (Scale) value is positive the data type is integer
if Scale < 0
    DataTypeFloat = true;
else
    DataTypeFloat = false;
end

% Word 9: the number of the first block of the 3D and analog data section
N_DataBlockStart = fread(C3Dfile, 1, 'int16');

% Word 10: number of analog samples per 3D frame
Analog_samps_3Dframe = fread(C3Dfile, 1, 'int16');

% Words 11-12: the 3D frame rate in Hz (floating point)
% aka POINT:RATE
FrameRate3D = fread(C3Dfile, 1, 'float32');

if N_AnalogMeasures_per3Dframe > 0
    % Number of analog channels sampled
    % ANALOG:USED
    N_AnalogChannels = N_AnalogMeasures_per3Dframe/Analog_samps_3Dframe;
else
    N_AnalogChannels = 0;
end

% Calculates the analog frame rate (Hz)
% aka ANALOG:RATE
AnalogFrameRate = FrameRate3D*Analog_samps_3Dframe;

% AnalogFrameRate = 2x's the Nyquist Frequency
Nyquist_Freq = AnalogFrameRate/2;


% Reads header's events
% places pointer before 150th word (bytes 299 and 300)
fseek(C3Dfile, 298, 'bof');

% word 150: key value (12345 decimal) will indicate if the
% C3D file supports labels with 4 char, doesn't indicate if data
% is actually stored
EventIndicator = fread(C3Dfile, 1, 'int16');

if EventIndicator == 12345,
    % Word 151: number of events in the c3d file's header
    N_events = fread(C3Dfile, 1, 'int16');
    
    % skips one position (2 bytes)
    fseek(C3Dfile, 2, 'cof');
    
    if N_events > 0
        % Builds a multidimensional array (Event) with three subarrays 
        % time, value, and name
        for n = 1:N_events
            % words 153-188 stores event times
            Event(n).time = fread(C3Dfile, 1, 'float');
        end

        fseek(C3Dfile, 188*2, 'bof');
        for n = 1:N_events
            % words 189-197 stores event values: 0x00 = 0N, 0x01 = OFF
            Event(n).value = fread(C3Dfile, 1, 'int8');
        end
        
        fseek(C3Dfile, 198*2, 'bof');
        for n = 1:N_events
            % words199-234 stores event names (4char)
            Event(n).name = cellstr(char(fread(C3Dfile, 4, 'char')'));
        end
    end

end

% Organizes header information into a structure array
HeaderInfo(1).name = 'First Parameter';
HeaderInfo(1).data = N_ParamBlockStart;
HeaderInfo(2).name = 'Number of Markers';
HeaderInfo(2).data = N_3Dtrajects;
HeaderInfo(3).name = 'Number of Analog Channels';
HeaderInfo(3).data = N_AnalogChannels;
HeaderInfo(4).name = 'First Frame';
HeaderInfo(4).data = FirstFrame;
HeaderInfo(5).name = 'Last Frame';
HeaderInfo(5).data = LastFrame;
HeaderInfo(6).name = 'Video Sampling Rate';
HeaderInfo(6).data = FrameRate3D;
HeaderInfo(7).name = 'Analog Sampling Rate (Hz)';
HeaderInfo(7).data = AnalogFrameRate;
HeaderInfo(8).name = 'Scale Factor';
HeaderInfo(8).data = Scale;
HeaderInfo(9).name = 'Start of Data Section';
HeaderInfo(9).data = N_DataBlockStart;
HeaderInfo(10).name = 'Maximum Interpolation Gap (mm)';
HeaderInfo(10).data = MaxInterpGap;
HeaderInfo(11).name = 'Number of Analog Samples per 3D Frame';
HeaderInfo(11).data = Analog_samps_3Dframe;
HeaderInfo(12).name = 'Nyquist Frequency (Hz)';
HeaderInfo(12).data = Nyquist_Freq;


%%               READS PARAMETER SECTION                          

% repositions the file position indicator 
fseek(C3Dfile, 512*(N_ParamBlockStart-1), 'bof');

% not sure what this is... *********************************************
dat1 = fread(C3Dfile, 1, 'int8');

% Identification byte, key = 80 if file is in C3D format
IDbyte2 = fread(C3Dfile, 1, 'int8');

% number of parameter blocks to follow (byte 3 in parameter section)
Total_number_paramblocks = fread(C3Dfile, 1, 'uint8');

% processor types: 1(INTEL-PC), 2(DEC-VAX), 3(MIPS-SUN/SGI)
% byte 4 = 83+processor type
processor_type = fread(C3Dfile, 1, 'int8')-83;

% number of characters in Group:ParamName
N_characters = fread(C3Dfile, 1, 'int8');

% ID number, where -value = group, +value = parameter
GroupNumber = fread(C3Dfile, 1, 'int8');

while N_characters > 0
    % End of the parameter record indicated by <0 characters 
    % for Group:ParamName
    
    % Organizes data into structural array 'ParameterGroup' with a name, description, and
    % data

    if GroupNumber < 0
        %%% Reads and organizes Group Data %%%
        GroupNumber = abs(GroupNumber);
        GroupName = fread(C3Dfile, [1, N_characters], 'char');
        
        % Group name
        ParameterGroup(GroupNumber).name = cellstr(char(GroupName));

        % outputs the current file position
        fileposition = ftell(C3Dfile);
        
        % offset is in bytes
        offset = fread(C3Dfile,1,'int16');
        
        % next_section is the position of the beginning of the next section
        next_section = fileposition + offset;

        % Group description
        descript_chars = fread(C3Dfile,1,'int8');
        GroupDescript = fread(C3Dfile,[1,descript_chars],'char');
        ParameterGroup(GroupNumber).description = cellstr(char(GroupDescript));
        ParamNumberIndex(GroupNumber) = 0;

        % positions pointer to next Group
        fseek(C3Dfile,next_section,'bof');

    
    else
        %%% Reads and organizes Parameter Data %%%
        clear dimension;
        
        ParamNumberIndex(GroupNumber) = ParamNumberIndex(GroupNumber) + 1;
        
        % indexes all paramaters within in group
        ParamNumber = ParamNumberIndex(GroupNumber);

        % assigns the parameter's name to array 'ParamName'
        ParamName = fread(C3Dfile, [1, N_characters], 'char');

        % reads and saves parameter's name
        if size(ParamName) > 0
            ParameterGroup(GroupNumber).Parameter(ParamNumber).name = cellstr(char(ParamName));
        end

        % outputs the current file position
        fileposition = ftell(C3Dfile);
        
        % offset is in bytes
        offset = fread(C3Dfile,1,'int16');
        
        % next_section is the position of the beginning of the next section
        next_section = fileposition + offset;

        % determines data type: -1=char, 1=byte, 2=integer*2, 4=real*4
        type = fread(C3Dfile,1,'int8');
        ParameterGroup(GroupNumber).Parameter(ParamNumber).datatype = type;

        % reads number of dimensions
        N_dimensions = fread(C3Dfile,1,'int8');

        if N_dimensions == 0
            % datalength = length of data section
            datalength = abs(type);
        else
            mult = 1;
            for k = 1:N_dimensions
                % change to "uint8" (was previously int8)
                dimension(k) = fread(C3Dfile,1,'uint8');
                mult = mult*dimension(k);
                
                % saves parameter dimension data
                ParameterGroup(GroupNumber).Parameter(ParamNumber).dim(k) = dimension(k);
            end
            
            % length of data section for multi-dimensional array
            datalength = abs(type)*mult;
        end
   
    

        if type == -1
            % for datatype == char
            
            % length of character word
            wordlength = dimension(1);

            if N_dimensions == 2 & datalength>0

                for j = 1:dimension(2)
                    % organizes character word data section into a 2D array
                    data = fread(C3Dfile,[1,wordlength],'char');
                    ParameterGroup(GroupNumber).Parameter(ParamNumber).data(j)...
                        = cellstr(char(data));
                end

            elseif N_dimensions==1 & datalength>0
                % organizes numerical data section into a 1D array
                data = fread(C3Dfile,[1,wordlength],'char');
                ParameterGroup(GroupNumber).Parameter(ParamNumber).data...
                    = cellstr(char(data));
            end

        elseif type == 1
            % for datatype == 'boolean'
            
            % Number of parameters
            Nparameters = datalength/abs(type);
            
            data = fread(C3Dfile,Nparameters,'int8');
            ParameterGroup(GroupNumber).Parameter(ParamNumber).data = data;

         
        elseif type == 2 & datalength > 0
            % for datatype == int
            
            Nparameters = datalength/abs(type);
            data = fread(C3Dfile,Nparameters,'int16');

            if N_dimensions > 1
                ParameterGroup(GroupNumber).Parameter(ParamNumber).data...
                    = reshape(data,dimension);
            else
                ParameterGroup(GroupNumber).Parameter(ParamNumber).data...
                    = data;
            end


        elseif type == 4 & datalength > 0
            % for datatype == float

            Nparameters = datalength/abs(type);
            data = fread(C3Dfile,Nparameters,'float');

            if N_dimensions > 1
                ParameterGroup(GroupNumber).Parameter(ParamNumber).data...
                    = reshape(data,dimension);
            else
                ParameterGroup(GroupNumber).Parameter(ParamNumber).data...
                    = data;
            end

        end

        % reads descriptions
        descript_chars = fread(C3Dfile,1,'int8');	

        if descript_chars > 0
            description = fread(C3Dfile,[1,descript_chars],'char');
            ParameterGroup(GroupNumber).Parameter(ParamNumber).description...
                = cellstr(char(description));
        end
        
        % moves pointer to next section
        fseek(C3Dfile,next_section,'bof');

    end

    % checks group:parameter characters and ID number to see if more
    % sections exist
    
    % number of characters in the next group:parameter name
    N_characters = fread(C3Dfile,1,'int8');
    % ID number, -ve = group, +ve = parameter
    GroupNumber = fread(C3Dfile,1,'int8');
end


%%            READS 3D DATA SECTION         
% Places the file position indicator past the single 512byte header section
% and at the 1st byte of the 3D/Analog data section
fseek(C3Dfile,(N_DataBlockStart-1)*512,'bof');
N_VidFrames = LastFrame - FirstFrame + 1;

if DataTypeFloat
    reps = [int2str(4*N_3Dtrajects), '*float32'];
    tempMarker = fread(C3Dfile, (4*N_3Dtrajects*N_VidFrames), reps, (4*Analog_samps_3Dframe*N_AnalogChannels));
else
    reps = [int2str(4*N_3Dtrajects), '*int16'];
    tempMarker = fread(C3Dfile, (4*N_3Dtrajects*N_VidFrames), reps, (2*Analog_samps_3Dframe*N_AnalogChannels));
end

% reshapes matrix 'tempMarker' from one column to 4xN_3DtrajectsxN_VidFrames
tempMarker = reshape(tempMarker, 4, N_3Dtrajects, N_VidFrames);

% removes camera contribution and 3D point residuals
Point_xyz = tempMarker(1:3,:,:);

% Calculates camera residual
adjMarker = squeeze(fix(tempMarker(4,:,:)));
highbyte = fix(adjMarker/256);
lowbyte = adjMarker-highbyte*256;
CameraInfo = highbyte;
ResidualError = lowbyte*abs(Scale);



%% Reads Analog Data
% N_VidFrames = LastFrame - FirstFrame + 1;
% if DataTypeFloat
%     % Set initial read to after first block of digital data (header data + 4
%     % sets of digital data * number of digital channels * number of bytes
%     %(32 bits = 4 bytes)
%     fseek(C3Dfile,(N_DataBlockStart-1)*512+4*N_3Dtrajects*4,'bof');
%     reps = [int2str(4*Analog_samps_3Dframe*N_AnalogChannels),'*float32'];
%     tempMarker = fread(C3Dfile, (4*Analog_samps_3Dframe*N_AnalogChannels*N_VidFrames), reps, 4*N_3Dtrajects*4);
% else
%     % Set initial read to after first block of digital data (header data + 4
%     % sets of digital data * number of digital channels * number of bytes
%     % (16 bits = 2 bytes)
%     fseek(C3Dfile,(N_DataBlockStart-1)*512+4*N_3Dtrajects*2,'bof');
%     reps = [int2str(Analog_samps_3Dframe*N_AnalogChannels),'*int16'];
%     tempMarker = fread(C3Dfile, (4*Analog_samps_3Dframe*N_AnalogChannels*N_VidFrames), reps, 4*N_3Dtrajects*2);
% end
% 
% % reshape matrix 'tempMarker' from one column to
% % N_AnalogChannelsxAnalog_samps_3DframxN_VidFrames (each row is set of
% % analog data for each 3D frame).
% tempMarker = reshape(tempMarker,N_AnalogChannels,Analog_samps_3Dframe,N_VidFrames);


%% Organizes 3D Data

% Find the group number for GROUP
% GroupName_group = [ParameterGroup.name];
% nameID_group = strcmp(GroupName_group,'GROUP');
% 
% % Find the parameter number for the 'LABELS' parameter
% ParamNumb_labels = [ParameterGroup(nameID_group).Parameter.name];
% numberID_grouplabels = strcmp(ParamNumb_labels,'LABELS');
% 
% % Get the labels
% GroupLabels = [ParameterGroup(nameID_group).Parameter(numberID_grouplabels).data];


%%%%% POINTS %%%%%
% Finds the group number for the 'POINT' group, stored in struct array ParameterGroup
GroupName_point = [ParameterGroup.name];
nameID_point = strcmp(GroupName_point,'POINT');

% Finds the parameter number for 'LABELS' for POINTS
GroupPoint_labels = [ParameterGroup(nameID_point).Parameter.name];
numberID_ptlabels = strcmp(GroupPoint_labels,'LABELS');

% Assigns subgroup 'LABELS' from group 'POINTS' to variable PointLabels
PointLabels = [ParameterGroup(nameID_point).Parameter(numberID_ptlabels).data];


%%%%% EVENTS %%%%%
% Find the group number for GROUP
GroupName_group = [ParameterGroup.name];
nameID_event = strcmp(GroupName_group, 'EVENT');
% Find the parameter number for EVENT:LABELS
ParamNumb_Elabels = [ParameterGroup(nameID_event).Parameter.name];
numberID_Elables = strcmp(ParamNumb_Elabels, 'LABELS');
% Get the labels
EventLabels = [ParameterGroup(nameID_event).Parameter(numberID_Elables).data];

% Find the parameter number for EVENT:CONTEXTS
ParamNumb_Econtexts = [ParameterGroup(nameID_event).Parameter.name];
numberID_Econtexts = strcmp(ParamNumb_Elabels, 'CONTEXTS');
% Get the labels
EventContexts = [ParameterGroup(nameID_event).Parameter(numberID_Econtexts).data];

% Find the parameter number for EVENT:TIMES
ParamNumb_Etimes = [ParameterGroup(nameID_event).Parameter.name];
numberID_Etimes = strcmp(ParamNumb_Etimes, 'TIMES');
% Get the times
EventTimes = [ParameterGroup(nameID_event).Parameter(numberID_Etimes).data];


% To obtain actual event times, add the 2 values together using double
% precision floating point storage and multiply by the rate (Hz) at which data
% was collected
EventTimes = (sum(EventTimes,'double'))*FrameRate3D;
EventTimes_int64 = int64(EventTimes);
EventTimes = num2cell(EventTimes_int64);

unsortedGaitEvents(:,:,1) = EventLabels;
unsortedGaitEvents(:,:,2) = EventContexts;
unsortedGaitEvents(:,:,3) = EventTimes;

% converts GaitEvents(:,:,3) into a matrix so sort function can be used
% GEtimes = cell2mat(unsortedGaitEvents(:,:,3));
% sortedGEtimes = sort(GEtimes);
% lengthGE = length(unsortedGaitEvents(:,:,3));

% converts GaitEvents(:,:,3) into a matrix so sort function can be used
GEtimes = cell2mat(unsortedGaitEvents(:,:,3));
%GEtimes = cell2mat(EventTimes);
[sortedGEtimes, index_sorted] = sort(GEtimes);

% using the index of sorted times, re-organizes gait events' (time,
% context, & label) so in chronological order of occurance
GaitEvents = unsortedGaitEvents(:,index_sorted,:);







