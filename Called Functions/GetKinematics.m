function  [Kinematics] = GetKinematics(filename)
%%%------------------------------------------------------------------%%%
% Ricky Pimentel   2017
% Center for Gait and Movement Analysis, Children's Hospital Colorado
%
% GetKinematics
% 
% Description
% This function will load and process a C3D file that has gait events and a
% kinematic model applied to it (PiG or GAIA). Various options can be
% selected using the GetKinematics Selector_File (see Selector_File.xlsx
% for more info). 
%
% Note: May not run on older C3D files!
%
% Input (optional)
% a single filename can be called as an input to this function. If no file
% is chosen when the function is called, the user will be prompted to
% select one or multiple files for analysis. 
% 
% Output
% Kinematics, a structural array containing:
%       RepTrial            Infomation on the representative trial of the
%                               selected files (Num), the name of the C3D file that
%                               contains the representative cycle, and
%                               RepCycles, the cycle iteration that contains the 
%                               left and right representative cycle. 
%       Trials                A structural array containing all the
%                               kinematic data within the file,
%                               interpolated to 0-100% of the gait cycle.
%                               The multiple columns represent the number
%                               of gait cycles in the trial. 
%       FileName          The name of the file loaded and processed. 
%       NumCycles       The number of valid cycles on the left and right
%                               sides, respectively. 
%       RepCycle          The representative cycle kinematics. Description
%                               of data orgainization can be found in
%                               ExportedKinematics.xlsx, cells A26:Y130. 
%       EnsembleAverage     The ensemble average kinematics. Organized
%                               similar to RepCycle
%       STD                 The standard deviation of the kinematics, also
%                               organized similar to RepCycle. 
%       GDI                   The gait deviation index scores for the left
%                                and right sides. 
%
%   Note: If multiple trials are imported, only one RepTrial, RepCycle,
%   EnsembleAverage, and STD will be output. 
%
% % DEBUG LOOP %%%%%%%%%
% close all
% clear
% clc
% %%% END DEBUG LOOP %%%%%

warning off;

%% Select data for IMPORT
if exist('filename','var') == 0  % test to see if data exists, if not -> load with uigetfile below    
[files2get, ~] = uigetfile('*.c3d', 'Please select C3D file(s) to be analyzed', 'MultiSelect', 'on');
% [FullFileName, ~] = uigetfile('*.c3d', 'Please select C3D file(s) to be analyzed', 'MultiSelect', 'on');
else
    files2get = filename; 
end

% Determine # of trials selected
if iscell(files2get) == 0
    NumTrials = 1;
else
    NumTrials = length(files2get);
end

%% Load Settings Data
[KinematicsOptions,~, ~, ~] = SelectOptions;

%% Define Input Settings
% Define Markerset Type
for i = 1:length(KinematicsOptions)
    Ind = strcmp(KinematicsOptions(i).Input,'Marker Set Type') == 1;
    if Ind == 1
        MarkerSetType = char(KinematicsOptions(i).Choice);
    end
end
% Define whether to delete first cycle
for i = 1:length(KinematicsOptions)
    Ind = strcmp(KinematicsOptions(i).Input,'Delete First Gait Cycle') == 1;
    if Ind == 1
        DeleteFirstCycle = char(KinematicsOptions(i).Choice);
    end
end
% Define whether to include TemporalDistanceParameters
for i = 1:length(KinematicsOptions)
    Ind = strcmp(KinematicsOptions(i).Input,'Temporal Distance') == 1;
    if Ind == 1
        TemporalDistance = char(KinematicsOptions(i).Choice);
    end
end
% GDI RC or EA (both will be processed, only the selected will be used
% for GDI calculation.
for i = 1:length(KinematicsOptions)
    Ind = strcmp(KinematicsOptions(i).Input,'GDI RC or EA') == 1;
    if Ind == 1
        GDI_RCorEA = char(KinematicsOptions(i).Choice);
    end
end
% Define whether to include Full Foot Global Angles, rather than just foot prog
for i = 1:length(KinematicsOptions)
    Ind = strcmp(KinematicsOptions(i).Input,'Foot Global Angles') == 1;
    if Ind == 1
        FootGlobAng = char(KinematicsOptions(i).Choice);
    end
end
% % if subject is hemiplegic, pick side
% for i = 1:length(KinematicsOptions)
%     Ind = strcmp(KinematicsOptions(i).Input,'Hemiplegic?') == 1;
%     if Ind == 1
%         Hemi = char(KinematicsOptions(i).Choice);
%     end
% end

%%  Determine Representative Trials and/or Cycles
[RepTrial, GDI] =  DetermineRepTrialCyc(files2get);
Kinematics.RepTrial = RepTrial;
% Kinematics.RepTrial.RepCycles = RepCycles (1,:);

%% Import C3D Data
for z = 1:NumTrials
    % Determine filename(s)
    if NumTrials == 1
        FullFileName = files2get;
    else
        FullFileName = char(files2get{z});
    end
    
    %% Opens and Extract C3D file's data
    [~, ~, ~, ~, kinematics, TempSpat,~] = GDI_AnalyTrial(FullFileName);

    % Save Temporal Spatial Data
    if strcmp(TemporalDistance,'Yes') ==1
        Kinematics(z).TempSpat = TempSpat;
    end
    % Save Kinematic Data
    Kinematics(z).Trials = kinematics;
    
    % Re-organize data if not in correct order
    %if strcmp(FootGlobAng,'No') == 1
        if strcmp(Kinematics(z).Trials(8).name,'AnkleDP') ==1
            % define un-organized data
            AnkleDP = Kinematics(z).Trials(8).data;
            FootP = Kinematics(z).Trials(9).data;
            KneeABAD = Kinematics(z).Trials(10).data;
            KneeR = Kinematics(z).Trials(11).data;
            FootR = Kinematics(z).Trials(12).data;
            % re-organize data
            Kinematics(z).Trials(8).name = 'KneeAA';
            Kinematics(z).Trials(8).data = KneeABAD;
            Kinematics(z).Trials(9).name = 'KneeR';
            Kinematics(z).Trials(9).data = KneeR;
            Kinematics(z).Trials(10).name = 'AnkleDP';
            Kinematics(z).Trials(10).data = AnkleDP;
            Kinematics(z).Trials(11).name = 'FootP';
            Kinematics(z).Trials(11).data = FootP;
            Kinematics(z).Trials(12).name = 'AnkleR';
            Kinematics(z).Trials(12).data = FootR;
        end
end

%%  save the file name(s) of imported data
if NumTrials == 1
    Kinematics(z).FileName = FullFileName;
else
    for z = 1:NumTrials
        Kinematics(z).FileName = files2get{z};
    end
end

%% Find Gait Cycles
% Find number of cycles for each leg for each trial
for z = 1:length(Kinematics)
    % LEFT
    [~,NL(z)] = size(Kinematics(z).Trials(1).data.left);
    NumCyclesL = NL;
    % RIGHT
    [~,NR(z)] = size(Kinematics(z).Trials(1).data.right);
    NumCyclesR = NR;
end
clearvars NL NR
% Find total number of cycles for all trials
TotCyclesL = sum(NumCyclesL);
TotCyclesR = sum(NumCyclesR);

% Delete First cycle?
% Option to delete the first gait cycle on each side. Gait cycles will
% still be present in the Kinematics.Trials.data, but the "deleted" first
% cycles will not be used in calculation of ensemble averages or be shown
% in number of cycles.
if strcmp(DeleteFirstCycle,'Yes') == 1
    Start = 2;
    TotCyclesL = TotCyclesL-(1*length(NumCyclesL));
    TotCyclesR = TotCyclesR-(1*length(NumCyclesR));
else
    Start = 1;
end

%% Combine all kinematic data for each cycle into one matrix
% LEFT loop
CyclesL = cell(1,TotCyclesL);
j = 1;
for z = 1:length(Kinematics)
    for i = Start:NumCyclesL(z)
        CombineCycles(:,1) = Kinematics(z).Trials(1).data.left(:,i);
        CombineCycles(:,2) = Kinematics(z).Trials(2).data.left(:,i);
        CombineCycles(:,3) = Kinematics(z).Trials(3).data.left(:,i);
        CombineCycles(:,7) = Kinematics(z).Trials(4).data.left(:,i);
        CombineCycles(:,8) = Kinematics(z).Trials(5).data.left(:,i);
        CombineCycles(:,9) = Kinematics(z).Trials(6).data.left(:,i);
        CombineCycles(:,13) = Kinematics(z).Trials(7).data.left(:,i);
        CombineCycles(:,14) = Kinematics(z).Trials(8).data.left(:,i);
        CombineCycles(:,15) = Kinematics(z).Trials(9).data.left(:,i);
        CombineCycles(:,19) = Kinematics(z).Trials(10).data.left(:,i);
        CombineCycles(:,20) = Kinematics(z).Trials(11).data.left(:,i);
        CombineCycles(:,21) = Kinematics(z).Trials(12).data.left(:,i);
        CombineCycles(:,22) = zeros(101,1);
        CombineCycles(:,23) = zeros(101,1);
        CombineCycles(:,24) = zeros(101,1);
        if strcmp(FootGlobAng,'Yes') == 1   % if foot global angles desired
        CombineCycles(:,25) = Kinematics(z).Trials(13).data.left(:,i);
        CombineCycles(:,27) = Kinematics(z).Trials(14).data.left(:,i);
        end
        
        CyclesL{j} = CombineCycles;
        j=j+1;
    end
end
clearvars CombineCycles

% RIGHT loop
CyclesR = cell(1,TotCyclesR);
j = 1;
for z = 1:length(Kinematics)
    for i = Start:NumCyclesR(z)
        CombineCycles(:,4) = Kinematics(z).Trials(1).data.right(:,i);
        CombineCycles(:,5) = Kinematics(z).Trials(2).data.right(:,i);
        CombineCycles(:,6) = Kinematics(z).Trials(3).data.right(:,i);
        CombineCycles(:,10) = Kinematics(z).Trials(4).data.right(:,i);
        CombineCycles(:,11) = Kinematics(z).Trials(5).data.right(:,i);
        CombineCycles(:,12) = Kinematics(z).Trials(6).data.right(:,i);
        CombineCycles(:,16) = Kinematics(z).Trials(7).data.right(:,i);
        CombineCycles(:,17) = Kinematics(z).Trials(8).data.right(:,i);
        CombineCycles(:,18) = Kinematics(z).Trials(9).data.right(:,i);
        CombineCycles(:,22) = Kinematics(z).Trials(10).data.right(:,i);
        CombineCycles(:,23) = Kinematics(z).Trials(11).data.right(:,i);
        CombineCycles(:,24) = Kinematics(z).Trials(12).data.right(:,i);
         if strcmp(FootGlobAng,'Yes') == 1   % if foot global angles desired
        CombineCycles(:,26) = Kinematics(z).Trials(13).data.right(:,i);
        CombineCycles(:,28) = Kinematics(z).Trials(14).data.right(:,i);
        end
        
        CyclesR{j} = CombineCycles;
        j=j+1;
    end
end
clearvars CombineCycles
%% Document # of cycles present in trials
% Get rid of first cycle in number of cycles tally
if strcmp(DeleteFirstCycle,'Yes') == 1
    NumCyclesL = NumCyclesL-1;
    NumCyclesR = NumCyclesR-1;
end

% Document number of cycles in Kinematics Structure
for z = 1:NumTrials
    Kinematics(z).NumCycles = [NumCyclesL(z), NumCyclesR(z)];
end

%% Save Representative Cycle
% for i = 1:length(Kinematics)
%     Ind = strcmp(Kinematics(i).FileName,Kinematics(1).RepTrial.Name) == 1;
%     if Ind == 1
%         z = i;
%         NumRepTrial = z;
%         break
%     end
% end
NumRepTrial = Kinematics(1).RepTrial.Num(1);
L = Kinematics(1).RepTrial.RepCycles(1);
R = Kinematics(1).RepTrial.RepCycles(2);
[~,N] = size(Kinematics(NumRepTrial).Trials(1).data.left);
while L > N
    L = L-1;
end

[~,N] = size(Kinematics(NumRepTrial).Trials(1).data.right);
while R > N
    R = R-1;
end

z = NumRepTrial; 
RepData(:,1) = Kinematics(z).Trials(1).data.left(:,L);
RepData(:,2) = Kinematics(z).Trials(2).data.left(:,L);
RepData(:,3) = Kinematics(z).Trials(3).data.left(:,L);
RepData(:,4) = Kinematics(z).Trials(1).data.right(:,R);
RepData(:,5) = Kinematics(z).Trials(2).data.right(:,R);
RepData(:,6) = Kinematics(z).Trials(3).data.right(:,R);
RepData(:,7) = Kinematics(z).Trials(4).data.left(:,L);
RepData(:,8) = Kinematics(z).Trials(5).data.left(:,L);
RepData(:,9) = Kinematics(z).Trials(6).data.left(:,L);
RepData(:,10) = Kinematics(z).Trials(4).data.right(:,R);
RepData(:,11) = Kinematics(z).Trials(5).data.right(:,R);
RepData(:,12) = Kinematics(z).Trials(6).data.right(:,R);
RepData(:,13) = Kinematics(z).Trials(7).data.left(:,L);
RepData(:,14) = Kinematics(z).Trials(8).data.left(:,L);
RepData(:,15) = Kinematics(z).Trials(9).data.left(:,L);
RepData(:,16) = Kinematics(z).Trials(7).data.right(:,R);
RepData(:,17) = Kinematics(z).Trials(8).data.right(:,R);
RepData(:,18) = Kinematics(z).Trials(9).data.right(:,R);
RepData(:,19) = Kinematics(z).Trials(10).data.left(:,L);
RepData(:,20) = Kinematics(z).Trials(11).data.left(:,L);
RepData(:,21) = Kinematics(z).Trials(12).data.left(:,L);
RepData(:,22) = Kinematics(z).Trials(10).data.right(:,R);
RepData(:,23) = Kinematics(z).Trials(11).data.right(:,R);
RepData(:,24) = Kinematics(z).Trials(12).data.right(:,R);
if strcmp(FootGlobAng,'Yes') == 1   % if foot global angles desired
    RepData(:,25) = Kinematics(z).Trials(13).data.left(:,L);
    RepData(:,26) = Kinematics(z).Trials(13).data.right(:,R);
    RepData(:,27) = Kinematics(z).Trials(14).data.left(:,L);
    RepData(:,28) = Kinematics(z).Trials(14).data.right(:,R);
end

Kinematics(1).RepCycle = RepData;

%% Calculate EA and STD
% Concatenate all the trials for each side
Kine_L = cat(3,CyclesL{:}); % Arrange the kinematic data into a 3D array for averaging
Kine_R = cat(3,CyclesR{:});

% Compute the Ensemble Average for each side
EA_L = mean(Kine_L,3,'omitnan'); % average the kinematic data, ignores all NaN values
EA_R = mean(Kine_R,3,'omitnan');
% add zeros to end of left side to make matricies similar dimensions
if strcmp(FootGlobAng,'Yes') == 1   % if foot global angles desired
    EA_L(:,28) = zeros(101,1);
end
% Combine Sides
Kinematics(1).EnsembleAverage = EA_L+EA_R;

% Compute the Standard Deviation for each side
STD_L = std(Kine_L,0,3);
STD_R = std(Kine_R,0,3);
% add zeros to end of left side to make matricies similar dimensions
if strcmp(FootGlobAng,'Yes') == 1   % if foot global angles desired
    STD_L(:,28) = zeros(101,1);
end
% Combine Sides
Kinematics(1).STD = STD_L+STD_R;

%% Compute GDI

if strcmp(GDI_RCorEA,'RC') == 1
    % GDI is computed in the function DetermineRepTrialCyc
    % UNCOMMENT TO CHOOSE HEMIPLEGIA SIDE EVERY TIME IN LOOPS
%     prompt = ('Is the subject hemiplegic?');
%     Hemi = questdlg(prompt, 'Hemi','Left','Right','No','No');
%     GDI = GDI_HemiInitial(char(Kinematics(1).RepTrial.Name), Hemi);

elseif strcmp(GDI_RCorEA,'EA') == 1
    GDI = GDI_Calculator_Fast(0,0,0,Kinematics(1).EnsembleAverage);
end
Kinematics(1).GDI = GDI;

end





