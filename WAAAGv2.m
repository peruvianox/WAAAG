%% WAAAG - Walking Ability At A Glance

% WAAAG will take a set of .c3d trials (Current Condition) either as full
% kinematic data or foot markers only (temporal spatial measures) and
% automatically analyze the data using the GetKinematics code. Using this
% initial information, WAAAG will make graphical comparisons to a target
% reference (normals/controls) and any additional comparisons (AFO trials,
% last gait analysis, Extra). This function requires at least one processed
% C3D file to run. All other inputs are C3Ds (either full 3D kinematics, or foot 
% markers only). The outputs are displayed as bars of the average metrics and 
% representative cycles shown as dots. 

% INPUTS
% one or multiple .C3D files for each condiiton that contain full 3D
% kinematics or foot marker trajectories. 

% OUTPUTS
% GaitMeasures.png figure : 
% a norm-referenced set of  plots of temporal spatial parameters for all desired walking conditions. 
% GDI-MAP.png figure :
% a set of GDI average and rep cycle plots for each condition. Also a movement analysis 
% profile of the subjects kinematic deviations from the norm for each condition. 
% WAAAG.xlsx :
% a series of spreadsheets. The first (WAAAG) contains numerical data for each of the
% conditions and an overview of the included walking trials.  Then a breakdown of each 
% condition (COND_Trial_Cyc) follows, with kinematic residual (difference between each 
% gait cycle and the condition average) rankings, and then a ranking of the cycles for the top 3 ranked trials.  
% .GCD files : 
% containing the average rotations of the kinematics for all trials and conditions. 
% WAAG_data.mat :
% containing all the measurements and inputs in a structured format. This
% can be used to re-run the analysis at a later date. 


% Version 1.1
% 1/17/2018
% Updates in this version include:
% 1) reversing the bars and dots so that the average is shown by the bars
%     and the rep trial is shown by a dot
% 2) reorgainzing the plots so that stride length is also included in main
%    temporal spatial plots GDI plot will now be included in a separate graph 
%       that also contains a movement analysis profile (MAP)
%       see Baker et al. 2009. Gait & Posture. 30 (3): 265-269 for more info
%       doi:10.1016/j.gaitpost.2009.05.020.
% 3) added trial and cycle count metrics to the excel spreadsheet
%     (next to Average Rotations)
% 4) updated image file type to PNG for simpler pasting into EPIC
% 5) changed to automatica excel copying and saving - no need to copy the
%     template and rename
% 6) added lists of the representative trials for each condition to the
%     main excel page
% 7) version # and date processed on excel spreadsheet

% Ricky Pimentel
% Center for Gait and Movement Analysis, Children's Hospital Colorado

clear;
clc;
warning off;
addpath(genpath('Called Functions'));

%% Load Data
prompt = 'What is the age of the subject?';
age = inputdlg(prompt);

% Rename current condition
prompt = 'What would you like to name the current condition?';
Name.Curr = inputdlg(prompt);

Type.Current = questdlg('What type of files are being imported?','Input Type', 'Full Kinematics','Temp-Spat Only','Full Kinematics');
if strcmp(Type.Current, 'Full Kinematics') == 1 % if full kinematics C3Ds
    KineData = GetKinematics;
    
    GDIRep.Current = GDI_Calculator_Fast(0,0,0,KineData(1).RepCycle); % compute RC GDI
    % compute GDIs for every cycle within every trial
    CycleCount = 1;
    for i = 1:length(KineData)% loop through all trials
        MinCycles = min(KineData(i).NumCycles); % find minimum # of cycles on each side for each trial
        for j = 1:MinCycles % loop through # of cycles
            KineMatrix(:,:,CycleCount) = [KineData(i).Trials(1).data.left(:,j), KineData(i).Trials(2).data.left(:,j), KineData(i).Trials(3).data.left(:,j),... % Left Pelvis
                KineData(i).Trials(1).data.right(:,j), KineData(i).Trials(2).data.right(:,j), KineData(i).Trials(3).data.right(:,j),... % Right Pelvis
                KineData(i).Trials(4).data.left(:,j), KineData(i).Trials(5).data.left(:,j), KineData(i).Trials(6).data.left(:,j),... % Left Hip
                KineData(i).Trials(4).data.right(:,j), KineData(i).Trials(5).data.right(:,j), KineData(i).Trials(6).data.right(:,j),... % Right Hip
                KineData(i).Trials(7).data.left(:,j), KineData(i).Trials(8).data.left(:,j), KineData(i).Trials(9).data.left(:,j),... % Left Knee
                KineData(i).Trials(7).data.right(:,j), KineData(i).Trials(8).data.right(:,j), KineData(i).Trials(9).data.right(:,j),... % Right Knee
                KineData(i).Trials(10).data.left(:,j), KineData(i).Trials(11).data.left(:,j), KineData(i).Trials(12).data.left(:,j),... % Left Ankle
                KineData(i).Trials(10).data.right(:,j), KineData(i).Trials(11).data.right(:,j), KineData(i).Trials(12).data.right(:,j)]; % Right Ankle
            CycleCount = CycleCount + 1; % advance to save next cycle in next spot
        end
    end
    [~,~,NumGDITrials] = size(KineMatrix); % find total # of cycles saved
    for i = 1:NumGDITrials
        TrialsGDI.Current(i,:) = GDI_Calculator_Fast(0,0,0,KineMatrix(:,:,i)); % compute GDI
    end
    KineData(1).TrialsGDI = TrialsGDI.Current;  % save GDI data into kinematics structure
    KineData(1).KineMatrix = KineMatrix;  % save cycle kinematics data into kinematics structure
    
else % if temp spat C3Ds
    [files2get, ~] = uigetfile('*.c3d*', 'Please select all C3D file(s) to be analyzed', 'MultiSelect', 'on');
    KineData = GAMSmini_BatchGaitMeasures(files2get,str2double(char(age)));
    GDIRep.Current = [0 0];
    KineData(1).TrialsGDI = [0 0 0];
    for i = 1:length(KineData)
        KineData(i).GDI = [0 0 0];
        if length(KineData) > 1
            KineData(i).FileName = files2get{i};
        else
            KineData(i).FileName = files2get;
        end
    end
end

NumConds = 1;

% Load Control Data
TSN = GetTempSpatNorms(age);
ControlKine = GetControlKinematics(age);
% TempSpat norms are organized as:
% [Cadence, Stride Length, Walking Speed, LStepLength,RStepLength,StepWidth]
% steps/min   m                   m/min                 cm              cm              cm

clearvars GDIplus GDIminus DiffPlus DiffMinus GDIavg Diff

%% Define mean variables to be plotted
for i = 1:length(KineData)
    % Overall Gait Performance
    Cadence.Current.Measures(i) = KineData(i).TempSpat(1).data;
    GaitSpeed.Current.Measures(i) = KineData(i).TempSpat(2).data;
    % GDI already exists in KineData(1).GDI
    % Spatial Data
    StrideLength.Current.Measures(i) = KineData(i).TempSpat(3).data;
    StepLength.Current.Left.Measures(i) = KineData(i).TempSpat(6).data;
    StepLength.Current.Right.Measures(i) = KineData(i).TempSpat(15).data;
    StepWidth.Current.AvgMeasures(i) = KineData(i).TempSpat(22).data;
    % Temporal Data
    Stance.Current.Left.Measures(i) = KineData(i).TempSpat(9).data;
    Stance.Current.Right.Measures(i) = KineData(i).TempSpat(18).data;
    Swing.Current.Left.Measures(i) = 100 - Stance.Current.Left.Measures(i);
    Swing.Current.Right.Measures(i) = 100 - Stance.Current.Right.Measures(i);
    DS1.Current.Left.Measures(i) = KineData(i).TempSpat(10).data;
    DS1.Current.Right.Measures(i) = KineData(i).TempSpat(19).data;
    SLS.Current.Left.Measures(i) = KineData(i).TempSpat(11).data;
    SLS.Current.Right.Measures(i) = KineData(i).TempSpat(20).data;
    DS2.Current.Left.Measures(i) = KineData(i).TempSpat(12).data;
    DS2.Current.Right.Measures(i) = KineData(i).TempSpat(21).data;
    % STDs
    StrideLength.Current.StdMeasures(i) = KineData(i).TempSpat(3).data;
    StepLength.Current.Left.StdMeasures(i) = KineData(i).TempSpat(24).data;
    StepLength.Current.Right.StdMeasures(i) = KineData(i).TempSpat(29).data;
    StepWidth.Current.StdMeasures(i) = KineData(i).TempSpat(23).data;
    % Temporal Data
    Stance.Current.Left.StdMeasures(i) = KineData(i).TempSpat(25).data;
    Stance.Current.Right.StdMeasures(i) = KineData(i).TempSpat(30).data;
    Swing.Current.Left.StdMeasures(i) = Stance.Current.Left.StdMeasures(i);
    Swing.Current.Right.StdMeasures(i) = Stance.Current.Right.StdMeasures(i);
    DS1.Current.Left.StdMeasures(i) = KineData(i).TempSpat(26).data;
    DS1.Current.Right.StdMeasures(i) = KineData(i).TempSpat(31).data;
    SLS.Current.Left.StdMeasures(i) = KineData(i).TempSpat(27).data;
    SLS.Current.Right.StdMeasures(i) = KineData(i).TempSpat(32).data;
    DS2.Current.Left.StdMeasures(i) = KineData(i).TempSpat(28).data;
    DS2.Current.Right.StdMeasures(i) = KineData(i).TempSpat(33).data;
end

% Calculate averages and STDs
% Overall Gait Performance
Cadence.Current.Avg = mean(Cadence.Current.Measures);
Cadence.Current.STD = std(Cadence.Current.Measures);
GaitSpeed.Current.Avg = mean(GaitSpeed.Current.Measures);
GaitSpeed.Current.STD = std(GaitSpeed.Current.Measures);
% Spatial Data
StrideLength.Current.Avg = mean(StrideLength.Current.Measures);
StrideLength.Current.STD = std(StrideLength.Current.Measures);
StepLength.Current.Left.Avg = mean(StepLength.Current.Left.Measures);
StepLength.Current.Left.STD = mean(StepLength.Current.Left.StdMeasures);
StepLength.Current.Right.Avg = mean(StepLength.Current.Right.Measures);
StepLength.Current.Right.STD = mean(StepLength.Current.Right.StdMeasures);
StepWidth.Current.Avg = mean(StepWidth.Current.AvgMeasures);
StepWidth.Current.STD = mean(StepWidth.Current.StdMeasures);
% Temporal Data
Stance.Current.Left.Avg = mean(Stance.Current.Left.Measures);
Stance.Current.Left.STD = mean(Stance.Current.Left.StdMeasures);
Stance.Current.Right.Avg = mean(Stance.Current.Right.Measures);
Stance.Current.Right.STD = mean(Stance.Current.Right.StdMeasures);
Swing.Current.Left.Avg = mean(Swing.Current.Left.Measures);
Swing.Current.Left.STD = mean(Swing.Current.Left.StdMeasures);
Swing.Current.Right.Avg = mean(Swing.Current.Right.Measures);
Swing.Current.Right.STD = mean(Swing.Current.Right.StdMeasures);

DS1.Current.Left.Avg = mean(DS1.Current.Left.Measures);
DS1.Current.Left.STD = mean(DS1.Current.Left.StdMeasures);
DS1.Current.Right.Avg = mean(DS1.Current.Right.Measures);
DS1.Current.Right.STD = mean(DS1.Current.Right.StdMeasures);
SLS.Current.Left.Avg = mean(SLS.Current.Left.Measures);
SLS.Current.Left.STD = mean(SLS.Current.Left.StdMeasures);
SLS.Current.Right.Avg = mean(SLS.Current.Right.Measures);
SLS.Current.Right.STD = mean(SLS.Current.Right.StdMeasures);
DS2.Current.Left.Avg = mean(DS2.Current.Left.Measures);
DS2.Current.Left.STD = mean(DS2.Current.Left.StdMeasures);
DS2.Current.Right.Avg = mean(DS2.Current.Right.Measures);
DS2.Current.Right.STD = mean(DS2.Current.Right.StdMeasures);

%% define rep trial variables
if strcmp(Type.Current, 'Full Kinematics') == 1
    [~, RepTrial] = min(KineData(1).RepTrial.Num);
    % Calculate averages and STDs
    % Overall Gait Performance
    Cadence.Current.Rep = KineData(RepTrial).TempSpat(1).data;
    GaitSpeed.Current.Rep = KineData(RepTrial).TempSpat(2).data;
    % Spatial Data
    StrideLength.Current.Rep = KineData(RepTrial).TempSpat(3).data;
    StepLength.Current.Left.Rep = KineData(RepTrial).TempSpat(6).data;
    StepLength.Current.Right.Rep = KineData(RepTrial).TempSpat(15).data;
    StepWidth.Current.Rep = KineData(RepTrial).TempSpat(22).data;
    % Temporal Data
    Stance.Current.Left.Rep = KineData(RepTrial).TempSpat(9).data;
    Stance.Current.Right.Rep = KineData(RepTrial).TempSpat(18).data;
    Swing.Current.Left.Rep = 100 - KineData(RepTrial).TempSpat(9).data;
    Swing.Current.Right.Rep = 100 - KineData(RepTrial).TempSpat(18).data;
    
    DS1.Current.Left.Rep = KineData(RepTrial).TempSpat(10).data;
    DS1.Current.Right.Rep = KineData(RepTrial).TempSpat(19).data;
    SLS.Current.Left.Rep = KineData(RepTrial).TempSpat(11).data;
    SLS.Current.Right.Rep = KineData(RepTrial).TempSpat(20).data;
    DS2.Current.Left.Rep = KineData(RepTrial).TempSpat(12).data;
    DS2.Current.Right.Rep = KineData(RepTrial).TempSpat(21).data;
end

%% Add AFO condition if desired
prompt = 'Would you like to add an AFO condition?';
AddAFO = questdlg(prompt,'AFO Condition','No','Yes: Full','Yes: Temp-Spat','No');
if strcmp(AddAFO, 'No') == 0
    % save type of AFO condition
    if strcmp(AddAFO,'Yes: Full') == 1
        Type.AFO = 'Full Kinematics';
    else
        Type.AFO = 'Temp-Spat Only';
    end
    AddAFO = 'Yes';
    NumConds = NumConds +1;
    % Rename AFO condition
    prompt = 'What would you like to name the AFO condition?';
    Name.AFO = inputdlg(prompt);
    if strcmp(Type.AFO,'Full Kinematics') == 1
        % Get AFO data
        AFOData = GetKinematics; % load AFO c3d data
        % GDI Computations
        GDIRep.AFO = GDI_Calculator_Fast(0,0,0,AFOData(1).RepCycle); % compute RC GDI
        % compute GDIs for every cycle within every trial
        CycleCount = 1;
        for i = 1:length(AFOData)% loop through all trials
            MinCycles = min(AFOData(i).NumCycles); % find minimum # of cycles on each side for each trial
            for j = 1:MinCycles % loop through # of cycles
                AFOMatrix(:,:,CycleCount) = [AFOData(i).Trials(1).data.left(:,j), AFOData(i).Trials(2).data.left(:,j), AFOData(i).Trials(3).data.left(:,j),... % Left Pelvis
                    AFOData(i).Trials(1).data.right(:,j), AFOData(i).Trials(2).data.right(:,j), AFOData(i).Trials(3).data.right(:,j),... % Right Pelvis
                    AFOData(i).Trials(4).data.left(:,j), AFOData(i).Trials(5).data.left(:,j), AFOData(i).Trials(6).data.left(:,j),... % Left Hip
                    AFOData(i).Trials(4).data.right(:,j), AFOData(i).Trials(5).data.right(:,j), AFOData(i).Trials(6).data.right(:,j),... % Right Hip
                    AFOData(i).Trials(7).data.left(:,j), AFOData(i).Trials(8).data.left(:,j), AFOData(i).Trials(9).data.left(:,j),... % Left Knee
                    AFOData(i).Trials(7).data.right(:,j), AFOData(i).Trials(8).data.right(:,j), AFOData(i).Trials(9).data.right(:,j),... % Right Knee
                    AFOData(i).Trials(10).data.left(:,j), AFOData(i).Trials(11).data.left(:,j), AFOData(i).Trials(12).data.left(:,j),... % Left Ankle
                    AFOData(i).Trials(10).data.right(:,j), AFOData(i).Trials(11).data.right(:,j), AFOData(i).Trials(12).data.right(:,j)]; % Right Ankle
                CycleCount = CycleCount + 1; % advance to save next cycle in next spot
            end
        end
        [~,~,NumGDITrials] = size(AFOMatrix); % find total # of cycles saved
        for i = 1:NumGDITrials
            TrialsGDI.AFO(i,:) = GDI_Calculator_Fast(0,0,0,AFOMatrix(:,:,i)); % compute GDI
        end
        AFOData(1).TrialsGDI = TrialsGDI.AFO;  % save GDI data into kinematics structure
        AFOData(1).KineMatrix = AFOMatrix;  % save cycle kinematics data into kinematics structure
        
    else
        [files2get, ~] = uigetfile('*.c3d*', 'Please select all C3D file(s) to be analyzed', 'MultiSelect', 'on');
        AFOData = GAMSmini_BatchGaitMeasures(files2get,str2double(char(age)));
        GDIRep.AFO = [0 0];
        AFOData(1).TrialsGDI = [0 0 0];
        for i = 1:length(AFOData)
            AFOData(i).GDI = [0 0 0];
            if length(AFOData) > 1
                AFOData(i).FileName = files2get{i};
            else
                AFOData(i).FileName = files2get;
            end
        end
    end
    
    % Define AFO variables to be plotted
    for i = 1:length(AFOData)
        % Overall Gait Performance
        Cadence.AFO.Measures(i) = AFOData(i).TempSpat(1).data;
        GaitSpeed.AFO.Measures(i) = AFOData(i).TempSpat(2).data;
        % GDI already exists in AFOData(1).GDI
        % Spatial Data
        StrideLength.AFO.Measures(i) = AFOData(i).TempSpat(3).data;
        StepLength.AFO.Left.Measures(i) = AFOData(i).TempSpat(6).data;
        StepLength.AFO.Right.Measures(i) = AFOData(i).TempSpat(15).data;
        StepWidth.AFO.AvgMeasures(i) = AFOData(i).TempSpat(22).data;
        % Temporal Data
        Stance.AFO.Left.Measures(i) = AFOData(i).TempSpat(9).data;
        Stance.AFO.Right.Measures(i) = AFOData(i).TempSpat(18).data;
        Swing.AFO.Left.Measures(i) = 100 - Stance.AFO.Left.Measures(i);
        Swing.AFO.Right.Measures(i) = 100 - Stance.AFO.Right.Measures(i);
        DS1.AFO.Left.Measures(i) = AFOData(i).TempSpat(10).data;
        DS1.AFO.Right.Measures(i) = AFOData(i).TempSpat(19).data;
        SLS.AFO.Left.Measures(i) = AFOData(i).TempSpat(11).data;
        SLS.AFO.Right.Measures(i) = AFOData(i).TempSpat(20).data;
        DS2.AFO.Left.Measures(i) = AFOData(i).TempSpat(12).data;
        DS2.AFO.Right.Measures(i) = AFOData(i).TempSpat(21).data;
        % STDs
        StrideLength.AFO.StdMeasures(i) = AFOData(i).TempSpat(3).data;
        StepLength.AFO.Left.StdMeasures(i) = AFOData(i).TempSpat(24).data;
        StepLength.AFO.Right.StdMeasures(i) = AFOData(i).TempSpat(29).data;
        StepWidth.AFO.StdMeasures(i) = AFOData(i).TempSpat(23).data;
        % Temporal Data
        Stance.AFO.Left.StdMeasures(i) = AFOData(i).TempSpat(25).data;
        Stance.AFO.Right.StdMeasures(i) = AFOData(i).TempSpat(30).data;
        Swing.AFO.Left.StdMeasures(i) = Stance.AFO.Left.StdMeasures(i);
        Swing.AFO.Right.StdMeasures(i) = Stance.AFO.Right.StdMeasures(i);
        DS1.AFO.Left.StdMeasures(i) = AFOData(i).TempSpat(26).data;
        DS1.AFO.Right.StdMeasures(i) = AFOData(i).TempSpat(31).data;
        SLS.AFO.Left.StdMeasures(i) = AFOData(i).TempSpat(27).data;
        SLS.AFO.Right.StdMeasures(i) = AFOData(i).TempSpat(32).data;
        DS2.AFO.Left.StdMeasures(i) = AFOData(i).TempSpat(28).data;
        DS2.AFO.Right.StdMeasures(i) = AFOData(i).TempSpat(33).data;
    end
    
    % Calculate averages and STDs
    % Overall Gait Performance
    Cadence.AFO.Avg = mean(Cadence.AFO.Measures);
    Cadence.AFO.STD = std(Cadence.AFO.Measures);
    GaitSpeed.AFO.Avg = mean(GaitSpeed.AFO.Measures);
    GaitSpeed.AFO.STD = std(GaitSpeed.AFO.Measures);
    % Spatial Data
    StrideLength.AFO.Avg = mean(StrideLength.AFO.Measures);
    StrideLength.AFO.STD = std(StrideLength.AFO.Measures);
    StepLength.AFO.Left.Avg = mean(StepLength.AFO.Left.Measures);
    StepLength.AFO.Left.STD = mean(StepLength.AFO.Left.StdMeasures);
    StepLength.AFO.Right.Avg = mean(StepLength.AFO.Right.Measures);
    StepLength.AFO.Right.STD = mean(StepLength.AFO.Right.StdMeasures);
    StepWidth.AFO.Avg = mean(StepWidth.AFO.AvgMeasures);
    StepWidth.AFO.STD = mean(StepWidth.AFO.StdMeasures);
    % Temporal Data
    Stance.AFO.Left.Avg = mean(Stance.AFO.Left.Measures);
    Stance.AFO.Left.STD = mean(Stance.AFO.Left.StdMeasures);
    Stance.AFO.Right.Avg = mean(Stance.AFO.Right.Measures);
    Stance.AFO.Right.STD = mean(Stance.AFO.Right.StdMeasures);
    Swing.AFO.Left.Avg = mean(Swing.AFO.Left.Measures);
    Swing.AFO.Left.STD = mean(Swing.AFO.Left.StdMeasures);
    Swing.AFO.Right.Avg = mean(Swing.AFO.Right.Measures);
    Swing.AFO.Right.STD = mean(Swing.AFO.Right.StdMeasures);
    
    DS1.AFO.Left.Avg = mean(DS1.AFO.Left.Measures);
    DS1.AFO.Left.STD = mean(DS1.AFO.Left.StdMeasures);
    DS1.AFO.Right.Avg = mean(DS1.AFO.Right.Measures);
    DS1.AFO.Right.STD = mean(DS1.AFO.Right.StdMeasures);
    SLS.AFO.Left.Avg = mean(SLS.AFO.Left.Measures);
    SLS.AFO.Left.STD = mean(SLS.AFO.Left.StdMeasures);
    SLS.AFO.Right.Avg = mean(SLS.AFO.Right.Measures);
    SLS.AFO.Right.STD = mean(SLS.AFO.Right.StdMeasures);
    DS2.AFO.Left.Avg = mean(DS2.AFO.Left.Measures);
    DS2.AFO.Left.STD = mean(DS2.AFO.Left.StdMeasures);
    DS2.AFO.Right.Avg = mean(DS2.AFO.Right.Measures);
    DS2.AFO.Right.STD = mean(DS2.AFO.Right.StdMeasures);
    
    if strcmp(Type.AFO, 'Full Kinematics') == 1
        % define rep trial variables
        RepTrial = AFOData(1).RepTrial.Num(1);
        % Calculate averages and STDs
        % Overall Gait Performance
        Cadence.AFO.Rep = AFOData(RepTrial).TempSpat(1).data;
        GaitSpeed.AFO.Rep = AFOData(RepTrial).TempSpat(2).data;
        % Spatial Data
        StrideLength.AFO.Rep = AFOData(RepTrial).TempSpat(3).data;
        StepLength.AFO.Left.Rep = AFOData(RepTrial).TempSpat(6).data;
        StepLength.AFO.Right.Rep = AFOData(RepTrial).TempSpat(15).data;
        StepWidth.AFO.Rep = AFOData(RepTrial).TempSpat(22).data;
        % Temporal Data
        Stance.AFO.Left.Rep = AFOData(RepTrial).TempSpat(9).data;
        Stance.AFO.Right.Rep = AFOData(RepTrial).TempSpat(18).data;
        Swing.AFO.Left.Rep = 100 - AFOData(RepTrial).TempSpat(9).data;
        Swing.AFO.Right.Rep = 100 - AFOData(RepTrial).TempSpat(18).data;
        
        DS1.AFO.Left.Rep = AFOData(RepTrial).TempSpat(10).data;
        DS1.AFO.Right.Rep = AFOData(RepTrial).TempSpat(19).data;
        SLS.AFO.Left.Rep = AFOData(RepTrial).TempSpat(11).data;
        SLS.AFO.Right.Rep = AFOData(RepTrial).TempSpat(20).data;
        DS2.AFO.Left.Rep = AFOData(RepTrial).TempSpat(12).data;
        DS2.AFO.Right.Rep = AFOData(RepTrial).TempSpat(21).data;
    end
end

clearvars GDIplus GDIminus DiffPlus DiffMinus GDIavg Diff

%% Add Last Condition if desired
prompt = 'Would you like to add a Last condition?';
AddLast = questdlg(prompt,'Last Condition','No','Yes: Full','Yes: Temp-Spat','No');
if strcmp(AddLast, 'No') == 0
    % save type of Last condition
    if strcmp(AddLast,'Yes: Full') == 1
        Type.Last = 'Full Kinematics';
    else
        Type.Last = 'Temp-Spat Only';
    end
    AddLast = 'Yes';
    NumConds = NumConds +1;
    % Rename AFO condition
    prompt = 'What would you like to name the Last condition?';
    Name.Last = inputdlg(prompt);
    
    % Get data for Last condition
    if strcmp(Type.Last,'Full Kinematics') == 1 % if full Kinematics
        % Get Last data
        LastData = GetKinematics;
        % GDI Computations
        GDIRep.Last = GDI_Calculator_Fast(0,0,0,LastData(1).RepCycle); % compute RC GDI
        % compute GDIs for every cycle within every trial
        CycleCount = 1;
        for i = 1:length(LastData)% loop through all trials
            MinCycles = min(LastData(i).NumCycles); % find minimum # of cycles on each side for each trial
            for j = 1:MinCycles % loop through # of cycles
                LastMatrix(:,:,CycleCount) = [LastData(i).Trials(1).data.left(:,j), LastData(i).Trials(2).data.left(:,j), LastData(i).Trials(3).data.left(:,j),... % Left Pelvis
                    LastData(i).Trials(1).data.right(:,j), LastData(i).Trials(2).data.right(:,j), LastData(i).Trials(3).data.right(:,j),... % Right Pelvis
                    LastData(i).Trials(4).data.left(:,j), LastData(i).Trials(5).data.left(:,j), LastData(i).Trials(6).data.left(:,j),... % Left Hip
                    LastData(i).Trials(4).data.right(:,j), LastData(i).Trials(5).data.right(:,j), LastData(i).Trials(6).data.right(:,j),... % Right Hip
                    LastData(i).Trials(7).data.left(:,j), LastData(i).Trials(8).data.left(:,j), LastData(i).Trials(9).data.left(:,j),... % Left Knee
                    LastData(i).Trials(7).data.right(:,j), LastData(i).Trials(8).data.right(:,j), LastData(i).Trials(9).data.right(:,j),... % Right Knee
                    LastData(i).Trials(10).data.left(:,j), LastData(i).Trials(11).data.left(:,j), LastData(i).Trials(12).data.left(:,j),... % Left Ankle
                    LastData(i).Trials(10).data.right(:,j), LastData(i).Trials(11).data.right(:,j), LastData(i).Trials(12).data.right(:,j)]; % Right Ankle
                CycleCount = CycleCount + 1; % advance to save next cycle in next spot
            end
        end
        [~,~,NumGDITrials] = size(LastMatrix); % find total # of cycles saved
        for i = 1:NumGDITrials
            TrialsGDI.Last(i,:) = GDI_Calculator_Fast(0,0,0,LastMatrix(:,:,i)); % compute GDI
        end
        LastData(1).TrialsGDI = TrialsGDI.Last;  % save GDI data into kinematics structure
        LastData(1).KineMatrix = LastMatrix;  % save cycle kinematics data into kinematics structure
    else
        [files2get, ~] = uigetfile('*.c3d*', 'Please select all C3D file(s) to be analyzed', 'MultiSelect', 'on');
        LastData = GAMSmini_BatchGaitMeasures(files2get,str2double(char(age)));
        GDIRep.Last = [0 0];
        LastData(1).TrialsGDI = [0 0 0];
        for i = 1:length(LastData)
            LastData(i).GDI = [0 0 0];
            if length(LastData) > 1
                LastData(i).FileName = files2get{i};
            else
                LastData(i).FileName = files2get;
            end
        end
    end
    
    
    % Define Last variables to be plotted
    for i = 1:length(LastData)
        % Overall Gait Performance
        Cadence.Last.Measures(i) = LastData(i).TempSpat(1).data;
        GaitSpeed.Last.Measures(i) = LastData(i).TempSpat(2).data;
        % GDI already exists in LastData(1).GDI
        % Spatial Data
        StrideLength.Last.Measures(i) = LastData(i).TempSpat(3).data;
        StepLength.Last.Left.Measures(i) = LastData(i).TempSpat(6).data;
        StepLength.Last.Right.Measures(i) = LastData(i).TempSpat(15).data;
        StepWidth.Last.AvgMeasures(i) = LastData(i).TempSpat(22).data;
        % Temporal Data
        Stance.Last.Left.Measures(i) = LastData(i).TempSpat(9).data;
        Stance.Last.Right.Measures(i) = LastData(i).TempSpat(18).data;
        Swing.Last.Left.Measures(i) = 100 - Stance.Last.Left.Measures(i);
        Swing.Last.Right.Measures(i) = 100 - Stance.Last.Right.Measures(i);
        DS1.Last.Left.Measures(i) = LastData(i).TempSpat(10).data;
        DS1.Last.Right.Measures(i) = LastData(i).TempSpat(19).data;
        SLS.Last.Left.Measures(i) = LastData(i).TempSpat(11).data;
        SLS.Last.Right.Measures(i) = LastData(i).TempSpat(20).data;
        DS2.Last.Left.Measures(i) = LastData(i).TempSpat(12).data;
        DS2.Last.Right.Measures(i) = LastData(i).TempSpat(21).data;
        % STDs
        StrideLength.Last.StdMeasures(i) = LastData(i).TempSpat(3).data;
        StepLength.Last.Left.StdMeasures(i) = LastData(i).TempSpat(24).data;
        StepLength.Last.Right.StdMeasures(i) = LastData(i).TempSpat(29).data;
        StepWidth.Last.StdMeasures(i) = LastData(i).TempSpat(23).data;
        % Temporal Data
        Stance.Last.Left.StdMeasures(i) = LastData(i).TempSpat(25).data;
        Stance.Last.Right.StdMeasures(i) = LastData(i).TempSpat(30).data;
        Swing.Last.Left.StdMeasures(i) = Stance.Last.Left.StdMeasures(i);
        Swing.Last.Right.StdMeasures(i) = Stance.Last.Right.StdMeasures(i);
        DS1.Last.Left.StdMeasures(i) = LastData(i).TempSpat(26).data;
        DS1.Last.Right.StdMeasures(i) = LastData(i).TempSpat(31).data;
        SLS.Last.Left.StdMeasures(i) = LastData(i).TempSpat(27).data;
        SLS.Last.Right.StdMeasures(i) = LastData(i).TempSpat(32).data;
        DS2.Last.Left.StdMeasures(i) = LastData(i).TempSpat(28).data;
        DS2.Last.Right.StdMeasures(i) = LastData(i).TempSpat(33).data;
    end
    
    % Calculate averages and STDs
    % Overall Gait Performance
    Cadence.Last.Avg = mean(Cadence.Last.Measures);
    Cadence.Last.STD = std(Cadence.Last.Measures);
    GaitSpeed.Last.Avg = mean(GaitSpeed.Last.Measures);
    GaitSpeed.Last.STD = std(GaitSpeed.Last.Measures);
    % Spatial Data
    StrideLength.Last.Avg = mean(StrideLength.Last.Measures);
    StrideLength.Last.STD = std(StrideLength.Last.Measures);
    StepLength.Last.Left.Avg = mean(StepLength.Last.Left.Measures);
    StepLength.Last.Left.STD = mean(StepLength.Last.Left.StdMeasures);
    StepLength.Last.Right.Avg = mean(StepLength.Last.Right.Measures);
    StepLength.Last.Right.STD = mean(StepLength.Last.Right.StdMeasures);
    StepWidth.Last.Avg = mean(StepWidth.Last.AvgMeasures);
    StepWidth.Last.STD = mean(StepWidth.Last.StdMeasures);
    % Temporal Data
    Stance.Last.Left.Avg = mean(Stance.Last.Left.Measures);
    Stance.Last.Left.STD = mean(Stance.Last.Left.StdMeasures);
    Stance.Last.Right.Avg = mean(Stance.Last.Right.Measures);
    Stance.Last.Right.STD = mean(Stance.Last.Right.StdMeasures);
    Swing.Last.Left.Avg = mean(Swing.Last.Left.Measures);
    Swing.Last.Left.STD = mean(Swing.Last.Left.StdMeasures);
    Swing.Last.Right.Avg = mean(Swing.Last.Right.Measures);
    Swing.Last.Right.STD = mean(Swing.Last.Right.StdMeasures);
    
    DS1.Last.Left.Avg = mean(DS1.Last.Left.Measures);
    DS1.Last.Left.STD = mean(DS1.Last.Left.StdMeasures);
    DS1.Last.Right.Avg = mean(DS1.Last.Right.Measures);
    DS1.Last.Right.STD = mean(DS1.Last.Right.StdMeasures);
    SLS.Last.Left.Avg = mean(SLS.Last.Left.Measures);
    SLS.Last.Left.STD = mean(SLS.Last.Left.StdMeasures);
    SLS.Last.Right.Avg = mean(SLS.Last.Right.Measures);
    SLS.Last.Right.STD = mean(SLS.Last.Right.StdMeasures);
    DS2.Last.Left.Avg = mean(DS2.Last.Left.Measures);
    DS2.Last.Left.STD = mean(DS2.Last.Left.StdMeasures);
    DS2.Last.Right.Avg = mean(DS2.Last.Right.Measures);
    DS2.Last.Right.STD = mean(DS2.Last.Right.StdMeasures);
    
    if strcmp(Type.Last, 'Full Kinematics') == 1
        % define rep trial variables
        RepTrial = LastData(1).RepTrial.Num(1);
        % Calculate averages and STDs
        % Overall Gait Performance
        Cadence.Last.Rep = LastData(RepTrial).TempSpat(1).data;
        GaitSpeed.Last.Rep = LastData(RepTrial).TempSpat(2).data;
        % Spatial Data
        StrideLength.Last.Rep = LastData(RepTrial).TempSpat(3).data;
        StepLength.Last.Left.Rep = LastData(RepTrial).TempSpat(6).data;
        StepLength.Last.Right.Rep = LastData(RepTrial).TempSpat(15).data;
        StepWidth.Last.Rep = LastData(RepTrial).TempSpat(22).data;
        % Temporal Data
        Stance.Last.Left.Rep = LastData(RepTrial).TempSpat(9).data;
        Stance.Last.Right.Rep = LastData(RepTrial).TempSpat(18).data;
        Swing.Last.Left.Rep = 100 - LastData(RepTrial).TempSpat(9).data;
        Swing.Last.Right.Rep = 100 - LastData(RepTrial).TempSpat(18).data;
        
        DS1.Last.Left.Rep = LastData(RepTrial).TempSpat(10).data;
        DS1.Last.Right.Rep = LastData(RepTrial).TempSpat(19).data;
        SLS.Last.Left.Rep = LastData(RepTrial).TempSpat(11).data;
        SLS.Last.Right.Rep = LastData(RepTrial).TempSpat(20).data;
        DS2.Last.Left.Rep = LastData(RepTrial).TempSpat(12).data;
        DS2.Last.Right.Rep = LastData(RepTrial).TempSpat(21).data;
    end
end
clearvars GDIplus GDIminus DiffPlus DiffMinus GDIavg Diff

%% Add Extra Condiditon if Desired
prompt = 'Would you like to add another condition?';
AddExtra = questdlg(prompt,'Extra Condition','No','Yes: Full','Yes: Temp-Spat','No');
if strcmp(AddExtra, 'No') == 0
    % save type of Last condition
    if strcmp(AddExtra,'Yes: Full') == 1
        Type.Extra = 'Full Kinematics';
    else
        Type.Extra = 'Temp-Spat Only';
    end
    AddExtra = 'Yes';
    NumConds = NumConds +1;
    % Rename Extra condition
    prompt = 'What would you like to name the Extra condition?';
    Name.Extra = inputdlg(prompt);
    
    % Get data for Extra condition
    if strcmp(Type.Extra,'Full Kinematics') == 1 % if full Kinematics
        % Get Extra data
        ExtraData = GetKinematics;
        % GDI Computations
        GDIRep.Extra = GDI_Calculator_Fast(0,0,0,ExtraData(1).RepCycle); % compute RC GDI
        % compute GDIs for every cycle within every trial
        CycleCount = 1;
        for i = 1:length(ExtraData)% loop through all trials
            MinCycles = min(ExtraData(i).NumCycles); % find minimum # of cycles on each side for each trial
            for j = 1:MinCycles % loop through # of cycles
                ExtraMatrix(:,:,CycleCount) = [ExtraData(i).Trials(1).data.left(:,j), ExtraData(i).Trials(2).data.left(:,j), ExtraData(i).Trials(3).data.left(:,j),... % Left Pelvis
                    ExtraData(i).Trials(1).data.right(:,j), ExtraData(i).Trials(2).data.right(:,j), ExtraData(i).Trials(3).data.right(:,j),... % Right Pelvis
                    ExtraData(i).Trials(4).data.left(:,j), ExtraData(i).Trials(5).data.left(:,j), ExtraData(i).Trials(6).data.left(:,j),... % Left Hip
                    ExtraData(i).Trials(4).data.right(:,j), ExtraData(i).Trials(5).data.right(:,j), ExtraData(i).Trials(6).data.right(:,j),... % Right Hip
                    ExtraData(i).Trials(7).data.left(:,j), ExtraData(i).Trials(8).data.left(:,j), ExtraData(i).Trials(9).data.left(:,j),... % Left Knee
                    ExtraData(i).Trials(7).data.right(:,j), ExtraData(i).Trials(8).data.right(:,j), ExtraData(i).Trials(9).data.right(:,j),... % Right Knee
                    ExtraData(i).Trials(10).data.left(:,j), ExtraData(i).Trials(11).data.left(:,j), ExtraData(i).Trials(12).data.left(:,j),... % Left Ankle
                    ExtraData(i).Trials(10).data.right(:,j), ExtraData(i).Trials(11).data.right(:,j), ExtraData(i).Trials(12).data.right(:,j)]; % Right Ankle
                CycleCount = CycleCount + 1; % advance to save next cycle in next spot
            end
        end
        [~,~,NumGDITrials] = size(ExtraMatrix); % find total # of cycles saved
        for i = 1:NumGDITrials
            TrialsGDI.Extra(i,:) = GDI_Calculator_Fast(0,0,0,ExtraMatrix(:,:,i)); % compute GDI
        end
        ExtraData(1).TrialsGDI = TrialsGDI.Extra;  % save GDI data into kinematics structure
        ExtraData(1).KineMatrix = ExtraMatrix;  % save cycle kinematics data into kinematics structure
        
    else
        [files2get, ~] = uigetfile('*.c3d*', 'Please select all C3D file(s) to be analyzed', 'MultiSelect', 'on');
        ExtraData = GAMSmini_BatchGaitMeasures(files2get,str2double(char(age)));
        GDIRep.Extra = [0 0];
        ExtraData(1).TrialsGDI = [0 0 0]; 
        for i = 1:length(ExtraData)
            ExtraData(i).GDI = [0 0 0 0 0 0 0];
            if length(ExtraData) > 1
                ExtraData(i).FileName = files2get{i};
            else
                ExtraData(i).FileName = files2get;
            end
        end
    end
    
    % Define Extra variables to be plotted
    for i = 1:length(ExtraData)
        % Overall Gait Performance
        Cadence.Extra.Measures(i) = ExtraData(i).TempSpat(1).data;
        GaitSpeed.Extra.Measures(i) = ExtraData(i).TempSpat(2).data;
        % GDI already exists in ExtraData(1).GDI
        % Spatial Data
        StrideLength.Extra.Measures(i) = ExtraData(i).TempSpat(3).data;
        StepLength.Extra.Left.Measures(i) = ExtraData(i).TempSpat(6).data;
        StepLength.Extra.Right.Measures(i) = ExtraData(i).TempSpat(15).data;
        StepWidth.Extra.AvgMeasures(i) = ExtraData(i).TempSpat(22).data;
        % Temporal Data
        Stance.Extra.Left.Measures(i) = ExtraData(i).TempSpat(9).data;
        Stance.Extra.Right.Measures(i) = ExtraData(i).TempSpat(18).data;
        Swing.Extra.Left.Measures(i) = 100 - Stance.Extra.Left.Measures(i);
        Swing.Extra.Right.Measures(i) = 100 - Stance.Extra.Right.Measures(i);
        DS1.Extra.Left.Measures(i) = ExtraData(i).TempSpat(10).data;
        DS1.Extra.Right.Measures(i) = ExtraData(i).TempSpat(19).data;
        SLS.Extra.Left.Measures(i) = ExtraData(i).TempSpat(11).data;
        SLS.Extra.Right.Measures(i) = ExtraData(i).TempSpat(20).data;
        DS2.Extra.Left.Measures(i) = ExtraData(i).TempSpat(12).data;
        DS2.Extra.Right.Measures(i) = ExtraData(i).TempSpat(21).data;
        % STDs
        StrideLength.Extra.StdMeasures(i) = ExtraData(i).TempSpat(3).data;
        StepLength.Extra.Left.StdMeasures(i) = ExtraData(i).TempSpat(24).data;
        StepLength.Extra.Right.StdMeasures(i) = ExtraData(i).TempSpat(29).data;
        StepWidth.Extra.StdMeasures(i) = ExtraData(i).TempSpat(23).data;
        % Temporal Data
        Stance.Extra.Left.StdMeasures(i) = ExtraData(i).TempSpat(25).data;
        Stance.Extra.Right.StdMeasures(i) = ExtraData(i).TempSpat(30).data;
        Swing.Extra.Left.StdMeasures(i) = Stance.Extra.Left.StdMeasures(i);
        Swing.Extra.Right.StdMeasures(i) = Stance.Extra.Right.StdMeasures(i);
        DS1.Extra.Left.StdMeasures(i) = ExtraData(i).TempSpat(26).data;
        DS1.Extra.Right.StdMeasures(i) = ExtraData(i).TempSpat(31).data;
        SLS.Extra.Left.StdMeasures(i) = ExtraData(i).TempSpat(27).data;
        SLS.Extra.Right.StdMeasures(i) = ExtraData(i).TempSpat(32).data;
        DS2.Extra.Left.StdMeasures(i) = ExtraData(i).TempSpat(28).data;
        DS2.Extra.Right.StdMeasures(i) = ExtraData(i).TempSpat(33).data;
    end
    
    % Calculate averages and STDs
    % Overall Gait Performance
    Cadence.Extra.Avg = mean(Cadence.Extra.Measures);
    Cadence.Extra.STD = std(Cadence.Extra.Measures);
    GaitSpeed.Extra.Avg = mean(GaitSpeed.Extra.Measures);
    GaitSpeed.Extra.STD = std(GaitSpeed.Extra.Measures);
    % Spatial Data
    StrideLength.Extra.Avg = mean(StrideLength.Extra.Measures);
    StrideLength.Extra.STD = std(StrideLength.Extra.Measures);
    StepLength.Extra.Left.Avg = mean(StepLength.Extra.Left.Measures);
    StepLength.Extra.Left.STD = mean(StepLength.Extra.Left.StdMeasures);
    StepLength.Extra.Right.Avg = mean(StepLength.Extra.Right.Measures);
    StepLength.Extra.Right.STD = mean(StepLength.Extra.Right.StdMeasures);
    StepWidth.Extra.Avg = mean(StepWidth.Extra.AvgMeasures);
    StepWidth.Extra.STD = mean(StepWidth.Extra.StdMeasures);
    % Temporal Data
    Stance.Extra.Left.Avg = mean(Stance.Extra.Left.Measures);
    Stance.Extra.Left.STD = mean(Stance.Extra.Left.StdMeasures);
    Stance.Extra.Right.Avg = mean(Stance.Extra.Right.Measures);
    Stance.Extra.Right.STD = mean(Stance.Extra.Right.StdMeasures);
    Swing.Extra.Left.Avg = mean(Swing.Extra.Left.Measures);
    Swing.Extra.Left.STD = mean(Swing.Extra.Left.StdMeasures);
    Swing.Extra.Right.Avg = mean(Swing.Extra.Right.Measures);
    Swing.Extra.Right.STD = mean(Swing.Extra.Right.StdMeasures);
    
    DS1.Extra.Left.Avg = mean(DS1.Extra.Left.Measures);
    DS1.Extra.Left.STD = mean(DS1.Extra.Left.StdMeasures);
    DS1.Extra.Right.Avg = mean(DS1.Extra.Right.Measures);
    DS1.Extra.Right.STD = mean(DS1.Extra.Right.StdMeasures);
    SLS.Extra.Left.Avg = mean(SLS.Extra.Left.Measures);
    SLS.Extra.Left.STD = mean(SLS.Extra.Left.StdMeasures);
    SLS.Extra.Right.Avg = mean(SLS.Extra.Right.Measures);
    SLS.Extra.Right.STD = mean(SLS.Extra.Right.StdMeasures);
    DS2.Extra.Left.Avg = mean(DS2.Extra.Left.Measures);
    DS2.Extra.Left.STD = mean(DS2.Extra.Left.StdMeasures);
    DS2.Extra.Right.Avg = mean(DS2.Extra.Right.Measures);
    DS2.Extra.Right.STD = mean(DS2.Extra.Right.StdMeasures);
    
    if strcmp(Type.Extra, 'Full Kinematics') == 1
    % define rep trial variables
    RepTrial = ExtraData(1).RepTrial.Num(1);
    % Calculate averages and STDs
    % Overall Gait Performance
    Cadence.Extra.Rep = ExtraData(RepTrial).TempSpat(1).data;
    GaitSpeed.Extra.Rep = ExtraData(RepTrial).TempSpat(2).data;
    % Spatial Data
    StrideLength.Extra.Rep = ExtraData(RepTrial).TempSpat(3).data;
    StepLength.Extra.Left.Rep = ExtraData(RepTrial).TempSpat(6).data;
    StepLength.Extra.Right.Rep = ExtraData(RepTrial).TempSpat(15).data;
    StepWidth.Extra.Rep = ExtraData(RepTrial).TempSpat(22).data;
    % Temporal Data
    Stance.Extra.Left.Rep = ExtraData(RepTrial).TempSpat(9).data;
    Stance.Extra.Right.Rep = ExtraData(RepTrial).TempSpat(18).data;
    Swing.Extra.Left.Rep = 100 - ExtraData(RepTrial).TempSpat(9).data;
    Swing.Extra.Right.Rep = 100 - ExtraData(RepTrial).TempSpat(18).data;
    
    DS1.Extra.Left.Rep = ExtraData(RepTrial).TempSpat(10).data;
    DS1.Extra.Right.Rep = ExtraData(RepTrial).TempSpat(19).data;
    SLS.Extra.Left.Rep = ExtraData(RepTrial).TempSpat(11).data;
    SLS.Extra.Right.Rep = ExtraData(RepTrial).TempSpat(20).data;
    DS2.Extra.Left.Rep = ExtraData(RepTrial).TempSpat(12).data;
    DS2.Extra.Right.Rep = ExtraData(RepTrial).TempSpat(21).data;
    end
end
clearvars GDIplus GDIminus DiffPlus DiffMinus GDIavg Diff

%% Pre-define colors and measures based on # of desired plots
if NumConds == 1 % if 1 condition
    XAxis = {char(Name.Curr)};
    Blue = rgb('DeepSkyBlue');
    RedGreen = [rgb('Crimson'); rgb('ForestGreen')];
    Gait_Speed = GaitSpeed.Current.Avg/TSN.TimeDist(3)*100;
    Cadence_ = Cadence.Current.Avg/TSN.TimeDist(1)*100;
    Stride_Length = StrideLength.Current.Avg/TSN.TimeDist(2)*100;
    Step_Length = [StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100];
    Step_Width = StepWidth.Current.Avg/TSN.TimeDist(6)*100;
    Stance_Phase = [Stance.Current.Left.Avg, Stance.Current.Right.Avg];
    Swing_Phase = [Swing.Current.Left.Avg, Swing.Current.Right.Avg];
    Initial_Double_Support = [DS1.Current.Left.Avg, DS1.Current.Right.Avg];
    Single_Limb_Support = [SLS.Current.Left.Avg, SLS.Current.Right.Avg];
    Secondary_Double_Support = [DS2.Current.Left.Avg, DS2.Current.Right.Avg];
elseif NumConds == 2 % if 2 conditions
    Blue = [rgb('RoyalBlue'); rgb('Cyan')];
    RedGreen = [rgb('Crimson'); rgb('ForestGreen'); rgb('LightCoral')*.85;rgb('LightGreen')*.85];
    if strcmp(AddAFO,'Yes') == 1 % if AFO is 2nd condition
        XAxis = {char(Name.Curr),char(Name.AFO)};
        Gait_Speed = [GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100];
        Cadence_ = [Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100];
        Stride_Length = [StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100];
        Step_Length = [StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100];
        Step_Width = [StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100];
        Stance_Phase = [Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg];
        Swing_Phase = [Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg];
        Initial_Double_Support = [DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg];
        Single_Limb_Support = [SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg];
        Secondary_Double_Support = [DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg];
    elseif strcmp(AddLast,'Yes') == 1 % if last is 2nd condition
        XAxis = {char(Name.Curr),char(Name.Last)};
        Gait_Speed = [GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Last.Avg/TSN.TimeDist(3)*100];
        Cadence_ = [Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100];
        Stride_Length = [StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100];
        Step_Length = [StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100];
        Step_Width = [StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100];
        Stance_Phase = [Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.Last.Left.Avg, Stance.Last.Right.Avg];
        Swing_Phase = [Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.Last.Left.Avg, Swing.Last.Right.Avg];
        Initial_Double_Support = [DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.Last.Left.Avg, DS1.Last.Right.Avg];
        Single_Limb_Support = [SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.Last.Left.Avg, SLS.Last.Right.Avg];
        Secondary_Double_Support = [DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.Last.Left.Avg, DS2.Last.Right.Avg];
    elseif strcmp(AddExtra,'Yes')  == 1 % if extra is 2nd condition
        XAxis = {char(Name.Curr),char(Name.Extra)};
        Gait_Speed = [GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100];
        Cadence_ = [Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100];
        Stride_Length = [StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100];
        Step_Length = [StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100];
        Step_Width = [StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100];
        Stance_Phase = [Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg];
        Swing_Phase = [Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg];
        Initial_Double_Support = [DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg];
        Single_Limb_Support = [SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg];
        Secondary_Double_Support = [DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg];
    end
elseif NumConds == 3 % if 3 conditions
    Blue = [rgb('RoyalBlue');rgb('DeepSkyBlue');rgb('Cyan')];
    RedGreen = [rgb('Crimson'); rgb('ForestGreen');rgb('Red');rgb('Green');rgb('LightCoral')*.85;rgb('LightGreen')*.85];
    if strcmp(AddAFO,'Yes') == 1 && strcmp(AddLast,'Yes') == 1 % if AFO and Last are 2nd and 3rd conditions
        XAxis = {char(Name.Curr), char(Name.AFO), char(Name.Last)};
        Gait_Speed = [GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100,GaitSpeed.Last.Avg/TSN.TimeDist(3)*100];
        Cadence_ = [Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100];
        Stride_Length = [StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100];
        Step_Length = [StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
            StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100];
        Step_Width = [StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100];
        Stance_Phase = [Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg,...
            Stance.Last.Left.Avg, Stance.Last.Right.Avg];
        Swing_Phase = [Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg...
            Swing.Last.Left.Avg, Swing.Last.Right.Avg];
        Initial_Double_Support = [DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg,...
            DS1.Last.Left.Avg, DS1.Last.Right.Avg];
        Single_Limb_Support = [SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg,...
            SLS.Last.Left.Avg, SLS.Last.Right.Avg];
        Secondary_Double_Support = [DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg,...
            DS2.Last.Left.Avg, DS2.Last.Right.Avg];
    elseif strcmp(AddAFO,'Yes') == 1 && strcmp(AddExtra,'Yes') == 1 % if AFO and Extra are 2nd and 3rd conditions
        XAxis = {char(Name.Curr), char(Name.AFO), char(Name.Extra)};
        Gait_Speed = [GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100,GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100];
        Cadence_ = [Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100];
        Stride_Length = [StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100];
        Step_Length = [StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
            StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100];
        Step_Width = [StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100];
        Stance_Phase = [Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg,...
            Stance.Extra.Left.Avg, Stance.Extra.Right.Avg];
        Swing_Phase = [Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg...
            Swing.Extra.Left.Avg, Swing.Extra.Right.Avg];
        Initial_Double_Support = [DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg,...
            DS1.Extra.Left.Avg, DS1.Extra.Right.Avg];
        Single_Limb_Support = [SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg,...
            SLS.Extra.Left.Avg, SLS.Extra.Right.Avg];
        Secondary_Double_Support = [DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg,...
            DS2.Extra.Left.Avg, DS2.Extra.Right.Avg];
    elseif strcmp(AddLast,'Yes') == 1 && strcmp(AddExtra,'Yes') == 1 % if Extra and Last are 2nd and 3rd conditions
        XAxis = {char(Name.Curr), char(Name.Last), char(Name.Extra)};
        Gait_Speed = [GaitSpeed.Current.Avg/TSN.TimeDist(3)*100,GaitSpeed.Last.Avg/TSN.TimeDist(3)*100,GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100];
        Cadence_ = [Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100];
        Stride_Length = [StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100];
        Step_Length = [StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100,...
            StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100,StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100];
        Step_Width = [StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100];
        Stance_Phase = [Stance.Current.Left.Avg, Stance.Current.Right.Avg,...
            Stance.Last.Left.Avg, Stance.Last.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg];
        Swing_Phase = [Swing.Current.Left.Avg, Swing.Current.Right.Avg,...
            Swing.Last.Left.Avg, Swing.Last.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg];
        Initial_Double_Support = [DS1.Current.Left.Avg, DS1.Current.Right.Avg,...
            DS1.Last.Left.Avg, DS1.Last.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg];
        Single_Limb_Support = [SLS.Current.Left.Avg, SLS.Current.Right.Avg,...
            SLS.Last.Left.Avg, SLS.Last.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg];
        Secondary_Double_Support = [DS2.Current.Left.Avg, DS2.Current.Right.Avg,...
            DS2.Last.Left.Avg, DS2.Last.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg];
    end
elseif NumConds == 4
    XAxis = {char(Name.Curr), char(Name.AFO), char(Name.Last), char(Name.Extra)};
    Blue = [rgb('RoyalBlue');rgb('DeepSkyBlue');rgb('Cyan');rgb('PowderBlue')];
    RedGreen = [rgb('Crimson'); rgb('ForestGreen');rgb('Red');rgb('Green');rgb('LightCoral')*.85;rgb('LightGreen')*.85; rgb('LightSalmon'); rgb('PaleGreen')];
    Gait_Speed = [GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100,GaitSpeed.Last.Avg/TSN.TimeDist(3)*100,GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100];
    Cadence_ = [Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100];
    Stride_Length = [StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100];
    Step_Length = [StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
        StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100,StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100];
    Step_Width = [StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100];
    Stance_Phase = [Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg,...
        Stance.Last.Left.Avg, Stance.Last.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg];
    Swing_Phase = [Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg...
        Swing.Last.Left.Avg, Swing.Last.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg];
    Initial_Double_Support = [DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg,...
        DS1.Last.Left.Avg, DS1.Last.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg];
    Single_Limb_Support = [SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg,...
        SLS.Last.Left.Avg, SLS.Last.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg];
    Secondary_Double_Support = [DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg,...
        DS2.Last.Left.Avg, DS2.Last.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg];
end


%% Plot  GAMS outputs with bullet graphs
GAMSplots = figure( 'Position', [100, 100, 1000, 500]);
% gait speed
subplot(251);
BulletGraph('V', 75, 90, Gait_Speed, 100, [0 150], Blue,NumConds, 'Yes', 'No');
grid on;    hold on;
% cadence
subplot(252);
BulletGraph('V', 75, 90, Cadence_ , 100, [0 150], Blue,NumConds, 'Yes', 'No');
hold on;  grid on;
% Stride Length
subplot(253);
BulletGraph('V', 75, 90, Stride_Length , 100, [0 150], Blue,NumConds, 'Yes', 'No');
hold on; grid on;
% step length
subplot(254);
BulletGraph('V', 75, 90, Step_Length, 100, [0 150], RedGreen,NumConds*2, 'Yes', 'Yes');
hold on; grid on;
% step width
subplot(255);
BulletGraph('ReverseV', 150, 125, Step_Width ,100, [0 300], Blue,NumConds,'Yes', 'No');
hold on; grid on;
% Stance Phase
subplot(256);
BulletGraph('V', 50, 58, Stance_Phase, 62, [0 100],RedGreen,NumConds*2, 'Yes', 'Yes');
hold on; grid on;
% Swing Phase
subplot(257);
BulletGraph('V', 25, 32, Swing_Phase, 38, [0 100],RedGreen,NumConds*2, 'Yes', 'Yes');
hold on; grid on;
% initial double support
subplot(258);
BulletGraph('ReverseV', 20, 15, Initial_Double_Support, 12, [0 50],RedGreen,NumConds*2, 'Yes', 'Yes');
hold on; grid on;
%  single limb support
subplot(259);
BulletGraph('V', 25, 32, Single_Limb_Support, 38, [0 50],RedGreen,NumConds*2, 'Yes', 'Yes');
hold on; grid on;
%  secondary double support
subplot(2,5,10);
BulletGraph('ReverseV', 20, 15, Secondary_Double_Support, 12, [0 50],RedGreen,NumConds*2, 'Yes', 'Yes');
hold on; grid on;

%% Add errorbars to bar plots
% prompt = 'Would you like to add errorbars?';
% ErrorBar = questdlg(prompt);
% if strcmp(ErrorBar,'Yes') == 1
Version = ver;
if datenum(Version(1).Date) > datenum(2017,2,16) % if the matlab version is after 2017a, plot errorbars with 'capsize' specifics
    if NumConds == 1
        X1 = 1;
        X2 = [0.9 1.1];
        subplot(251); % gait speed
        e = errorbar(1,GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Current.STD/TSN.TimeDist(3)*100, '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
        subplot(252); % cadence
        errorbar(1,Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Current.STD/TSN.TimeDist(1)*100, '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
        subplot(253); % Stride Length
        errorbar(1,StrideLength.Current.Avg/TSN.TimeDist(2)*100,StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Current.STD/TSN.TimeDist(2)*100, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(254); % step length
        errorbar(0.9,StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Left.STD/TSN.TimeDist(4)*100, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        errorbar(1.1,StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(255); % step width
        errorbar(1,StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Current.STD/TSN.TimeDist(6)*100, '.k','LineWidth',1,'CapSize',10,'Marker','none');
        subplot(256); % Stance Phase
        errorbar(0.9,Stance.Current.Left.Avg, Stance.Current.Left.STD, Stance.Current.Left.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        errorbar(1.1,Stance.Current.Right.Avg, Stance.Current.Right.STD, Stance.Current.Right.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(257); % Swing Phase
        errorbar(0.9,Swing.Current.Left.Avg, Swing.Current.Left.STD, Swing.Current.Left.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        errorbar(1.1,Swing.Current.Right.Avg, Swing.Current.Right.STD, Swing.Current.Right.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(258); % initial double support
        errorbar(0.9,DS1.Current.Left.Avg, DS1.Current.Left.STD, DS1.Current.Left.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        errorbar(1.1,DS1.Current.Right.Avg, DS1.Current.Right.STD, DS1.Current.Right.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(259); %  single limb support
        errorbar(0.9,SLS.Current.Left.Avg, SLS.Current.Left.STD, SLS.Current.Left.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        errorbar(1.1,SLS.Current.Right.Avg, SLS.Current.Right.STD, SLS.Current.Right.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(2,5,10); %  secondary double support
        errorbar(0.9,DS2.Current.Left.Avg, DS2.Current.Left.STD, DS2.Current.Left.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        errorbar(1.1,DS2.Current.Right.Avg, DS2.Current.Right.STD, DS2.Current.Right.STD, '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
    elseif NumConds == 2
        % define X locations
        X1 =  [0.85 1.15];
        X2 = [0.75 0.9 1.1 1.25];
        if strcmp(AddAFO,'Yes') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.AFO.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        elseif strcmp(AddLast,'Yes') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Last.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.Last.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.Last.Left.Avg, Stance.Last.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.Last.Left.Avg, Swing.Last.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.Last.Left.Avg, DS1.Last.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.Last.Left.Avg, SLS.Last.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.Last.Left.Avg, DS2.Last.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        elseif strcmp(AddExtra,'Yes') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        end
    elseif NumConds == 3
        % define X locations
        X1 =  [0.75 1 1.25];
        X2 = [0.6 0.725 0.9375 1.0625 1.275 1.4];
        if strcmp(AddExtra,'No') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100,GaitSpeed.Last.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.AFO.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100],....
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
                StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg, Stance.Last.Left.Avg, Stance.Last.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg, Swing.Last.Left.Avg, Swing.Last.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg, DS1.Last.Left.Avg, DS1.Last.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg, SLS.Last.Left.Avg, SLS.Last.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg, DS2.Last.Left.Avg, DS2.Last.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        elseif strcmp(AddAFO,'No') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Last.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.Last.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],....
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100,...
                StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.Last.Left.Avg, Stance.Last.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.Last.Left.Avg, Swing.Last.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.Last.Left.Avg, DS1.Last.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.Last.Left.Avg, SLS.Last.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.Last.Left.Avg, DS2.Last.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        elseif strcmp(AddLast,'No') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.AFO.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],....
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
                StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        end
    elseif NumConds == 4
        % define X locations
        X1 = [0.7 0.9 1.1 1.3];
        X2 = [0.49 0.61 0.79 0.91 1.09 1.21 1.39 1.51];
        subplot(251); % gait speed
        errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100,...
            GaitSpeed.Last.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100],...
            [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100,...
            GaitSpeed.Last.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
            [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100,...
            GaitSpeed.Last.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
            '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
        subplot(252); % cadence
        errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100,...
            Cadence.Last.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100],...
            [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100,...
            Cadence.Last.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
            [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100,...
            Cadence.Last.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
            '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
        subplot(253); % Stride Length
        errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100,...
            StrideLength.Last.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100],...
            [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100,...
            StrideLength.Last.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
            [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.AFO.STD/TSN.TimeDist(2)*100,...
            StrideLength.Last.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
            '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(254); % step length
        errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
            StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
            StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100,...
            StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100],...
            [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
            [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
            '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(255); % step width
        errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100,...
            StepWidth.Last.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100],...
            [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100,...
            StepWidth.Last.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
            [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100,...
            StepWidth.Last.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
            '.k','LineWidth',1, 'CapSize',10, 'Marker','none');
        subplot(256); % Stance Phase
        errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg,...
            Stance.Last.Left.Avg, Stance.Last.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg],...
            [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, ...
            Stance.Last.Left.STD, Stance.Last.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
            [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, ...
            Stance.Last.Left.STD, Stance.Last.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
            '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(257); % Swing Phase
        errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg,...
            Swing.Last.Left.Avg, Swing.Last.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg],...
            [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, ...
            Swing.Last.Left.STD, Swing.Last.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
            [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, ...
            Swing.Last.Left.STD, Swing.Last.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
            '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(258); % initial double support
        errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg,...
            DS1.Last.Left.Avg, DS1.Last.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg],...
            [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, ...
            DS1.Last.Left.STD, DS1.Last.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
            [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, ...
            DS1.Last.Left.STD, DS1.Last.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
            '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(259); %  single limb support
        errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg,...
            SLS.Last.Left.Avg, SLS.Last.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg],...
            [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, ...
            SLS.Last.Left.STD, SLS.Last.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
            [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, ...
            SLS.Last.Left.STD, SLS.Last.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
            '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
        subplot(2,5,10); %  secondary double support
        errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg,...
            DS2.Last.Left.Avg, DS2.Last.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg],...
            [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, ...
            DS2.Last.Left.STD, DS2.Last.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
            [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, ...
            DS2.Last.Left.STD, DS2.Last.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
            '.k','LineWidth',1, 'CapSize',8, 'Marker','none');
    end
else % if version if before 2017a, plot with no capsize edits
    %         if strcmp(ErrorBar,'Yes') == 1
    if NumConds == 1
        subplot(251); % gait speed
        errorbar(1,GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Current.STD/TSN.TimeDist(3)*100, '.k','LineWidth',1,  'Marker','none');
        subplot(252); % cadence
        errorbar(1,Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Current.STD/TSN.TimeDist(1)*100, '.k','LineWidth',1,  'Marker','none');
        subplot(253); % Stride Length
        errorbar(1,StrideLength.Current.Avg/TSN.TimeDist(2)*100,StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Current.STD/TSN.TimeDist(2)*100, '.k','LineWidth',1, 'Marker','none');
        subplot(254); % step length
        errorbar(0.9,StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Left.STD/TSN.TimeDist(4)*100, '.k','LineWidth',1,  'Marker','none');
        errorbar(1.1,StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, '.k','LineWidth',1,  'Marker','none');
        subplot(255); % step width
        errorbar(1,StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Current.STD/TSN.TimeDist(6)*100, '.k','LineWidth',1,'Marker','none');
        subplot(256); % Stance Phase
        errorbar(0.9,Stance.Current.Left.Avg, Stance.Current.Left.STD, Stance.Current.Left.STD, '.k','LineWidth',1,  'Marker','none');
        errorbar(1.1,Stance.Current.Right.Avg, Stance.Current.Right.STD, Stance.Current.Right.STD, '.k','LineWidth',1,  'Marker','none');
        subplot(257); % Swing Phase
        errorbar(0.9,Swing.Current.Left.Avg, Swing.Current.Left.STD, Swing.Current.Left.STD, '.k','LineWidth',1,  'Marker','none');
        errorbar(1.1,Swing.Current.Right.Avg, Swing.Current.Right.STD, Swing.Current.Right.STD, '.k','LineWidth',1,  'Marker','none');
        subplot(258); % initial double support
        errorbar(0.9,DS1.Current.Left.Avg, DS1.Current.Left.STD, DS1.Current.Left.STD, '.k','LineWidth',1,  'Marker','none');
        errorbar(1.1,DS1.Current.Right.Avg, DS1.Current.Right.STD, DS1.Current.Right.STD, '.k','LineWidth',1,  'Marker','none');
        subplot(259); %  single limb support
        errorbar(0.9,SLS.Current.Left.Avg, SLS.Current.Left.STD, SLS.Current.Left.STD, '.k','LineWidth',1,  'Marker','none');
        errorbar(1.1,SLS.Current.Right.Avg, SLS.Current.Right.STD, SLS.Current.Right.STD, '.k','LineWidth',1,  'Marker','none');
        subplot(2,5,10); %  secondary double support
        errorbar(0.9,DS2.Current.Left.Avg, DS2.Current.Left.STD, DS2.Current.Left.STD, '.k','LineWidth',1,  'Marker','none');
        errorbar(1.1,DS2.Current.Right.Avg, DS2.Current.Right.STD, DS2.Current.Right.STD, '.k','LineWidth',1,  'Marker','none');
    elseif NumConds == 2
        % define X locations
        X1 =  [0.85 1.15];
        X2 = [0.75 0.9 1.1 1.25];
        if strcmp(AddAFO,'Yes') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.AFO.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
        elseif strcmp(AddLast,'Yes') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Last.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.Last.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.Last.Left.Avg, Stance.Last.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.Last.Left.Avg, Swing.Last.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.Last.Left.Avg, DS1.Last.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.Last.Left.Avg, SLS.Last.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.Last.Left.Avg, DS2.Last.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
        elseif strcmp(AddExtra,'Yes') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(253); % Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100,StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
        end
    elseif NumConds == 3
        % define X locations
        X1 =  [0.75 1 1.25];
        X2 = [0.6 0.725 0.9375 1.0625 1.275 1.4];
        if strcmp(AddExtra,'No') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100,GaitSpeed.Last.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(253); %Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
                StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg, Stance.Last.Left.Avg, Stance.Last.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg, Swing.Last.Left.Avg, Swing.Last.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg, DS1.Last.Left.Avg, DS1.Last.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg, SLS.Last.Left.Avg, SLS.Last.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg, DS2.Last.Left.Avg, DS2.Last.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
        elseif strcmp(AddAFO,'No') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.Last.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.Last.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.Last.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.Last.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(253); %Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.Last.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.Last.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100,...
                StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.Last.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.Last.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.Last.Left.Avg, Stance.Last.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.Last.Left.STD, Stance.Last.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.Last.Left.Avg, Swing.Last.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.Last.Left.STD, Swing.Last.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.Last.Left.Avg, DS1.Last.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.Last.Left.STD, DS1.Last.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.Last.Left.Avg, SLS.Last.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.Last.Left.STD, SLS.Last.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.Last.Left.Avg, DS2.Last.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.Last.Left.STD, DS2.Last.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
        elseif strcmp(AddLast,'No') == 1
            subplot(251); % gait speed
            errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(252); % cadence
            errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(253); %Stride Length
            errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
                '.k','LineWidth',1, 'Marker','none');
            subplot(254); % step length
            errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
                StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100],...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
                StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(255); % step width
            errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
                '.k','LineWidth',1,  'Marker','none');
            subplot(256); % Stance Phase
            errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg],...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(257); % Swing Phase
            errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg],...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(258); % initial double support
            errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg],...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(259); %  single limb support
            errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg],...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
            subplot(2,5,10); %  secondary double support
            errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg],...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
                '.k','LineWidth',1,  'Marker','none');
        end
    elseif NumConds == 4
        % define X locations
        X1 = [0.7 0.9 1.1 1.3];
        X2 = [0.49 0.61 0.79 0.91 1.09 1.21 1.39 1.51];
        subplot(251); % gait speed
        errorbar(X1,[GaitSpeed.Current.Avg/TSN.TimeDist(3)*100, GaitSpeed.AFO.Avg/TSN.TimeDist(3)*100,...
            GaitSpeed.Last.Avg/TSN.TimeDist(3)*100, GaitSpeed.Extra.Avg/TSN.TimeDist(3)*100],...
            [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100,...
            GaitSpeed.Last.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
            [GaitSpeed.Current.STD/TSN.TimeDist(3)*100, GaitSpeed.AFO.STD/TSN.TimeDist(3)*100,...
            GaitSpeed.Last.STD/TSN.TimeDist(3)*100, GaitSpeed.Extra.STD/TSN.TimeDist(3)*100],...
            '.k','LineWidth',1,  'Marker','none');
        subplot(252); % cadence
        errorbar(X1,[Cadence.Current.Avg/TSN.TimeDist(1)*100, Cadence.AFO.Avg/TSN.TimeDist(1)*100,...
            Cadence.Last.Avg/TSN.TimeDist(1)*100, Cadence.Extra.Avg/TSN.TimeDist(1)*100],...
            [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100,...
            Cadence.Last.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
            [Cadence.Current.STD/TSN.TimeDist(1)*100, Cadence.AFO.STD/TSN.TimeDist(1)*100,...
            Cadence.Last.STD/TSN.TimeDist(1)*100, Cadence.Extra.STD/TSN.TimeDist(1)*100],...
            '.k','LineWidth',1,  'Marker','none');
        subplot(253); %Stride Length
        errorbar(X1,[StrideLength.Current.Avg/TSN.TimeDist(2)*100, StrideLength.AFO.Avg/TSN.TimeDist(2)*100,...
            StrideLength.Last.Avg/TSN.TimeDist(2)*100, StrideLength.Extra.Avg/TSN.TimeDist(2)*100],...
            [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100,...
            StrideLength.Last.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
            [StrideLength.Current.STD/TSN.TimeDist(2)*100, StrideLength.AFO.STD/TSN.TimeDist(2)*100,...
            StrideLength.Last.STD/TSN.TimeDist(2)*100, StrideLength.Extra.STD/TSN.TimeDist(2)*100],...
            '.k','LineWidth',1, 'Marker','none');
        subplot(254); % step length
        errorbar(X2,[StepLength.Current.Left.Avg/TSN.TimeDist(4)*100, StepLength.Current.Right.Avg/TSN.TimeDist(4)*100, ...
            StepLength.AFO.Left.Avg/TSN.TimeDist(4)*100, StepLength.AFO.Right.Avg/TSN.TimeDist(4)*100,...
            StepLength.Last.Left.Avg/TSN.TimeDist(4)*100, StepLength.Last.Right.Avg/TSN.TimeDist(4)*100,...
            StepLength.Extra.Left.Avg/TSN.TimeDist(4)*100, StepLength.Extra.Right.Avg/TSN.TimeDist(4)*100],...
            [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
            [StepLength.Current.Left.STD/TSN.TimeDist(4)*100, StepLength.Current.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.AFO.Left.STD/TSN.TimeDist(4)*100, StepLength.AFO.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.Last.Left.STD/TSN.TimeDist(4)*100, StepLength.Last.Right.STD/TSN.TimeDist(4)*100, ...
            StepLength.Extra.Left.STD/TSN.TimeDist(4)*100, StepLength.Extra.Right.STD/TSN.TimeDist(4)*100], ...
            '.k','LineWidth',1,  'Marker','none');
        subplot(255); % step width
        errorbar(X1,[StepWidth.Current.Avg/TSN.TimeDist(6)*100, StepWidth.AFO.Avg/TSN.TimeDist(6)*100,...
            StepWidth.Last.Avg/TSN.TimeDist(6)*100, StepWidth.Extra.Avg/TSN.TimeDist(6)*100],...
            [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100,...
            StepWidth.Last.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
            [StepWidth.Current.STD/TSN.TimeDist(6)*100, StepWidth.AFO.STD/TSN.TimeDist(6)*100,...
            StepWidth.Last.STD/TSN.TimeDist(6)*100, StepWidth.Extra.STD/TSN.TimeDist(6)*100],...
            '.k','LineWidth',1,  'Marker','none');
        subplot(256); % Stance Phase
        errorbar(X2,[Stance.Current.Left.Avg, Stance.Current.Right.Avg, Stance.AFO.Left.Avg, Stance.AFO.Right.Avg,...
            Stance.Last.Left.Avg, Stance.Last.Right.Avg, Stance.Extra.Left.Avg, Stance.Extra.Right.Avg],...
            [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, ...
            Stance.Last.Left.STD, Stance.Last.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
            [Stance.Current.Left.STD, Stance.Current.Right.STD, Stance.AFO.Left.STD, Stance.AFO.Right.STD, ...
            Stance.Last.Left.STD, Stance.Last.Right.STD, Stance.Extra.Left.STD, Stance.Extra.Right.STD], ...
            '.k','LineWidth',1,  'Marker','none');
        subplot(257); % Swing Phase
        errorbar(X2,[Swing.Current.Left.Avg, Swing.Current.Right.Avg, Swing.AFO.Left.Avg, Swing.AFO.Right.Avg,...
            Swing.Last.Left.Avg, Swing.Last.Right.Avg, Swing.Extra.Left.Avg, Swing.Extra.Right.Avg],...
            [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, ...
            Swing.Last.Left.STD, Swing.Last.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
            [Swing.Current.Left.STD, Swing.Current.Right.STD, Swing.AFO.Left.STD, Swing.AFO.Right.STD, ...
            Swing.Last.Left.STD, Swing.Last.Right.STD, Swing.Extra.Left.STD, Swing.Extra.Right.STD], ...
            '.k','LineWidth',1,  'Marker','none');
        subplot(258); % initial double support
        errorbar(X2,[DS1.Current.Left.Avg, DS1.Current.Right.Avg, DS1.AFO.Left.Avg, DS1.AFO.Right.Avg,...
            DS1.Last.Left.Avg, DS1.Last.Right.Avg, DS1.Extra.Left.Avg, DS1.Extra.Right.Avg],...
            [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, ...
            DS1.Last.Left.STD, DS1.Last.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
            [DS1.Current.Left.STD, DS1.Current.Right.STD, DS1.AFO.Left.STD, DS1.AFO.Right.STD, ...
            DS1.Last.Left.STD, DS1.Last.Right.STD, DS1.Extra.Left.STD, DS1.Extra.Right.STD], ...
            '.k','LineWidth',1,  'Marker','none');
        subplot(259); %  single limb support
        errorbar(X2,[SLS.Current.Left.Avg, SLS.Current.Right.Avg, SLS.AFO.Left.Avg, SLS.AFO.Right.Avg,...
            SLS.Last.Left.Avg, SLS.Last.Right.Avg, SLS.Extra.Left.Avg, SLS.Extra.Right.Avg],...
            [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, ...
            SLS.Last.Left.STD, SLS.Last.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
            [SLS.Current.Left.STD, SLS.Current.Right.STD, SLS.AFO.Left.STD, SLS.AFO.Right.STD, ...
            SLS.Last.Left.STD, SLS.Last.Right.STD, SLS.Extra.Left.STD, SLS.Extra.Right.STD], ...
            '.k','LineWidth',1,  'Marker','none');
        subplot(2,5,10); %  secondary double support
        errorbar(X2,[DS2.Current.Left.Avg, DS2.Current.Right.Avg, DS2.AFO.Left.Avg, DS2.AFO.Right.Avg,...
            DS2.Last.Left.Avg, DS2.Last.Right.Avg, DS2.Extra.Left.Avg, DS2.Extra.Right.Avg],...
            [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, ...
            DS2.Last.Left.STD, DS2.Last.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
            [DS2.Current.Left.STD, DS2.Current.Right.STD, DS2.AFO.Left.STD, DS2.AFO.Right.STD, ...
            DS2.Last.Left.STD, DS2.Last.Right.STD, DS2.Extra.Left.STD, DS2.Extra.Right.STD], ...
            '.k','LineWidth',1,  'Marker','none');
    end
end

%% If type is full kinematics,  define and plot rep trial values as a dot
% if trials are not full kinematics, set rep values to 0
Absent = -10;
if strcmp(Type.Current, 'Full Kinematics') == 0 % current condition 
    GaitSpeed.Current.Rep = Absent;
    Cadence.Current.Rep = Absent;
    StrideLength.Current.Rep = Absent;
    StepLength.Current.Left.Rep = Absent; StepLength.Current.Right.Rep = Absent;
    StepWidth.Current.Rep = Absent;
    Stance.Current.Left.Rep = Absent; Stance.Current.Right.Rep = Absent;
    Swing.Current.Left.Rep = Absent; Swing.Current.Right.Rep = Absent;
    DS1.Current.Left.Rep = Absent; DS1.Current.Right.Rep = Absent;
    SLS.Current.Left.Rep = Absent; SLS.Current.Right.Rep = Absent;
    DS2.Current.Left.Rep = Absent; DS2.Current.Right.Rep = Absent;
end
if strcmp(AddAFO, 'Yes') == 1
    if strcmp(Type.AFO, 'Full Kinematics') == 0 % AFO condition
        GaitSpeed.AFO.Rep = Absent;
        Cadence.AFO.Rep = Absent;
        StrideLength.AFO.Rep = Absent;
        StepLength.AFO.Left.Rep = Absent; StepLength.AFO.Right.Rep = Absent;
        StepWidth.AFO.Rep = Absent;
        Stance.AFO.Left.Rep = Absent; Stance.AFO.Right.Rep = Absent;
        Swing.AFO.Left.Rep = Absent; Swing.AFO.Right.Rep = Absent;
        DS1.AFO.Left.Rep = Absent; DS1.AFO.Right.Rep = Absent;
        SLS.AFO.Left.Rep = Absent; SLS.AFO.Right.Rep = Absent;
        DS2.AFO.Left.Rep = Absent; DS2.AFO.Right.Rep = Absent;
    end
end
if strcmp(AddLast, 'Yes') == 1
    if strcmp(Type.Last, 'Full Kinematics') == 0 % Last condition
        GaitSpeed.Last.Rep = Absent;
        Cadence.Last.Rep = Absent;
        StrideLength.Last.Rep = Absent;
        StepLength.Last.Left.Rep = Absent; StepLength.Last.Right.Rep = Absent;
        StepWidth.Last.Rep = Absent;
        Stance.Last.Left.Rep = Absent; Stance.Last.Right.Rep = Absent;
        Swing.Last.Left.Rep = Absent; Swing.Last.Right.Rep = Absent;
        DS1.Last.Left.Rep = Absent; DS1.Last.Right.Rep = Absent;
        SLS.Last.Left.Rep = Absent; SLS.Last.Right.Rep = Absent;
        DS2.Last.Left.Rep = Absent; DS2.Last.Right.Rep = Absent;
    end
end
if strcmp(AddExtra, 'Yes') == 1
    if strcmp(Type.Extra, 'Full Kinematics') == 0 % Extra condition
        GaitSpeed.Extra.Rep = Absent;
        Cadence.Extra.Rep = Absent;
        StrideLength.Extra.Rep = Absent;
        StepLength.Extra.Left.Rep = Absent; StepLength.Extra.Right.Rep = Absent;
        StepWidth.Extra.Rep = Absent;
        Stance.Extra.Left.Rep = Absent; Stance.Extra.Right.Rep = Absent;
        Swing.Extra.Left.Rep = Absent; Swing.Extra.Right.Rep = Absent;
        DS1.Extra.Left.Rep = Absent; DS1.Extra.Right.Rep = Absent;
        SLS.Extra.Left.Rep = Absent; SLS.Extra.Right.Rep = Absent;
        DS2.Extra.Left.Rep = Absent; DS2.Extra.Right.Rep = Absent;
    end
end

% define REP series
if NumConds == 1 % if 1 condition
    Rep.Gait_Speed = GaitSpeed.Current.Rep/TSN.TimeDist(3)*100;
    Rep.Cadence_ = Cadence.Current.Rep/TSN.TimeDist(1)*100;
    Rep.Stride_Length = StrideLength.Current.Rep/TSN.TimeDist(2)*100;
    Rep.Step_Length = [StepLength.Current.Left.Rep/TSN.TimeDist(4)*100, StepLength.Current.Right.Rep/TSN.TimeDist(4)*100];
    Rep.Step_Width = StepWidth.Current.Rep/TSN.TimeDist(6)*100;
    Rep.Stance_Phase = [Stance.Current.Left.Rep, Stance.Current.Right.Rep];
    Rep.Swing_Phase = [Swing.Current.Left.Rep, Swing.Current.Right.Rep];
    Rep.Initial_Double_Support = [DS1.Current.Left.Rep, DS1.Current.Right.Rep];
    Rep.Single_Limb_Support = [SLS.Current.Left.Rep, SLS.Current.Right.Rep];
    Rep.Secondary_Double_Support = [DS2.Current.Left.Rep, DS2.Current.Right.Rep];
elseif NumConds == 2 % if 2 conditions
    if strcmp(AddAFO,'Yes') == 1 % if AFO is 2nd condition
        Rep.Gait_Speed = [GaitSpeed.Current.Rep/TSN.TimeDist(3)*100, GaitSpeed.AFO.Rep/TSN.TimeDist(3)*100];
        Rep.Cadence_ = [Cadence.Current.Rep/TSN.TimeDist(1)*100, Cadence.AFO.Rep/TSN.TimeDist(1)*100];
        Rep.Stride_Length = [StrideLength.Current.Rep/TSN.TimeDist(2)*100, StrideLength.AFO.Rep/TSN.TimeDist(2)*100];
        Rep.Step_Length = [StepLength.Current.Left.Rep/TSN.TimeDist(4)*100, StepLength.Current.Right.Rep/TSN.TimeDist(4)*100, StepLength.AFO.Left.Rep/TSN.TimeDist(4)*100, StepLength.AFO.Right.Rep/TSN.TimeDist(4)*100];
        Rep.Step_Width = [StepWidth.Current.Rep/TSN.TimeDist(6)*100, StepWidth.AFO.Rep/TSN.TimeDist(6)*100];
        Rep.Stance_Phase = [Stance.Current.Left.Rep, Stance.Current.Right.Rep, Stance.AFO.Left.Rep, Stance.AFO.Right.Rep];
        Rep.Swing_Phase = [Swing.Current.Left.Rep, Swing.Current.Right.Rep, Swing.AFO.Left.Rep, Swing.AFO.Right.Rep];
        Rep.Initial_Double_Support = [DS1.Current.Left.Rep, DS1.Current.Right.Rep, DS1.AFO.Left.Rep, DS1.AFO.Right.Rep];
        Rep.Single_Limb_Support = [SLS.Current.Left.Rep, SLS.Current.Right.Rep, SLS.AFO.Left.Rep, SLS.AFO.Right.Rep];
        Rep.Secondary_Double_Support = [DS2.Current.Left.Rep, DS2.Current.Right.Rep, DS2.AFO.Left.Rep, DS2.AFO.Right.Rep];
    elseif strcmp(AddLast,'Yes') == 1 % if last is 2nd condition
        Rep.Gait_Speed = [GaitSpeed.Current.Rep/TSN.TimeDist(3)*100, GaitSpeed.Last.Rep/TSN.TimeDist(3)*100];
        Rep.Cadence_ = [Cadence.Current.Rep/TSN.TimeDist(1)*100, Cadence.Last.Rep/TSN.TimeDist(1)*100];
        Rep.Stride_Length = [StrideLength.Current.Rep/TSN.TimeDist(2)*100, StrideLength.Last.Rep/TSN.TimeDist(2)*100];
        Rep.Step_Length = [StepLength.Current.Left.Rep/TSN.TimeDist(4)*100, StepLength.Current.Right.Rep/TSN.TimeDist(4)*100, StepLength.Last.Left.Rep/TSN.TimeDist(4)*100, StepLength.Last.Right.Rep/TSN.TimeDist(4)*100];
        Rep.Step_Width = [StepWidth.Current.Rep/TSN.TimeDist(6)*100, StepWidth.Last.Rep/TSN.TimeDist(6)*100];
        Rep.Stance_Phase = [Stance.Current.Left.Rep, Stance.Current.Right.Rep, Stance.Last.Left.Rep, Stance.Last.Right.Rep];
        Rep.Swing_Phase = [Swing.Current.Left.Rep, Swing.Current.Right.Rep, Swing.Last.Left.Rep, Swing.Last.Right.Rep];
        Rep.Initial_Double_Support = [DS1.Current.Left.Rep, DS1.Current.Right.Rep, DS1.Last.Left.Rep, DS1.Last.Right.Rep];
        Rep.Single_Limb_Support = [SLS.Current.Left.Rep, SLS.Current.Right.Rep, SLS.Last.Left.Rep, SLS.Last.Right.Rep];
        Rep.Secondary_Double_Support = [DS2.Current.Left.Rep, DS2.Current.Right.Rep, DS2.Last.Left.Rep, DS2.Last.Right.Rep];
    elseif strcmp(AddExtra,'Yes')  == 1 % if extra is 2nd condition
        Rep.Gait_Speed = [GaitSpeed.Current.Rep/TSN.TimeDist(3)*100, GaitSpeed.Extra.Rep/TSN.TimeDist(3)*100];
        Rep.Cadence_ = [Cadence.Current.Rep/TSN.TimeDist(1)*100, Cadence.Extra.Rep/TSN.TimeDist(1)*100];
        Rep.Stride_Length = [StrideLength.Current.Rep/TSN.TimeDist(2)*100, StrideLength.Extra.Rep/TSN.TimeDist(2)*100];
        Rep.Step_Length = [StepLength.Current.Left.Rep/TSN.TimeDist(4)*100, StepLength.Current.Right.Rep/TSN.TimeDist(4)*100, StepLength.Extra.Left.Rep/TSN.TimeDist(4)*100, StepLength.Extra.Right.Rep/TSN.TimeDist(4)*100];
        Rep.Step_Width = [StepWidth.Current.Rep/TSN.TimeDist(6)*100, StepWidth.Extra.Rep/TSN.TimeDist(6)*100];
        Rep.Stance_Phase = [Stance.Current.Left.Rep, Stance.Current.Right.Rep, Stance.Extra.Left.Rep, Stance.Extra.Right.Rep];
        Rep.Swing_Phase = [Swing.Current.Left.Rep, Swing.Current.Right.Rep, Swing.Extra.Left.Rep, Swing.Extra.Right.Rep];
        Rep.Initial_Double_Support = [DS1.Current.Left.Rep, DS1.Current.Right.Rep, DS1.Extra.Left.Rep, DS1.Extra.Right.Rep];
        Rep.Single_Limb_Support = [SLS.Current.Left.Rep, SLS.Current.Right.Rep, SLS.Extra.Left.Rep, SLS.Extra.Right.Rep];
        Rep.Secondary_Double_Support = [DS2.Current.Left.Rep, DS2.Current.Right.Rep, DS2.Extra.Left.Rep, DS2.Extra.Right.Rep];
    end
elseif NumConds == 3 % if 3 conditions
    if strcmp(AddAFO,'Yes') + strcmp(AddLast,'Yes') == 2 % if AFO and Last are 2nd and 3rd conditions
        Rep.Gait_Speed = [GaitSpeed.Current.Rep/TSN.TimeDist(3)*100, GaitSpeed.AFO.Rep/TSN.TimeDist(3)*100,GaitSpeed.Last.Rep/TSN.TimeDist(3)*100];
        Rep.Cadence_ = [Cadence.Current.Rep/TSN.TimeDist(1)*100, Cadence.AFO.Rep/TSN.TimeDist(1)*100, Cadence.Last.Rep/TSN.TimeDist(1)*100];
        Rep.Stride_Length = [StrideLength.Current.Rep/TSN.TimeDist(2)*100, StrideLength.AFO.Rep/TSN.TimeDist(2)*100, StrideLength.Last.Rep/TSN.TimeDist(2)*100];
        Rep.Step_Length = [StepLength.Current.Left.Rep/TSN.TimeDist(4)*100, StepLength.Current.Right.Rep/TSN.TimeDist(4)*100, StepLength.AFO.Left.Rep/TSN.TimeDist(4)*100, StepLength.AFO.Right.Rep/TSN.TimeDist(4)*100,...
            StepLength.Last.Left.Rep/TSN.TimeDist(4)*100, StepLength.Last.Right.Rep/TSN.TimeDist(4)*100];
        Rep.Step_Width = [StepWidth.Current.Rep/TSN.TimeDist(6)*100, StepWidth.AFO.Rep/TSN.TimeDist(6)*100, StepWidth.Last.Rep/TSN.TimeDist(6)*100];
        Rep.Stance_Phase = [Stance.Current.Left.Rep, Stance.Current.Right.Rep, Stance.AFO.Left.Rep, Stance.AFO.Right.Rep,...
            Stance.Last.Left.Rep, Stance.Last.Right.Rep];
        Rep.Swing_Phase = [Swing.Current.Left.Rep, Swing.Current.Right.Rep, Swing.AFO.Left.Rep, Swing.AFO.Right.Rep...
            Swing.Last.Left.Rep, Swing.Last.Right.Rep];
        Rep.Initial_Double_Support = [DS1.Current.Left.Rep, DS1.Current.Right.Rep, DS1.AFO.Left.Rep, DS1.AFO.Right.Rep,...
            DS1.Last.Left.Rep, DS1.Last.Right.Rep];
        Rep.Single_Limb_Support = [SLS.Current.Left.Rep, SLS.Current.Right.Rep, SLS.AFO.Left.Rep, SLS.AFO.Right.Rep,...
            SLS.Last.Left.Rep, SLS.Last.Right.Rep];
        Rep.Secondary_Double_Support = [DS2.Current.Left.Rep, DS2.Current.Right.Rep, DS2.AFO.Left.Rep, DS2.AFO.Right.Rep,...
            DS2.Last.Left.Rep, DS2.Last.Right.Rep];
    elseif strcmp(AddAFO,'Yes') + strcmp(AddExtra,'Yes') == 2 % if AFO and Extra are 2nd and 3rd conditions
        Rep.Gait_Speed = [GaitSpeed.Current.Rep/TSN.TimeDist(3)*100, GaitSpeed.AFO.Rep/TSN.TimeDist(3)*100,GaitSpeed.Extra.Rep/TSN.TimeDist(3)*100];
        Rep.Cadence_ = [Cadence.Current.Rep/TSN.TimeDist(1)*100, Cadence.AFO.Rep/TSN.TimeDist(1)*100, Cadence.Extra.Rep/TSN.TimeDist(1)*100];
        Rep.Stride_Length = [StrideLength.Current.Rep/TSN.TimeDist(2)*100, StrideLength.AFO.Rep/TSN.TimeDist(2)*100, StrideLength.Extra.Rep/TSN.TimeDist(2)*100];
        Rep.Step_Length = [StepLength.Current.Left.Rep/TSN.TimeDist(4)*100, StepLength.Current.Right.Rep/TSN.TimeDist(4)*100, StepLength.AFO.Left.Rep/TSN.TimeDist(4)*100, StepLength.AFO.Right.Rep/TSN.TimeDist(4)*100,...
            StepLength.Extra.Left.Rep/TSN.TimeDist(4)*100, StepLength.Extra.Right.Rep/TSN.TimeDist(4)*100];
        Rep.Step_Width = [StepWidth.Current.Rep/TSN.TimeDist(6)*100, StepWidth.AFO.Rep/TSN.TimeDist(6)*100, StepWidth.Extra.Rep/TSN.TimeDist(6)*100];
        Rep.Stance_Phase = [Stance.Current.Left.Rep, Stance.Current.Right.Rep, Stance.AFO.Left.Rep, Stance.AFO.Right.Rep,...
            Stance.Extra.Left.Rep, Stance.Extra.Right.Rep];
        Rep.Swing_Phase = [Swing.Current.Left.Rep, Swing.Current.Right.Rep, Swing.AFO.Left.Rep, Swing.AFO.Right.Rep...
            Swing.Extra.Left.Rep, Swing.Extra.Right.Rep];
        Rep.Initial_Double_Support = [DS1.Current.Left.Rep, DS1.Current.Right.Rep, DS1.AFO.Left.Rep, DS1.AFO.Right.Rep,...
            DS1.Extra.Left.Rep, DS1.Extra.Right.Rep];
        Rep.Single_Limb_Support = [SLS.Current.Left.Rep, SLS.Current.Right.Rep, SLS.AFO.Left.Rep, SLS.AFO.Right.Rep,...
            SLS.Extra.Left.Rep, SLS.Extra.Right.Rep];
        Rep.Secondary_Double_Support = [DS2.Current.Left.Rep, DS2.Current.Right.Rep, DS2.AFO.Left.Rep, DS2.AFO.Right.Rep,...
            DS2.Extra.Left.Rep, DS2.Extra.Right.Rep];
    elseif strcmp(AddLast,'Yes') + strcmp(AddExtra,'Yes') == 2 % if Extra and Last are 2nd and 3rd conditions
        Rep.Gait_Speed = [GaitSpeed.Current.Rep/TSN.TimeDist(3)*100,GaitSpeed.Last.Rep/TSN.TimeDist(3)*100,GaitSpeed.Extra.Rep/TSN.TimeDist(3)*100];
        Rep.Cadence_ = [Cadence.Current.Rep/TSN.TimeDist(1)*100, Cadence.Last.Rep/TSN.TimeDist(1)*100, Cadence.Extra.Rep/TSN.TimeDist(1)*100];
        Rep.Stride_Length = [StrideLength.Current.Rep/TSN.TimeDist(2)*100, StrideLength.Last.Rep/TSN.TimeDist(2)*100, StrideLength.Extra.Rep/TSN.TimeDist(2)*100];
        Rep.Step_Length = [StepLength.Current.Left.Rep/TSN.TimeDist(4)*100, StepLength.Current.Right.Rep/TSN.TimeDist(4)*100,...
            StepLength.Last.Left.Rep/TSN.TimeDist(4)*100, StepLength.Last.Right.Rep/TSN.TimeDist(4)*100,StepLength.Extra.Left.Rep/TSN.TimeDist(4)*100, StepLength.Extra.Right.Rep/TSN.TimeDist(4)*100];
        Rep.Step_Width = [StepWidth.Current.Rep/TSN.TimeDist(6)*100, StepWidth.Last.Rep/TSN.TimeDist(6)*100, StepWidth.Extra.Rep/TSN.TimeDist(6)*100];
        Rep.Stance_Phase = [Stance.Current.Left.Rep, Stance.Current.Right.Rep,...
            Stance.Last.Left.Rep, Stance.Last.Right.Rep, Stance.Extra.Left.Rep, Stance.Extra.Right.Rep];
        Rep.Swing_Phase = [Swing.Current.Left.Rep, Swing.Current.Right.Rep,...
            Swing.Last.Left.Rep, Swing.Last.Right.Rep, Swing.Extra.Left.Rep, Swing.Extra.Right.Rep];
        Rep.Initial_Double_Support = [DS1.Current.Left.Rep, DS1.Current.Right.Rep,...
            DS1.Last.Left.Rep, DS1.Last.Right.Rep, DS1.Extra.Left.Rep, DS1.Extra.Right.Rep];
        Rep.Single_Limb_Support = [SLS.Current.Left.Rep, SLS.Current.Right.Rep,...
            SLS.Last.Left.Rep, SLS.Last.Right.Rep, SLS.Extra.Left.Rep, SLS.Extra.Right.Rep];
        Rep.Secondary_Double_Support = [DS2.Current.Left.Rep, DS2.Current.Right.Rep,...
            DS2.Last.Left.Rep, DS2.Last.Right.Rep, DS2.Extra.Left.Rep, DS2.Extra.Right.Rep];
    end
elseif NumConds == 4
    Rep.Gait_Speed = [GaitSpeed.Current.Rep/TSN.TimeDist(3)*100, GaitSpeed.AFO.Rep/TSN.TimeDist(3)*100,GaitSpeed.Last.Rep/TSN.TimeDist(3)*100,GaitSpeed.Extra.Rep/TSN.TimeDist(3)*100];
    Rep.Cadence_ = [Cadence.Current.Rep/TSN.TimeDist(1)*100, Cadence.AFO.Rep/TSN.TimeDist(1)*100, Cadence.Last.Rep/TSN.TimeDist(1)*100, Cadence.Extra.Rep/TSN.TimeDist(1)*100];
    Rep.Stride_Length = [StrideLength.Current.Rep/TSN.TimeDist(2)*100, StrideLength.AFO.Rep/TSN.TimeDist(2)*100, StrideLength.Last.Rep/TSN.TimeDist(2)*100, StrideLength.Extra.Rep/TSN.TimeDist(2)*100];
    Rep.Step_Length = [StepLength.Current.Left.Rep/TSN.TimeDist(4)*100, StepLength.Current.Right.Rep/TSN.TimeDist(4)*100, StepLength.AFO.Left.Rep/TSN.TimeDist(4)*100, StepLength.AFO.Right.Rep/TSN.TimeDist(4)*100,...
        StepLength.Last.Left.Rep/TSN.TimeDist(4)*100, StepLength.Last.Right.Rep/TSN.TimeDist(4)*100,StepLength.Extra.Left.Rep/TSN.TimeDist(4)*100, StepLength.Extra.Right.Rep/TSN.TimeDist(4)*100];
    Rep.Step_Width = [StepWidth.Current.Rep/TSN.TimeDist(6)*100, StepWidth.AFO.Rep/TSN.TimeDist(6)*100, StepWidth.Last.Rep/TSN.TimeDist(6)*100, StepWidth.Extra.Rep/TSN.TimeDist(6)*100];
    Rep.Stance_Phase = [Stance.Current.Left.Rep, Stance.Current.Right.Rep, Stance.AFO.Left.Rep, Stance.AFO.Right.Rep,...
        Stance.Last.Left.Rep, Stance.Last.Right.Rep, Stance.Extra.Left.Rep, Stance.Extra.Right.Rep];
    Rep.Swing_Phase = [Swing.Current.Left.Rep, Swing.Current.Right.Rep, Swing.AFO.Left.Rep, Swing.AFO.Right.Rep...
        Swing.Last.Left.Rep, Swing.Last.Right.Rep, Swing.Extra.Left.Rep, Swing.Extra.Right.Rep];
    Rep.Initial_Double_Support = [DS1.Current.Left.Rep, DS1.Current.Right.Rep, DS1.AFO.Left.Rep, DS1.AFO.Right.Rep,...
        DS1.Last.Left.Rep, DS1.Last.Right.Rep, DS1.Extra.Left.Rep, DS1.Extra.Right.Rep];
    Rep.Single_Limb_Support = [SLS.Current.Left.Rep, SLS.Current.Right.Rep, SLS.AFO.Left.Rep, SLS.AFO.Right.Rep,...
        SLS.Last.Left.Rep, SLS.Last.Right.Rep, SLS.Extra.Left.Rep, SLS.Extra.Right.Rep];
    Rep.Secondary_Double_Support = [DS2.Current.Left.Rep, DS2.Current.Right.Rep, DS2.AFO.Left.Rep, DS2.AFO.Right.Rep,...
        DS2.Last.Left.Rep, DS2.Last.Right.Rep, DS2.Extra.Left.Rep, DS2.Extra.Right.Rep];
end

DotSize = 16;
% Plot the dots
if NumConds == 1
    X1 = 1;
    X2 = [0.9 1.1];
    % plots
    subplot(251); % gait speed
    plot(X1,Rep.Gait_Speed, '.k', 'MarkerSize',DotSize);
    subplot(252); % cadence
    plot(X1,Rep.Cadence_, '.k', 'MarkerSize',DotSize);
    subplot(253); % Stride Length
    plot(X1,Rep.Stride_Length, '.k', 'MarkerSize',DotSize);
    subplot(254); % step length
    plot(X2,Rep.Step_Length, '.k',  'MarkerSize',DotSize);
    %         plot(1.1,Rep.Step_Length(2), '.k',  'MarkerSize',DotSize);
    subplot(255); % step width
    plot(X1,Rep.Step_Width, '.k','MarkerSize',DotSize);
    subplot(256); % Stance Phase
    plot(X2,Rep.Stance_Phase, '.k',  'MarkerSize',DotSize);
    %         plot(1.1,Rep.Stance_Phase(2), '.k',  'MarkerSize',DotSize);
    subplot(257); % Swing Phase
    plot(X2,Rep.Swing_Phase, '.k',  'MarkerSize',DotSize);
    %         plot(1.1,Rep.Swing_Phase(2), '.k',  'MarkerSize',DotSize);
    subplot(258); % initial double support
    plot(X2,Rep.Initial_Double_Support, '.k', 'MarkerSize',DotSize);
    %         plot(1.1,Rep.Initial_Double_Support(2), '.k', 'MarkerSize',DotSize);
    subplot(259); %  single limb support
    plot(X2,Rep.Single_Limb_Support, '.k', 'MarkerSize',DotSize);
    %         plot(1.1,Rep.Single_Limb_Support(2), '.k', 'MarkerSize',DotSize);
    subplot(2,5,10); %  secondary double support
    plot(X2,Rep.Secondary_Double_Support, '.k', 'MarkerSize',DotSize);
    %         plot(1.1,Rep.Secondary_Double_Support(2), '.k', 'MarkerSize',DotSize);
    
elseif NumConds == 2
    % define X locations
    X1 =  [0.85 1.15];
    X2 = [0.75 0.9 1.1 1.25];
    % plots
    subplot(251); % gait speed
    plot(X1,Rep.Gait_Speed, '.k', 'MarkerSize',DotSize);
    subplot(252); % cadence
    plot(X1,Rep.Cadence_, '.k', 'MarkerSize',DotSize);
    subplot(253); % Stride Length
    plot(X1,Rep.Stride_Length, '.k', 'MarkerSize',DotSize);
    subplot(254); % step length
    plot(X2,Rep.Step_Length, '.k',  'MarkerSize',DotSize);
    subplot(255); % step width
    plot(X1,Rep.Step_Width, '.k','MarkerSize',DotSize);
    subplot(256); % Stance Phase
    plot(X2,Rep.Stance_Phase, '.k',  'MarkerSize',DotSize);
    subplot(257); % Swing Phase
    plot(X2,Rep.Swing_Phase, '.k',  'MarkerSize',DotSize);
    subplot(258); % initial double support
    plot(X2,Rep.Initial_Double_Support, '.k', 'MarkerSize',DotSize);
    subplot(259); %  single limb support
    plot(X2,Rep.Single_Limb_Support, '.k', 'MarkerSize',DotSize);
    subplot(2,5,10); %  secondary double support
    plot(X2,Rep.Secondary_Double_Support, '.k', 'MarkerSize',DotSize);
    
elseif NumConds == 3
    % define X locations
    X1 =  [0.75 1 1.25];
    X2 = [0.6 0.725 0.9375 1.0625 1.275 1.4];
    % plots
    subplot(251); % gait speed
    plot(X1,Rep.Gait_Speed, '.k', 'MarkerSize',DotSize);
    subplot(252); % cadence
    plot(X1,Rep.Cadence_, '.k', 'MarkerSize',DotSize);
    subplot(253); % Stride Length
    plot(X1,Rep.Stride_Length, '.k', 'MarkerSize',DotSize);
    subplot(254); % step length
    plot(X2,Rep.Step_Length, '.k',  'MarkerSize',DotSize);
    subplot(255); % step width
    plot(X1,Rep.Step_Width, '.k','MarkerSize',DotSize);
    subplot(256); % Stance Phase
    plot(X2,Rep.Stance_Phase, '.k',  'MarkerSize',DotSize);
    subplot(257); % Swing Phase
    plot(X2,Rep.Swing_Phase, '.k',  'MarkerSize',DotSize);
    subplot(258); % initial double support
    plot(X2,Rep.Initial_Double_Support, '.k', 'MarkerSize',DotSize);
    subplot(259); %  single limb support
    plot(X2,Rep.Single_Limb_Support, '.k', 'MarkerSize',DotSize);
    subplot(2,5,10); %  secondary double support
    plot(X2,Rep.Secondary_Double_Support, '.k', 'MarkerSize',DotSize);
    
elseif NumConds == 4
    % define X locations
    X1 = [0.7 0.9 1.1 1.3];
    X2 = [0.49 0.61 0.79 0.91 1.09 1.21 1.39 1.51];
    % plots
    subplot(251); % gait speed
    plot(X1,Rep.Gait_Speed, '.k', 'MarkerSize',DotSize);
    subplot(252); % cadence
    plot(X1,Rep.Cadence_, '.k', 'MarkerSize',DotSize);
    subplot(253); % Stride Length
    plot(X1,Rep.Stride_Length, '.k', 'MarkerSize',DotSize);
    subplot(254); % step length
    plot(X2,Rep.Step_Length, '.k',  'MarkerSize',DotSize);
    subplot(255); % step width
    plot(X1,Rep.Step_Width, '.k','MarkerSize',DotSize);
    subplot(256); % Stance Phase
    plot(X2,Rep.Stance_Phase, '.k',  'MarkerSize',DotSize);
    subplot(257); % Swing Phase
    plot(X2,Rep.Swing_Phase, '.k',  'MarkerSize',DotSize);
    subplot(258); % initial double support
    plot(X2,Rep.Initial_Double_Support, '.k', 'MarkerSize',DotSize);
    subplot(259); %  single limb support
    plot(X2,Rep.Single_Limb_Support, '.k', 'MarkerSize',DotSize);
    subplot(2,5,10); %  secondary double support
    plot(X2,Rep.Secondary_Double_Support, '.k', 'MarkerSize',DotSize);
    
end


%% Set Graph Properties and labels
% gait speed
subplot(251);
title('Gait Speed', 'FontSize',8);
ylabel('% of Normal', 'FontSize',8);
set(gca, 'FontSize',7);
% cadence
subplot(252);
title('Cadence', 'FontSize',8);
set(gca, 'FontSize',7);
% Stride Length
subplot(253);
title('Stride Length', 'FontSize',8);
set(gca, 'FontSize',7);
% step length
subplot(254);
title('Step Length', 'FontSize',8);
set(gca, 'FontSize',7);
% step width
subplot(255);
title('Step Width', 'FontSize',8);
set(gca, 'FontSize',7);
% Stance Phase
subplot(256);
title('Stance', 'FontSize',8);
ylabel('% of Gait Cycle', 'FontSize',8);
set(gca, 'FontSize',7);
% Swing Phase
subplot(257);
title('Swing', 'FontSize',8);
set(gca, 'FontSize',7);
% initial double support
subplot(258);
title('Initial Double Support', 'FontSize',8);
set(gca, 'FontSize',7);
%  single limb support
subplot(259);
title('Single Limb Support', 'FontSize',8);
set(gca, 'FontSize',7);
%  secondary double support
subplot(2,5,10);
title('Final Double Support', 'FontSize',8);
set(gca, 'FontSize',7);

%% Create legend and bar labels based on # of bars present and save figure
if NumConds == 2
    for i = [1 2 3 5] % blue bars
        subplot(2,5,i);
        ax = gca;
        ax.XTick =  [0.85 1.15];
        ax.XTickLabel = XAxis;
        ax.XTickLabelRotation = 45;
    end
    for i = [4 6 7 8 9 10] % red and green bars
        subplot(2, 5, i);
        ax = gca;
        ax.XTick =  [0.825 1.175];
        ax.XTickLabel = XAxis;
        ax.XTickLabelRotation = 45;
    end
elseif NumConds == 3 % if 3 conditions
    for i = [1 2 3 5] % blue bars
        subplot(2,5,i);
        ax = gca;
        ax.XTick =  [0.75 1 1.25];
        ax.XTickLabel = XAxis;
        ax.XTickLabelRotation = 45;
    end
    for i = [4 6 7 8 9 10] % red and green bars
        subplot(2, 5, i);
        ax = gca;
        ax.XTick =  [0.6625 1 1.3375];
        ax.XTickLabel = XAxis;
        ax.XTickLabelRotation = 45;
    end
elseif NumConds == 4 % if 4 conditions
    for i = [1 2 3 5] % blue bars
        subplot(2,5,i);
        ax = gca;
        ax.XTick =  [0.7 0.9 1.1 1.3];
        ax.XTickLabel = {char(Name.Curr), char(Name.AFO), char(Name.Last), char(Name.Extra)};
        ax.XTickLabelRotation = 45;
    end
    for i = [4 6 7 8 9 10] % red and green bars
        subplot(2, 5, i);
        ax = gca;
        ax.XTick =  [0.55 0.85 1.15 1.45];
        ax.XTickLabel = {char(Name.Curr), char(Name.AFO), char(Name.Last), char(Name.Extra)};
        ax.XTickLabelRotation = 45;
    end
end

% legend
subplot(2, 5, 7);
plot(0.5,95, '.k', 'MarkerSize',DotSize);
text(0.53,95, ' = Rep Trial','Color','k','FontSize',7);
text(0.5, 85, 'Blue = Bilateral', 'Color', Blue(1,:), 'FontSize',7);
text(0.5, 75, 'Red = Left', 'Color', RedGreen(1,:), 'FontSize',7);
text(0.5, 65, 'Green = Right', 'Color', RedGreen(2,:), 'FontSize',7);
% squeeze plots and save figure   
subplotsqueeze(GAMSplots,1.05);

saveas(GAMSplots,'GaitMeasures.png');

%% Create GDI and MAP plot
% if any of the trials contain full kinematics, plot the GDI and MAPs
if strcmp(Type.Current, 'Full Kinematics') == 1
    GDIMAP = 'Yes';
else
    if strcmp(AddAFO,'Yes') == 1
        if strcmp(Type.AFO, 'Full Kinematics') == 1
            GDIMAP = 'Yes';
        else
            GDIMAP = 'No';
        end
    elseif strcmp(AddLast,'Yes') == 1
        if strcmp(Type.Last, 'Full Kinematics') == 1
            GDIMAP = 'Yes';
        else
            GDIMAP = 'No';
        end
    elseif strcmp(AddExtra,'Yes') == 1
        if strcmp(Type.Extra, 'Full Kinematics') == 1
            GDIMAP = 'Yes';
        else
            GDIMAP = 'No';
        end
    elseif strcmp(AddAFO,'Yes') == 0 && strcmp(AddLast,'Yes') == 0 && strcmp(AddExtra,'Yes') == 0
        GDIMAP = 'No';
    end
end

if strcmp(GDIMAP, 'Yes') == 1    
    % compute GDIs
    if NumConds == 1
        RepGDIs = [GDIRep.Current(1:2)];
        GDIs = [KineData(1).TrialsGDI(:,1), KineData(1).TrialsGDI(:,2)];
    elseif NumConds == 2
        if strcmp(AddAFO,'Yes') == 1
            RepGDIs = [GDIRep.Current(1:2), GDIRep.AFO(1:2) ];
            MaxCyc = max([length(KineData(1).TrialsGDI(:,1)), length(AFOData(1).TrialsGDI(:,1))]); % find max # of cycles
            KineData(1).TrialsGDI(end+1:MaxCyc,:) = nan; % add nans to make all columns the same length
            AFOData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            GDIs = [KineData(1).TrialsGDI(:,1), KineData(1).TrialsGDI(:,2), AFOData(1).TrialsGDI(:,1), AFOData(1).TrialsGDI(:,2)];
        elseif strcmp(AddLast,'Yes') == 1
            RepGDIs = [GDIRep.Current(1:2), GDIRep.Last(1:2) ];
            MaxCyc = max([length(KineData(1).TrialsGDI(:,1)), length(LastData(1).TrialsGDI(:,1))]); % find max # of cycles
            KineData(1).TrialsGDI(end+1:MaxCyc,:) = nan; % add nans to make all columns the same length
            LastData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            GDIs = [KineData(1).TrialsGDI(:,1), KineData(1).TrialsGDI(:,2), LastData(1).TrialsGDI(:,1), LastData(1).TrialsGDI(:,2)];
        elseif strcmp(AddExtra,'Yes') == 1
            RepGDIs = [GDIRep.Current(1:2), GDIRep.Extra(1:2) ];
            MaxCyc = max([length(KineData(1).TrialsGDI(:,1)), length(ExtraData(1).TrialsGDI(:,1))]); % find max # of cycles
            KineData(1).TrialsGDI(end+1:MaxCyc,:) = nan; % add nans to make all columns the same length
            ExtraData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            GDIs = [KineData(1).TrialsGDI(:,1), KineData(1).TrialsGDI(:,2), ExtraData(1).TrialsGDI(:,1), ExtraData(1).TrialsGDI(:,2)];
        end
    elseif NumConds == 3
        if strcmp(AddExtra,'No') == 1
            RepGDIs = [GDIRep.Current(1:2) GDIRep.AFO(1:2) GDIRep.Last(1:2)];
            MaxCyc = max([length(KineData(1).TrialsGDI(:,1)), length(AFOData(1).TrialsGDI(:,1)), length(LastData(1).TrialsGDI(:,1))]); % find max # of cycles
            KineData(1).TrialsGDI(end+1:MaxCyc,:) = nan; % add nans to make all columns the same length
            AFOData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            LastData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            GDIs = [KineData(1).TrialsGDI(:,1), KineData(1).TrialsGDI(:,2), AFOData(1).TrialsGDI(:,1), AFOData(1).TrialsGDI(:,2), LastData(1).TrialsGDI(:,1), LastData(1).TrialsGDI(:,2)];
        elseif strcmp(AddLast,'No') == 1
            RepGDIs = [GDIRep.Current(1:2) GDIRep.AFO(1:2) GDIRep.Extra(1:2)];
            MaxCyc = max([length(KineData(1).TrialsGDI(:,1)), length(AFOData(1).TrialsGDI(:,1)), length(ExtraData(1).TrialsGDI(:,1))]); % find max # of cycles
            KineData(1).TrialsGDI(end+1:MaxCyc,:) = nan; % add nans to make all columns the same length
            AFOData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            ExtraData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            GDIs = [KineData(1).TrialsGDI(:,1), KineData(1).TrialsGDI(:,2), AFOData(1).TrialsGDI(:,1), AFOData(1).TrialsGDI(:,2), ExtraData(1).TrialsGDI(:,1), ExtraData(1).TrialsGDI(:,2)];
       elseif strcmp(AddAFO,'No') == 1
            RepGDIs = [GDIRep.Current(1:2) GDIRep.Last(1:2) GDIRep.Extra(1:2)];
            MaxCyc = max([length(KineData(1).TrialsGDI(:,1)), length(LastData(1).TrialsGDI(:,1)), length(ExtraData(1).TrialsGDI(:,1))]); % find max # of cycles
            KineData(1).TrialsGDI(end+1:MaxCyc,:) = nan; % add nans to make all columns the same length
            LastData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            ExtraData(1).TrialsGDI(end+1:MaxCyc,:) = nan; 
            GDIs = [KineData(1).TrialsGDI(:,1), KineData(1).TrialsGDI(:,2), LastData(1).TrialsGDI(:,1), LastData(1).TrialsGDI(:,2), ExtraData(1).TrialsGDI(:,1), ExtraData(1).TrialsGDI(:,2)];
        end
    elseif NumConds == 4
        RepGDIs = [GDIRep.Current(1:2), GDIRep.AFO(1:2), GDIRep.Last(1:2), GDIRep.Extra(1:2)];
        MaxCyc = max([length(KineData(1).TrialsGDI(:,1)), length(AFOData(1).TrialsGDI(:,1)), length(LastData(1).TrialsGDI(:,1)), length(ExtraData(1).TrialsGDI(:,1))]); % find max # of cycles
        KineData(1).TrialsGDI(end+1:MaxCyc,:) = nan; % add nans to make all columns the same length
        AFOData(1).TrialsGDI(end+1:MaxCyc,:) = nan;
        LastData(1).TrialsGDI(end+1:MaxCyc,:) = nan;
        ExtraData(1).TrialsGDI(end+1:MaxCyc,:) = nan;
        GDIs = [KineData(1).TrialsGDI(:,1), KineData(1).TrialsGDI(:,2), AFOData(1).TrialsGDI(:,1), AFOData(1).TrialsGDI(:,2),...
            LastData(1).TrialsGDI(:,1), LastData(1).TrialsGDI(:,2), ExtraData(1).TrialsGDI(:,1), ExtraData(1).TrialsGDI(:,2)];
    end
    
    % create figure
    GDIandMAPplot = figure('Position', [100, 100, 1000, 500]);
    DotSize = 16;
    subplot(3,5, [1 2 6 7 11 12]);    
    grid on; 
    
    % find mean and STD of GDIs
    meanGDIs = nanmean(GDIs,1);
    stdGDIs = nanstd(GDIs, 1);
    % Bullet graph for GDIs
    BulletGraph('V', 75, 90, meanGDIs, 0, [0 125] ,RedGreen,NumConds*2, 'Yes', 'Yes');
    % plot errorbars
    errorbar(X2, meanGDIs,stdGDIs,'.k','Marker','none');
    plot(X2, RepGDIs, 'k.','MarkerSize',DotSize);
    
    plot(0.5,122, '.k', 'MarkerSize',DotSize);
    text(0.53,122, ' = Rep Cycle','Color','k','FontSize',7);
    
    % set limits and edits axes
    ylim([0 125]);
    ylabel('GDI Score');
   % title('GDI Averages and RC');
    % Set axes
    ax = gca; 
    if NumConds == 1
        Xname = [1];
        ax.XTick = Xname;
        ax.XTickLabel = XAxis;
        ax.XTickLabelRotation = 45;
    elseif NumConds == 2
        Xname = [0.825 1.175];
        ax.XTick =  Xname; 
        ax.XTickLabel = XAxis;
        ax.XTickLabelRotation = 45;
    elseif NumConds == 3 % if 3 conditions
        Xname = [0.6625 1 1.3375];
        ax.XTick =  Xname;
        ax.XTickLabel = XAxis;
        ax.XTickLabelRotation = 45;
    elseif NumConds == 4 % if 4 conditions
        Xname = [0.55 0.85 1.15 1.45];
        ax.XTick =  Xname; 
        ax.XTickLabel = XAxis;
        ax.XTickLabelRotation = 45;
    end
    ax.FontSize = 8; 
    h = hline(100,'k');
    set(h, 'LineWidth',2); 
    
    %% Movement Analysis Profile (MAP)
    % The MAP shows how far the participant's kinematics deviate from the
    % normal ranges in terms of RMS difference. The black bars refer to 1
    % standard deviation from the inter-subject norm.
    % Deviation from the norm can be defined as the participant's average being
    % outside the inter-subject norm. Greater deviation is worse. All joints
    % are summed to the overall Gait Profile Score to the far right.
    % R. Baker et al. Gait & Posture 30 (2009) 265–269
    
    % compute control data
    ControlKine.LandRAvg = [ControlKine.Avg(:,1:3), ControlKine.Avg(:,1:3), ControlKine.Avg(:,4:6), ControlKine.Avg(:,4:6),...
        ControlKine.Avg(:,7:9), ControlKine.Avg(:,7:9), ControlKine.Avg(:,10:12), ControlKine.Avg(:,10:12)];
    ControlKine.LandRStd = [ControlKine.Std(:,1:3), ControlKine.Std(:,1:3), ControlKine.Std(:,4:6), ControlKine.Std(:,4:6),...
        ControlKine.Std(:,7:9), ControlKine.Std(:,7:9), ControlKine.Std(:,10:12), ControlKine.Std(:,10:12)];
    
    ControlKine.RmsDiffs = rms(ControlKine.Std);
    %ControlKine.StdRMS = rms(ControlKine.Std);
    x = 2:2:100; % indexing to drop from 101 points to 51
    
    % For Current Condition
    if strcmp(Type.Current, 'Full Kinematics') == 1
        MAP.Current.All = KineData(1).KineMatrix(:,:,:); % save all the kinematic cycles
         MAP.Current.All(x,:,:)  =[]; % drop to 51 points
    NumTrials = size(MAP.Current.All, 3); % find # of cycles
    for i = 1:NumTrials
        MAP.Current.Diffs(:,:,i) = abs(MAP.Current.All(:,:,i) - ControlKine.LandRAvg);
        MAP.Current.RmsDiffs(i,:) = rms(MAP.Current.Diffs(:,:,i),1); % compute RMS error for each trial
    end
    else % if it is a temp-spat study
        MAP.Last.All = zeros(51,24);  % set all kinematics to 0
        MAP.Current.RmsDiffs = zeros(51,24);
    end
    MAP.Current.AvgRMS = mean(MAP.Current.RmsDiffs); % compute overall average kineamtics
    MAP.Current.StdRMS = std(MAP.Current.RmsDiffs); % compute overall STD kinematics
    
    % for AFO Condition
    if strcmp(AddAFO,'Yes') == 1
        if strcmp(Type.AFO, 'Full Kinematics') == 1
            MAP.AFO.All = AFOData(1).KineMatrix(:,:,:); % save all the kinematic cycles
            MAP.AFO.All(x,:,:)  =[]; % drop to 51 points
            NumTrials = size(MAP.AFO.All, 3); % find # of cycles
            for i = 1:NumTrials
                MAP.AFO.Diffs(:,:,i) = abs(MAP.AFO.All(:,:,i) - ControlKine.LandRAvg);
                MAP.AFO.RmsDiffs(i,:) = rms(MAP.AFO.Diffs(:,:,i),1); % compute RMS error for each trial
            end
        else % if it is a temp-spat study
            MAP.AFO.All = zeros(51,24);  % set all kinematics to 0
             MAP.AFO.RmsDiffs =  zeros(51,24);
        end
        MAP.AFO.AvgRMS = mean(MAP.AFO.RmsDiffs); % compute overall average kineamtics
        MAP.AFO.StdRMS = std(MAP.AFO.RmsDiffs); % compute overall STD kinematics
    end
    
    % for Last Condition
    if strcmp(AddLast,'Yes') == 1
        if strcmp(Type.Last, 'Full Kinematics') == 1
            MAP.Last.All = LastData(1).KineMatrix(:,:,:); % save all the kinematic cycles
            MAP.Last.All(x,:,:)  =[]; % drop to 51 points
            NumTrials = size(MAP.Last.All, 3); % find # of cycles
            for i = 1:NumTrials
                MAP.Last.Diffs(:,:,i) = abs(MAP.Last.All(:,:,i) - ControlKine.LandRAvg);
                MAP.Last.RmsDiffs(i,:) = rms(MAP.Last.Diffs(:,:,i),1); % compute RMS error for each trial
            end
        else % if it is a temp-spat study
            MAP.Last.All = zeros(51,24);  % set all kinematics to 0
             MAP.Last.RmsDiffs =  zeros(51,24); 
        end
        MAP.Last.AvgRMS = mean(MAP.Last.RmsDiffs); % compute overall average kineamtics
        MAP.Last.StdRMS = std(MAP.Last.RmsDiffs); % compute overall STD kinematics
    end
    
    % for Extra Condition
    if strcmp(AddExtra,'Yes') == 1
        if strcmp(Type.Extra, 'Full Kinematics') == 1
            MAP.Extra.All = ExtraData(1).KineMatrix(:,:,:); % save all the kinematic cycles
            MAP.Extra.All(x,:,:)  =[]; % drop to 51 points
            NumTrials = size(MAP.Extra.All, 3); % find # of cycles
            for i = 1:NumTrials
                MAP.Extra.Diffs(:,:,i) = abs(MAP.Extra.All(:,:,i) - ControlKine.LandRAvg);
                MAP.Extra.RmsDiffs(i,:) = rms(MAP.Extra.Diffs(:,:,i),1); % compute RMS error for each trial
            end
        else % if it is a temp-spat study
            MAP.Extra.All = zeros(51,24);  % set all kinematics to 0
            MAP.Extra.RmsDiffs = zeros(51,24); % compute RMS error for each trial
        end
        MAP.Extra.AvgRMS = mean(MAP.Extra.RmsDiffs); % compute overall average kineamtics
        MAP.Extra.StdRMS = std(MAP.Extra.RmsDiffs); % compute overall STD kinematics
    end
    
    % redefine values in structure format
    if NumConds == 1
        Pelvis.Tilt.Avg = [MAP.Current.AvgRMS(1), MAP.Current.AvgRMS(4)]; % Pelvis
        Pelvis.Tilt.Std = [MAP.Current.StdRMS(1), MAP.Current.StdRMS(4)];
        Pelvis.Obl.Avg = [MAP.Current.AvgRMS(2), MAP.Current.AvgRMS(5)];
        Pelvis.Obl.Std = [MAP.Current.StdRMS(2), MAP.Current.StdRMS(5)];
        Pelvis.Rot.Avg = [MAP.Current.AvgRMS(3), MAP.Current.AvgRMS(6)];
        Pelvis.Rot.Std = [MAP.Current.StdRMS(4), MAP.Current.StdRMS(6)];
        Hip.FE.Avg = [MAP.Current.AvgRMS(7), MAP.Current.AvgRMS(10)]; % Hip
        Hip.FE.Std = [MAP.Current.StdRMS(7), MAP.Current.StdRMS(10)];
        Hip.AbAd.Avg = [MAP.Current.AvgRMS(8), MAP.Current.AvgRMS(11)];
        Hip.AbAd.Std = [MAP.Current.StdRMS(8), MAP.Current.StdRMS(11)];
        Hip.Rot.Avg = [MAP.Current.AvgRMS(9), MAP.Current.AvgRMS(12)];
        Hip.Rot.Std = [MAP.Current.StdRMS(9), MAP.Current.StdRMS(12)];
        Knee.FE.Avg = [MAP.Current.AvgRMS(13), MAP.Current.AvgRMS(16)]; % Knee
        Knee.FE.Std = [MAP.Current.StdRMS(13), MAP.Current.StdRMS(16)];
        Ankle.FE.Avg = [MAP.Current.AvgRMS(19), MAP.Current.AvgRMS(22)]; % Ankle
        Ankle.FE.Std = [MAP.Current.StdRMS(19), MAP.Current.StdRMS(22)];
        Ankle.Prog.Avg = [MAP.Current.AvgRMS(20), MAP.Current.AvgRMS(23)];
        Ankle.Prog.Std = [MAP.Current.StdRMS(20), MAP.Current.StdRMS(23)];
    elseif NumConds == 2
        if strcmp(AddAFO,'Yes') == 1 %&& strcmp(Type.AFO,'Full Kinematics') == 1
            Pelvis.Tilt.Avg = [MAP.Current.AvgRMS(1), MAP.Current.AvgRMS(4), MAP.AFO.AvgRMS(1), MAP.AFO.AvgRMS(4)]; % Pelvis
            Pelvis.Tilt.Std = [MAP.Current.StdRMS(1), MAP.Current.StdRMS(4), MAP.AFO.StdRMS(1), MAP.AFO.StdRMS(4)];
            Pelvis.Obl.Avg = [MAP.Current.AvgRMS(2), MAP.Current.AvgRMS(5), MAP.AFO.AvgRMS(2), MAP.AFO.AvgRMS(5)];
            Pelvis.Obl.Std = [MAP.Current.StdRMS(2), MAP.Current.StdRMS(5), MAP.AFO.StdRMS(2), MAP.AFO.StdRMS(5)];
            Pelvis.Rot.Avg = [MAP.Current.AvgRMS(3), MAP.Current.AvgRMS(6), MAP.AFO.AvgRMS(3), MAP.AFO.AvgRMS(6)];
            Pelvis.Rot.Std = [MAP.Current.StdRMS(4), MAP.Current.StdRMS(6), MAP.AFO.StdRMS(4), MAP.AFO.StdRMS(6)];
            Hip.FE.Avg = [MAP.Current.AvgRMS(7), MAP.Current.AvgRMS(10), MAP.AFO.AvgRMS(7), MAP.AFO.AvgRMS(10)]; % Hip
            Hip.FE.Std = [MAP.Current.StdRMS(7), MAP.Current.StdRMS(10), MAP.AFO.StdRMS(7), MAP.AFO.StdRMS(10)];
            Hip.AbAd.Avg = [MAP.Current.AvgRMS(8), MAP.Current.AvgRMS(11), MAP.AFO.AvgRMS(8), MAP.AFO.AvgRMS(11)];
            Hip.AbAd.Std = [MAP.Current.StdRMS(8), MAP.Current.StdRMS(11), MAP.AFO.StdRMS(8), MAP.AFO.StdRMS(11)];
            Hip.Rot.Avg = [MAP.Current.AvgRMS(9), MAP.Current.AvgRMS(12), MAP.AFO.AvgRMS(9), MAP.AFO.AvgRMS(12)];
            Hip.Rot.Std = [MAP.Current.StdRMS(9), MAP.Current.StdRMS(12), MAP.AFO.StdRMS(9), MAP.AFO.StdRMS(12)];
            Knee.FE.Avg = [MAP.Current.AvgRMS(13), MAP.Current.AvgRMS(16), MAP.AFO.AvgRMS(13), MAP.AFO.AvgRMS(16)]; % Knee
            Knee.FE.Std = [MAP.Current.StdRMS(13), MAP.Current.StdRMS(16), MAP.AFO.StdRMS(13), MAP.AFO.StdRMS(16)];
            Ankle.FE.Avg = [MAP.Current.AvgRMS(19), MAP.Current.AvgRMS(22), MAP.AFO.AvgRMS(19), MAP.AFO.AvgRMS(22)]; % Ankle
            Ankle.FE.Std = [MAP.Current.StdRMS(19), MAP.Current.StdRMS(22), MAP.AFO.StdRMS(19), MAP.AFO.StdRMS(22)];
            Ankle.Prog.Avg = [MAP.Current.AvgRMS(20), MAP.Current.AvgRMS(23), MAP.AFO.AvgRMS(20), MAP.AFO.AvgRMS(23)];
            Ankle.Prog.Std = [MAP.Current.StdRMS(20), MAP.Current.StdRMS(23), MAP.AFO.StdRMS(20), MAP.AFO.StdRMS(23)];
        elseif strcmp(AddLast,'Yes') == 1 % && strcmp(Type.Last,'Full Kinematics') == 1
            Pelvis.Tilt.Avg = [MAP.Current.AvgRMS(1), MAP.Current.AvgRMS(4), MAP.Last.AvgRMS(1), MAP.Last.AvgRMS(4)]; % Pelvis
            Pelvis.Tilt.Std = [MAP.Current.StdRMS(1), MAP.Current.StdRMS(4), MAP.Last.StdRMS(1), MAP.Last.StdRMS(4)];
            Pelvis.Obl.Avg = [MAP.Current.AvgRMS(2), MAP.Current.AvgRMS(5), MAP.Last.AvgRMS(2), MAP.Last.AvgRMS(5)];
            Pelvis.Obl.Std = [MAP.Current.StdRMS(2), MAP.Current.StdRMS(5), MAP.Last.StdRMS(2), MAP.Last.StdRMS(5)];
            Pelvis.Rot.Avg = [MAP.Current.AvgRMS(3), MAP.Current.AvgRMS(6), MAP.Last.AvgRMS(3), MAP.Last.AvgRMS(6)];
            Pelvis.Rot.Std = [MAP.Current.StdRMS(4), MAP.Current.StdRMS(6), MAP.Last.StdRMS(4), MAP.Last.StdRMS(6)];
            Hip.FE.Avg = [MAP.Current.AvgRMS(7), MAP.Current.AvgRMS(10), MAP.Last.AvgRMS(7), MAP.Last.AvgRMS(10)]; % Hip
            Hip.FE.Std = [MAP.Current.StdRMS(7), MAP.Current.StdRMS(10), MAP.Last.StdRMS(7), MAP.Last.StdRMS(10)];
            Hip.AbAd.Avg = [MAP.Current.AvgRMS(8), MAP.Current.AvgRMS(11), MAP.Last.AvgRMS(8), MAP.Last.AvgRMS(11)];
            Hip.AbAd.Std = [MAP.Current.StdRMS(8), MAP.Current.StdRMS(11), MAP.Last.StdRMS(8), MAP.Last.StdRMS(11)];
            Hip.Rot.Avg = [MAP.Current.AvgRMS(9), MAP.Current.AvgRMS(12), MAP.Last.AvgRMS(9), MAP.Last.AvgRMS(12)];
            Hip.Rot.Std = [MAP.Current.StdRMS(9), MAP.Current.StdRMS(12), MAP.Last.StdRMS(9), MAP.Last.StdRMS(12)];
            Knee.FE.Avg = [MAP.Current.AvgRMS(13), MAP.Current.AvgRMS(16), MAP.Last.AvgRMS(13), MAP.Last.AvgRMS(16)]; % Knee
            Knee.FE.Std = [MAP.Current.StdRMS(13), MAP.Current.StdRMS(16), MAP.Last.StdRMS(13), MAP.Last.StdRMS(16)];
            Ankle.FE.Avg = [MAP.Current.AvgRMS(19), MAP.Current.AvgRMS(22), MAP.Last.AvgRMS(19), MAP.Last.AvgRMS(22)]; % Ankle
            Ankle.FE.Std = [MAP.Current.StdRMS(19), MAP.Current.StdRMS(22), MAP.Last.StdRMS(19), MAP.Last.StdRMS(22)];
            Ankle.Prog.Avg = [MAP.Current.AvgRMS(20), MAP.Current.AvgRMS(23), MAP.Last.AvgRMS(20), MAP.Last.AvgRMS(23)];
            Ankle.Prog.Std = [MAP.Current.StdRMS(20), MAP.Current.StdRMS(23), MAP.Last.StdRMS(20), MAP.Last.StdRMS(23)];
        elseif strcmp(AddExtra,'Yes') == 1 %&& strcmp(Type.Extra,'Full Kinematics') == 1
            Pelvis.Tilt.Avg = [MAP.Current.AvgRMS(1), MAP.Current.AvgRMS(4), MAP.Extra.AvgRMS(1), MAP.Extra.AvgRMS(4)]; % Pelvis
            Pelvis.Tilt.Std = [MAP.Current.StdRMS(1), MAP.Current.StdRMS(4), MAP.Extra.StdRMS(1), MAP.Extra.StdRMS(4)];
            Pelvis.Obl.Avg = [MAP.Current.AvgRMS(2), MAP.Current.AvgRMS(5), MAP.Extra.AvgRMS(2), MAP.Extra.AvgRMS(5)];
            Pelvis.Obl.Std = [MAP.Current.StdRMS(2), MAP.Current.StdRMS(5), MAP.Extra.StdRMS(2), MAP.Extra.StdRMS(5)];
            Pelvis.Rot.Avg = [MAP.Current.AvgRMS(3), MAP.Current.AvgRMS(6), MAP.Extra.AvgRMS(3), MAP.Extra.AvgRMS(6)];
            Pelvis.Rot.Std = [MAP.Current.StdRMS(4), MAP.Current.StdRMS(6), MAP.Extra.StdRMS(4), MAP.Extra.StdRMS(6)];
            Hip.FE.Avg = [MAP.Current.AvgRMS(7), MAP.Current.AvgRMS(10), MAP.Extra.AvgRMS(7), MAP.Extra.AvgRMS(10)]; % Hip
            Hip.FE.Std = [MAP.Current.StdRMS(7), MAP.Current.StdRMS(10), MAP.Extra.StdRMS(7), MAP.Extra.StdRMS(10)];
            Hip.AbAd.Avg = [MAP.Current.AvgRMS(8), MAP.Current.AvgRMS(11), MAP.Extra.AvgRMS(8), MAP.Extra.AvgRMS(11)];
            Hip.AbAd.Std = [MAP.Current.StdRMS(8), MAP.Current.StdRMS(11), MAP.Extra.StdRMS(8), MAP.Extra.StdRMS(11)];
            Hip.Rot.Avg = [MAP.Current.AvgRMS(9), MAP.Current.AvgRMS(12), MAP.Extra.AvgRMS(9), MAP.Extra.AvgRMS(12)];
            Hip.Rot.Std = [MAP.Current.StdRMS(9), MAP.Current.StdRMS(12), MAP.Extra.StdRMS(9), MAP.Extra.StdRMS(12)];
            Knee.FE.Avg = [MAP.Current.AvgRMS(13), MAP.Current.AvgRMS(16), MAP.Extra.AvgRMS(13), MAP.Extra.AvgRMS(16)]; % Knee
            Knee.FE.Std = [MAP.Current.StdRMS(13), MAP.Current.StdRMS(16), MAP.Extra.StdRMS(13), MAP.Extra.StdRMS(16)];
            Ankle.FE.Avg = [MAP.Current.AvgRMS(19), MAP.Current.AvgRMS(22), MAP.Extra.AvgRMS(19), MAP.Extra.AvgRMS(22)]; % Ankle
            Ankle.FE.Std = [MAP.Current.StdRMS(19), MAP.Current.StdRMS(22), MAP.Extra.StdRMS(19), MAP.Extra.StdRMS(22)];
            Ankle.Prog.Avg = [MAP.Current.AvgRMS(20), MAP.Current.AvgRMS(23), MAP.Extra.AvgRMS(20), MAP.Extra.AvgRMS(23)];
            Ankle.Prog.Std = [MAP.Current.StdRMS(20), MAP.Current.StdRMS(23), MAP.Extra.StdRMS(20), MAP.Extra.StdRMS(23)];
        end
    elseif NumConds == 3
        if strcmp(AddExtra,'No') == 1 %&& strcmp(Type.AFO,'Full Kinematics') == 1 && strcmp(Type.Last,'Full Kinematics') == 1
            Pelvis.Tilt.Avg = [MAP.Current.AvgRMS(1), MAP.Current.AvgRMS(4), MAP.AFO.AvgRMS(1), MAP.AFO.AvgRMS(4), MAP.Last.AvgRMS(1), MAP.Last.AvgRMS(4)]; % Pelvis
            Pelvis.Tilt.Std = [MAP.Current.StdRMS(1), MAP.Current.StdRMS(4), MAP.AFO.StdRMS(1), MAP.AFO.StdRMS(4), MAP.Last.StdRMS(1), MAP.Last.StdRMS(4)];
            Pelvis.Obl.Avg = [MAP.Current.AvgRMS(2), MAP.Current.AvgRMS(5), MAP.AFO.AvgRMS(2), MAP.AFO.AvgRMS(5), MAP.Last.AvgRMS(2), MAP.Last.AvgRMS(5)];
            Pelvis.Obl.Std = [MAP.Current.StdRMS(2), MAP.Current.StdRMS(5), MAP.AFO.StdRMS(2), MAP.AFO.StdRMS(5), MAP.Last.StdRMS(2), MAP.Last.StdRMS(5)];
            Pelvis.Rot.Avg = [MAP.Current.AvgRMS(3), MAP.Current.AvgRMS(6), MAP.AFO.AvgRMS(3), MAP.AFO.AvgRMS(6), MAP.Last.AvgRMS(3), MAP.Last.AvgRMS(6)];
            Pelvis.Rot.Std = [MAP.Current.StdRMS(4), MAP.Current.StdRMS(6), MAP.AFO.StdRMS(4), MAP.AFO.StdRMS(6), MAP.Last.StdRMS(4), MAP.Last.StdRMS(6)];
            Hip.FE.Avg = [MAP.Current.AvgRMS(7), MAP.Current.AvgRMS(10), MAP.AFO.AvgRMS(7), MAP.AFO.AvgRMS(10), MAP.Last.AvgRMS(7), MAP.Last.AvgRMS(10)]; % Hip
            Hip.FE.Std = [MAP.Current.StdRMS(7), MAP.Current.StdRMS(10), MAP.AFO.StdRMS(7), MAP.AFO.StdRMS(10), MAP.Last.StdRMS(7), MAP.Last.StdRMS(10)];
            Hip.AbAd.Avg = [MAP.Current.AvgRMS(8), MAP.Current.AvgRMS(11), MAP.AFO.AvgRMS(8), MAP.AFO.AvgRMS(11), MAP.Last.AvgRMS(8), MAP.Last.AvgRMS(11)];
            Hip.AbAd.Std = [MAP.Current.StdRMS(8), MAP.Current.StdRMS(11), MAP.AFO.StdRMS(8), MAP.AFO.StdRMS(11), MAP.Last.StdRMS(8), MAP.Last.StdRMS(11)];
            Hip.Rot.Avg = [MAP.Current.AvgRMS(9), MAP.Current.AvgRMS(12), MAP.AFO.AvgRMS(9), MAP.AFO.AvgRMS(12), MAP.Last.AvgRMS(9), MAP.Last.AvgRMS(12)];
            Hip.Rot.Std = [MAP.Current.StdRMS(9), MAP.Current.StdRMS(12), MAP.AFO.StdRMS(9), MAP.AFO.StdRMS(12), MAP.Last.StdRMS(9), MAP.Last.StdRMS(12)];
            Knee.FE.Avg = [MAP.Current.AvgRMS(13), MAP.Current.AvgRMS(16), MAP.AFO.AvgRMS(13), MAP.AFO.AvgRMS(16), MAP.Last.AvgRMS(13), MAP.Last.AvgRMS(16)]; % Knee
            Knee.FE.Std = [MAP.Current.StdRMS(13), MAP.Current.StdRMS(16), MAP.AFO.StdRMS(13), MAP.AFO.StdRMS(16), MAP.Last.StdRMS(13), MAP.Last.StdRMS(16)];
            Ankle.FE.Avg = [MAP.Current.AvgRMS(19), MAP.Current.AvgRMS(22), MAP.AFO.AvgRMS(19), MAP.AFO.AvgRMS(22), MAP.Last.AvgRMS(19), MAP.Last.AvgRMS(22)]; % Ankle
            Ankle.FE.Std = [MAP.Current.StdRMS(19), MAP.Current.StdRMS(22), MAP.AFO.StdRMS(19), MAP.AFO.StdRMS(22), MAP.Last.StdRMS(19), MAP.Last.StdRMS(22)];
            Ankle.Prog.Avg = [MAP.Current.AvgRMS(20), MAP.Current.AvgRMS(23), MAP.AFO.AvgRMS(20), MAP.AFO.AvgRMS(23), MAP.Last.AvgRMS(20), MAP.Last.AvgRMS(23)];
            Ankle.Prog.Std = [MAP.Current.StdRMS(20), MAP.Current.StdRMS(23), MAP.AFO.StdRMS(20), MAP.AFO.StdRMS(23), MAP.Last.StdRMS(20), MAP.Last.StdRMS(23)];
        elseif strcmp(AddLast,'No') == 1 %&& strcmp(Type.Last,'Full Kinematics') == 1  && strcmp(Type.Extra,'Full Kinematics') == 1
            Pelvis.Tilt.Avg = [MAP.Current.AvgRMS(1), MAP.Current.AvgRMS(4), MAP.AFO.AvgRMS(1), MAP.AFO.AvgRMS(4), MAP.Extra.AvgRMS(1), MAP.Extra.AvgRMS(4)]; % Pelvis
            Pelvis.Tilt.Std = [MAP.Current.StdRMS(1), MAP.Current.StdRMS(4), MAP.AFO.StdRMS(1), MAP.AFO.StdRMS(4), MAP.Extra.StdRMS(1), MAP.Extra.StdRMS(4)];
            Pelvis.Obl.Avg = [MAP.Current.AvgRMS(2), MAP.Current.AvgRMS(5), MAP.AFO.AvgRMS(2), MAP.AFO.AvgRMS(5), MAP.Extra.AvgRMS(2), MAP.Extra.AvgRMS(5)];
            Pelvis.Obl.Std = [MAP.Current.StdRMS(2), MAP.Current.StdRMS(5), MAP.AFO.StdRMS(2), MAP.AFO.StdRMS(5), MAP.Extra.StdRMS(2), MAP.Extra.StdRMS(5)];
            Pelvis.Rot.Avg = [MAP.Current.AvgRMS(3), MAP.Current.AvgRMS(6), MAP.AFO.AvgRMS(3), MAP.AFO.AvgRMS(6), MAP.Extra.AvgRMS(3), MAP.Extra.AvgRMS(6)];
            Pelvis.Rot.Std = [MAP.Current.StdRMS(4), MAP.Current.StdRMS(6), MAP.AFO.StdRMS(4), MAP.AFO.StdRMS(6), MAP.Extra.StdRMS(4), MAP.Extra.StdRMS(6)];
            Hip.FE.Avg = [MAP.Current.AvgRMS(7), MAP.Current.AvgRMS(10), MAP.AFO.AvgRMS(7), MAP.AFO.AvgRMS(10), MAP.Extra.AvgRMS(7), MAP.Extra.AvgRMS(10)]; % Hip
            Hip.FE.Std = [MAP.Current.StdRMS(7), MAP.Current.StdRMS(10), MAP.AFO.StdRMS(7), MAP.AFO.StdRMS(10), MAP.Extra.StdRMS(7), MAP.Extra.StdRMS(10)];
            Hip.AbAd.Avg = [MAP.Current.AvgRMS(8), MAP.Current.AvgRMS(11), MAP.AFO.AvgRMS(8), MAP.AFO.AvgRMS(11), MAP.Extra.AvgRMS(8), MAP.Extra.AvgRMS(11)];
            Hip.AbAd.Std = [MAP.Current.StdRMS(8), MAP.Current.StdRMS(11), MAP.AFO.StdRMS(8), MAP.AFO.StdRMS(11), MAP.Extra.StdRMS(8), MAP.Extra.StdRMS(11)];
            Hip.Rot.Avg = [MAP.Current.AvgRMS(9), MAP.Current.AvgRMS(12), MAP.AFO.AvgRMS(9), MAP.AFO.AvgRMS(12), MAP.Extra.AvgRMS(9), MAP.Extra.AvgRMS(12)];
            Hip.Rot.Std = [MAP.Current.StdRMS(9), MAP.Current.StdRMS(12), MAP.AFO.StdRMS(9), MAP.AFO.StdRMS(12), MAP.Extra.StdRMS(9), MAP.Extra.StdRMS(12)];
            Knee.FE.Avg = [MAP.Current.AvgRMS(13), MAP.Current.AvgRMS(16), MAP.AFO.AvgRMS(13), MAP.AFO.AvgRMS(16), MAP.Extra.AvgRMS(13), MAP.Extra.AvgRMS(16)]; % Knee
            Knee.FE.Std = [MAP.Current.StdRMS(13), MAP.Current.StdRMS(16), MAP.AFO.StdRMS(13), MAP.AFO.StdRMS(16), MAP.Extra.StdRMS(13), MAP.Extra.StdRMS(16)];
            Ankle.FE.Avg = [MAP.Current.AvgRMS(19), MAP.Current.AvgRMS(22), MAP.AFO.AvgRMS(19), MAP.AFO.AvgRMS(22), MAP.Extra.AvgRMS(19), MAP.Extra.AvgRMS(22)]; % Ankle
            Ankle.FE.Std = [MAP.Current.StdRMS(19), MAP.Current.StdRMS(22), MAP.AFO.StdRMS(19), MAP.AFO.StdRMS(22), MAP.Extra.StdRMS(19), MAP.Extra.StdRMS(22)];
            Ankle.Prog.Avg = [MAP.Current.AvgRMS(20), MAP.Current.AvgRMS(23), MAP.AFO.AvgRMS(20), MAP.AFO.AvgRMS(23), MAP.Extra.AvgRMS(20), MAP.Extra.AvgRMS(23)];
            Ankle.Prog.Std = [MAP.Current.StdRMS(20), MAP.Current.StdRMS(23), MAP.AFO.StdRMS(20), MAP.AFO.StdRMS(23), MAP.Extra.StdRMS(20), MAP.Extra.StdRMS(23)];
        elseif strcmp(AddAFO,'No') == 1 % && strcmp(Type.Extra,'Full Kinematics') == 1
            Pelvis.Tilt.Avg = [MAP.Current.AvgRMS(1), MAP.Current.AvgRMS(4), MAP.Last.AvgRMS(1), MAP.Last.AvgRMS(4), MAP.Extra.AvgRMS(1), MAP.Extra.AvgRMS(4)]; % Pelvis
            Pelvis.Tilt.Std = [MAP.Current.StdRMS(1), MAP.Current.StdRMS(4), MAP.Last.StdRMS(1), MAP.Last.StdRMS(4), MAP.Extra.StdRMS(1), MAP.Extra.StdRMS(4)];
            Pelvis.Obl.Avg = [MAP.Current.AvgRMS(2), MAP.Current.AvgRMS(5), MAP.Last.AvgRMS(2), MAP.Last.AvgRMS(5), MAP.Extra.AvgRMS(2), MAP.Extra.AvgRMS(5)];
            Pelvis.Obl.Std = [MAP.Current.StdRMS(2), MAP.Current.StdRMS(5), MAP.Last.StdRMS(2), MAP.Last.StdRMS(5), MAP.Extra.StdRMS(2), MAP.Extra.StdRMS(5)];
            Pelvis.Rot.Avg = [MAP.Current.AvgRMS(3), MAP.Current.AvgRMS(6), MAP.Last.AvgRMS(3), MAP.Last.AvgRMS(6), MAP.Extra.AvgRMS(3), MAP.Extra.AvgRMS(6)];
            Pelvis.Rot.Std = [MAP.Current.StdRMS(4), MAP.Current.StdRMS(6), MAP.Last.StdRMS(4), MAP.Last.StdRMS(6), MAP.Extra.StdRMS(4), MAP.Extra.StdRMS(6)];
            Hip.FE.Avg = [MAP.Current.AvgRMS(7), MAP.Current.AvgRMS(10), MAP.Last.AvgRMS(7), MAP.Last.AvgRMS(10), MAP.Extra.AvgRMS(7), MAP.Extra.AvgRMS(10)]; % Hip
            Hip.FE.Std = [MAP.Current.StdRMS(7), MAP.Current.StdRMS(10), MAP.Last.StdRMS(7), MAP.Last.StdRMS(10), MAP.Extra.StdRMS(7), MAP.Extra.StdRMS(10)];
            Hip.AbAd.Avg = [MAP.Current.AvgRMS(8), MAP.Current.AvgRMS(11), MAP.Last.AvgRMS(8), MAP.Last.AvgRMS(11), MAP.Extra.AvgRMS(8), MAP.Extra.AvgRMS(11)];
            Hip.AbAd.Std = [MAP.Current.StdRMS(8), MAP.Current.StdRMS(11), MAP.Last.StdRMS(8), MAP.Last.StdRMS(11), MAP.Extra.StdRMS(8), MAP.Extra.StdRMS(11)];
            Hip.Rot.Avg = [MAP.Current.AvgRMS(9), MAP.Current.AvgRMS(12), MAP.Last.AvgRMS(9), MAP.Last.AvgRMS(12), MAP.Extra.AvgRMS(9), MAP.Extra.AvgRMS(12)];
            Hip.Rot.Std = [MAP.Current.StdRMS(9), MAP.Current.StdRMS(12), MAP.Last.StdRMS(9), MAP.Last.StdRMS(12), MAP.Extra.StdRMS(9), MAP.Extra.StdRMS(12)];
            Knee.FE.Avg = [MAP.Current.AvgRMS(13), MAP.Current.AvgRMS(16), MAP.Last.AvgRMS(13), MAP.Last.AvgRMS(16), MAP.Extra.AvgRMS(13), MAP.Extra.AvgRMS(16)]; % Knee
            Knee.FE.Std = [MAP.Current.StdRMS(13), MAP.Current.StdRMS(16), MAP.Last.StdRMS(13), MAP.Last.StdRMS(16), MAP.Extra.StdRMS(13), MAP.Extra.StdRMS(16)];
            Ankle.FE.Avg = [MAP.Current.AvgRMS(19), MAP.Current.AvgRMS(22), MAP.Last.AvgRMS(19), MAP.Last.AvgRMS(22), MAP.Extra.AvgRMS(19), MAP.Extra.AvgRMS(22)]; % Ankle
            Ankle.FE.Std = [MAP.Current.StdRMS(19), MAP.Current.StdRMS(22), MAP.Last.StdRMS(19), MAP.Last.StdRMS(22), MAP.Extra.StdRMS(19), MAP.Extra.StdRMS(22)];
            Ankle.Prog.Avg = [MAP.Current.AvgRMS(20), MAP.Current.AvgRMS(23), MAP.Last.AvgRMS(20), MAP.Last.AvgRMS(23), MAP.Extra.AvgRMS(20), MAP.Extra.AvgRMS(23)];
            Ankle.Prog.Std = [MAP.Current.StdRMS(20), MAP.Current.StdRMS(23), MAP.Last.StdRMS(20), MAP.Last.StdRMS(23), MAP.Extra.StdRMS(20), MAP.Extra.StdRMS(23)];
        end
    elseif NumConds == 4
        Pelvis.Tilt.Avg = [MAP.Current.AvgRMS(1), MAP.Current.AvgRMS(4), MAP.AFO.AvgRMS(1), MAP.AFO.AvgRMS(4), MAP.Last.AvgRMS(1), MAP.Last.AvgRMS(4), MAP.Extra.AvgRMS(1), MAP.Extra.AvgRMS(4)]; % Pelvis
        Pelvis.Tilt.Std = [MAP.Current.StdRMS(1), MAP.Current.StdRMS(4), MAP.AFO.StdRMS(1), MAP.AFO.StdRMS(4), MAP.Last.StdRMS(1), MAP.Last.StdRMS(4), MAP.Extra.StdRMS(1), MAP.Extra.StdRMS(4)];
        Pelvis.Obl.Avg = [MAP.Current.AvgRMS(2), MAP.Current.AvgRMS(5), MAP.AFO.AvgRMS(2), MAP.AFO.AvgRMS(5), MAP.Last.AvgRMS(2), MAP.Last.AvgRMS(5), MAP.Extra.AvgRMS(2), MAP.Extra.AvgRMS(5)];
        Pelvis.Obl.Std = [MAP.Current.StdRMS(2), MAP.Current.StdRMS(5), MAP.AFO.StdRMS(2), MAP.AFO.StdRMS(5), MAP.Last.StdRMS(2), MAP.Last.StdRMS(5), MAP.Extra.StdRMS(2), MAP.Extra.StdRMS(5)];
        Pelvis.Rot.Avg = [MAP.Current.AvgRMS(3), MAP.Current.AvgRMS(6), MAP.AFO.AvgRMS(3), MAP.AFO.AvgRMS(6), MAP.Last.AvgRMS(3), MAP.Last.AvgRMS(6), MAP.Extra.AvgRMS(3), MAP.Extra.AvgRMS(6)];
        Pelvis.Rot.Std = [MAP.Current.StdRMS(4), MAP.Current.StdRMS(6), MAP.AFO.StdRMS(4), MAP.AFO.StdRMS(6), MAP.Last.StdRMS(4), MAP.Last.StdRMS(6), MAP.Extra.StdRMS(4), MAP.Extra.StdRMS(6)];
        Hip.FE.Avg = [MAP.Current.AvgRMS(7), MAP.Current.AvgRMS(10), MAP.AFO.AvgRMS(7), MAP.AFO.AvgRMS(10), MAP.Last.AvgRMS(7), MAP.Last.AvgRMS(10), MAP.Extra.AvgRMS(7), MAP.Extra.AvgRMS(10)]; % Hip
        Hip.FE.Std = [MAP.Current.StdRMS(7), MAP.Current.StdRMS(10), MAP.AFO.StdRMS(7), MAP.AFO.StdRMS(10), MAP.Last.StdRMS(7), MAP.Last.StdRMS(10), MAP.Extra.StdRMS(7), MAP.Extra.StdRMS(10)];
        Hip.AbAd.Avg = [MAP.Current.AvgRMS(8), MAP.Current.AvgRMS(11), MAP.AFO.AvgRMS(8), MAP.AFO.AvgRMS(11), MAP.Last.AvgRMS(8), MAP.Last.AvgRMS(11), MAP.Extra.AvgRMS(8), MAP.Extra.AvgRMS(11)];
        Hip.AbAd.Std = [MAP.Current.StdRMS(8), MAP.Current.StdRMS(11), MAP.AFO.StdRMS(8), MAP.AFO.StdRMS(11), MAP.Last.StdRMS(8), MAP.Last.StdRMS(11), MAP.Extra.StdRMS(8), MAP.Extra.StdRMS(11)];
        Hip.Rot.Avg = [MAP.Current.AvgRMS(9), MAP.Current.AvgRMS(12), MAP.AFO.AvgRMS(9), MAP.AFO.AvgRMS(12), MAP.Last.AvgRMS(9), MAP.Last.AvgRMS(12), MAP.Extra.AvgRMS(9), MAP.Extra.AvgRMS(12)];
        Hip.Rot.Std = [MAP.Current.StdRMS(9), MAP.Current.StdRMS(12), MAP.AFO.StdRMS(9), MAP.AFO.StdRMS(12), MAP.Last.StdRMS(9), MAP.Last.StdRMS(12), MAP.Extra.StdRMS(9), MAP.Extra.StdRMS(12)];
        Knee.FE.Avg = [MAP.Current.AvgRMS(13), MAP.Current.AvgRMS(16), MAP.AFO.AvgRMS(13), MAP.AFO.AvgRMS(16), MAP.Last.AvgRMS(13), MAP.Last.AvgRMS(16), MAP.Extra.AvgRMS(13), MAP.Extra.AvgRMS(16)]; % Knee
        Knee.FE.Std = [MAP.Current.StdRMS(13), MAP.Current.StdRMS(16), MAP.AFO.StdRMS(13), MAP.AFO.StdRMS(16), MAP.Last.StdRMS(13), MAP.Last.StdRMS(16), MAP.Extra.StdRMS(13), MAP.Extra.StdRMS(16)];
        Ankle.FE.Avg = [MAP.Current.AvgRMS(19), MAP.Current.AvgRMS(22), MAP.AFO.AvgRMS(19), MAP.AFO.AvgRMS(22), MAP.Last.AvgRMS(19), MAP.Last.AvgRMS(22), MAP.Extra.AvgRMS(19), MAP.Extra.AvgRMS(22)]; % Ankle
        Ankle.FE.Std = [MAP.Current.StdRMS(19), MAP.Current.StdRMS(22), MAP.AFO.StdRMS(19), MAP.AFO.StdRMS(22), MAP.Last.StdRMS(19), MAP.Last.StdRMS(22), MAP.Extra.StdRMS(19), MAP.Extra.StdRMS(22)];
        Ankle.Prog.Avg = [MAP.Current.AvgRMS(20), MAP.Current.AvgRMS(23), MAP.AFO.AvgRMS(20), MAP.AFO.AvgRMS(23), MAP.Last.AvgRMS(20), MAP.Last.AvgRMS(23), MAP.Extra.AvgRMS(20), MAP.Extra.AvgRMS(23)];
        Ankle.Prog.Std = [MAP.Current.StdRMS(20), MAP.Current.StdRMS(23), MAP.AFO.StdRMS(20), MAP.AFO.StdRMS(23), MAP.Last.StdRMS(20), MAP.Last.StdRMS(23), MAP.Extra.StdRMS(20), MAP.Extra.StdRMS(23)];
    end
    
    x = [0 0 2 2];
    % Pelvis
    subplot(3,5,3); grid on; hold on;
    BulletGraph('V', 0, 0, Pelvis.Tilt.Avg, 0, [0 (max(Pelvis.Tilt.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Pelvis.Tilt.Avg,Pelvis.Tilt.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(1) ControlKine.RmsDiffs(1) 0];  f = fill(x, y, 'k'); alpha(f,0.7);    % then the control over top in black
    ylabel('RMS Diff (Degrees)'); 
    if max(Pelvis.Tilt.Avg) > 25
        title('Pelvic Tilt *', 'FontSize',8); ylim([0 ceil(max(Pelvis.Tilt.Avg))+5]);
    else
        title('Pelvic Tilt', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    subplot(3,5,4); grid on; hold on;
    BulletGraph('V', 0, 0, Pelvis.Obl.Avg, 0, [0 (max(Pelvis.Obl.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Pelvis.Obl.Avg,Pelvis.Obl.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(2) ControlKine.RmsDiffs(2) 0];  f = fill(x, y, 'k'); alpha(f,0.7);    % then the control over top in black
    if max(Pelvis.Obl.Avg) > 25
        title('Pelvic Obl *', 'FontSize',8); ylim([0 ceil(max(Pelvis.Obl.Avg))+5]);
    else
        title('Pelvic Obl', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    subplot(3,5,5); grid on; hold on;
    BulletGraph('V', 0, 0, Pelvis.Rot.Avg, 0, [0 (max(Pelvis.Rot.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Pelvis.Rot.Avg,Pelvis.Rot.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(3) ControlKine.RmsDiffs(3) 0];  f = fill(x, y, 'k'); alpha(f,0.7);    % then the control over top in black
    if max(Pelvis.Rot.Avg) > 25
        title('Pelvic Rot *', 'FontSize',8); ylim([0 ceil(max(Pelvis.Rot.Avg))+5]);
    else
        title('Pelvic Rot', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    % Hip
    subplot(3,5,8); grid on; hold on;
    BulletGraph('V', 0, 0, Hip.FE.Avg, 0, [0 (max(Hip.FE.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Hip.FE.Avg,Hip.FE.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(4) ControlKine.RmsDiffs(4) 0];  f = fill(x, y, 'k'); alpha(f,0.7);    % then the control over top in black
    ylabel('RMS Diff (Degrees)');
    if max(Hip.FE.Avg) > 25
        title('Hip F/E *', 'FontSize',8); ylim([0 ceil(max(Hip.FE.Avg))+5]);
    else
        title('Hip F/E', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    subplot(3,5,9); grid on; hold on;
    BulletGraph('V', 0, 0, Hip.AbAd.Avg, 0, [0 (max(Hip.AbAd.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Hip.AbAd.Avg,Hip.AbAd.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(5) ControlKine.RmsDiffs(5) 0];  f = fill(x, y, 'k'); alpha(f,0.7);    % then the control over top in black
    if max(Hip.AbAd.Avg) > 25
        title('Hip Ab/Ad *', 'FontSize',8); ylim([0 ceil(max(Hip.AbAd.Avg))+5]);
    else
        title('Hip Ab/Ad', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    subplot(3,5,10); grid on; hold on;
    BulletGraph('V', 0, 0, Hip.Rot.Avg, 0, [0 (max(Hip.Rot.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Hip.Rot.Avg,Hip.Rot.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(6) ControlKine.RmsDiffs(6) 0]; f = fill(x, y, 'k'); alpha(f,0.7);   % then the control over top in black
    if max(Hip.Rot.Avg) > 25
        title('Hip Rot *', 'FontSize',8); ylim([0 ceil(max(Hip.Rot.Avg))+5]);
    else
        title('Hip Rot', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    % Knee
    subplot(3,5,13); grid on; hold on;
    BulletGraph('V', 0, 0, Knee.FE.Avg, 0, [0 (max(Knee.FE.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Knee.FE.Avg,Knee.FE.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(7) ControlKine.RmsDiffs(7) 0];  f = fill(x, y, 'k'); alpha(f,0.7);    % then the control over top in black
    ylabel('RMS Diff (Degrees)');
    if max(Knee.FE.Avg) > 25
        title('Knee F/E *', 'FontSize',8); ylim([0 ceil(max(Knee.FE.Avg))+5]);
    else
        title('Knee F/E', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    % Ankle
    subplot(3,5,14); grid on; hold on;
    BulletGraph('V', 0, 0, Ankle.FE.Avg, 0, [0 (max(Ankle.FE.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Ankle.FE.Avg,Ankle.FE.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(10) ControlKine.RmsDiffs(10) 0];  f = fill(x, y, 'k'); alpha(f,0.7); % then the control over top in black
    if max(Ankle.FE.Avg) > 25
        title('Ankle D/P *', 'FontSize',8); ylim([0 ceil(max(Ankle.FE.Avg))+5]);
    else
        title('Ankle D/P', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    subplot(3,5,15); grid on; hold on;
    BulletGraph('V', 0, 0, Ankle.Prog.Avg, 0, [0 (max(Ankle.Prog.Avg) + 10)],RedGreen,NumConds*2, 'Yes', 'Yes');
    errorbar(X2,Ankle.Prog.Avg,Ankle.Prog.Std,'k.','Marker','none'); % errorbars
    y = [0 ControlKine.RmsDiffs(11) ControlKine.RmsDiffs(11) 0];  f = fill(x, y, 'k'); alpha(f,0.7);   % then the control over top in black
    if max(Ankle.Prog.Avg) > 25
        title('Foot Prog *', 'FontSize',8); ylim([0 ceil(max(Ankle.Prog.Avg))+5]);
    else
        title('Foot Prog', 'FontSize',8); ylim([0 30]);
    end
    ax = gca; ax.XTick = Xname; ax.XTickLabel = XAxis; ax.XTickLabelRotation = 45; ax.FontSize = 7;
    
    supertitle({'GDI                                                                          Movement Analysis Profile (MAP)   ',' ' }, 'FontSize',10);
    
    saveas(GDIandMAPplot, 'GDI-MAP.png'); 
end


%% Generate GCD of Average Rotations
if strcmp(Type.Current,'Full Kinematics') == 1
    % Calculate Average Rotations
    % uses the averages from every gait cycle in the representative trial, just like GAMS did
    [~, ave_rotations, ~] = GAMS_AnalyTrial(char(KineData(1).RepTrial.Name));
    AveRot = ave_rotations;
    for i = 1:51
        AveRot(1).data(i) = ave_rotations(1).data; % Pelvis
        AveRot(2).data(i) = ave_rotations(2).data;
        AveRot(3).data(i) = ave_rotations(3).data; % Hip
        AveRot(4).data(i) = ave_rotations(4).data;
        if length(ave_rotations) == 12 % if shank rotations are included aka GAIA model
            AveRot(5).data(i) = ave_rotations(11).data; % distal shank rotation
            AveRot(6).data(i) = ave_rotations(12).data;
        else % not a GAIA model - no shank rotations
            AveRot(5).data(i) = ave_rotations(5).data; % Knee
            AveRot(6).data(i) = ave_rotations(6).data;
        end
        AveRot(7).data(i) = ave_rotations(7).data; % Ankle
        AveRot(8).data(i) = ave_rotations(8).data;
    end
    for i = 1:31
        AveRot(9).data(i) = ave_rotations(9).data;
        AveRot(10).data(i) = ave_rotations(10).data;
    end
    % export average rotations in GCD format
    CreateGCD('Ave_Rotations','ave_rot', AveRot);
end

if strcmp(AddAFO,'Yes') == 1 && strcmp(Type.AFO,'Full Kinematics') == 1 % create average rotations for AFO if desired
    % uses the averages from every gait cycle in the representative trial, just like GAMS did
    [~, AFO_ave_rotations, ~] = GAMS_AnalyTrial(char(AFOData(1).RepTrial.Name));
    AFOAveRot = AFO_ave_rotations;
    for i = 1:51
        AFOAveRot(1).data(i) = AFO_ave_rotations(1).data;
        AFOAveRot(2).data(i) = AFO_ave_rotations(2).data;
        AFOAveRot(3).data(i) = AFO_ave_rotations(3).data;
        AFOAveRot(4).data(i) = AFO_ave_rotations(4).data;
        if length(AFO_ave_rotations) == 12 % if shank rotations are included aka GAIA model
            AFOAveRot(5).data(i) = AFO_ave_rotations(11).data; % distal shank rotation
            AFOAveRot(6).data(i) = AFO_ave_rotations(12).data;
        else % not a GAIA model - no shank rotations
            AFOAveRot(5).data(i) = AFO_ave_rotations(5).data; % Knee
            AFOAveRot(6).data(i) = AFO_ave_rotations(6).data;
        end
        AFOAveRot(7).data(i) = AFO_ave_rotations(7).data;
        AFOAveRot(8).data(i) = AFO_ave_rotations(8).data;
    end
    for i = 1:31
        AFOAveRot(9).data(i) = AFO_ave_rotations(9).data;
        AFOAveRot(10).data(i) = AFO_ave_rotations(10).data;
    end
    % export average rotations in GCD format
    CreateGCD('AFO_Ave_Rotations','ave_rot', AFOAveRot);
end


if strcmp(AddLast,'Yes') == 1 && strcmp(Type.Last,'Full Kinematics') == 1 % create average rotations for Last if desired
    % uses the averages from every gait cycle in the representative trial, just like GAMS did
    [~, Last_ave_rotations, ~] = GAMS_AnalyTrial(char(LastData(1).RepTrial.Name));
    Last_AveRot = Last_ave_rotations;
    for i = 1:51
        Last_AveRot(1).data(i) = Last_ave_rotations(1).data;
        Last_AveRot(2).data(i) = Last_ave_rotations(2).data;
        Last_AveRot(3).data(i) = Last_ave_rotations(3).data;
        Last_AveRot(4).data(i) = Last_ave_rotations(4).data;
        if length(Last_ave_rotations) == 12 % if shank rotations are included aka GAIA model
            Last_AveRot(5).data(i) = Last_ave_rotations(11).data; % distal shank rotation
            Last_AveRot(6).data(i) = Last_ave_rotations(12).data;
        else % not a GAIA model - no shank rotations
            Last_AveRot(5).data(i) = Last_ave_rotations(5).data; % Knee
            Last_AveRot(6).data(i) = Last_ave_rotations(6).data;
        end
        Last_AveRot(7).data(i) = Last_ave_rotations(7).data;
        Last_AveRot(8).data(i) = Last_ave_rotations(8).data;
    end
    for i = 1:31
        Last_AveRot(9).data(i) = Last_ave_rotations(9).data;
        Last_AveRot(10).data(i) = Last_ave_rotations(10).data;
    end
    % export average rotations in GCD format
    CreateGCD('Last_Ave_Rotations','ave_rot', Last_AveRot);
end

if strcmp(AddExtra,'Yes') == 1 && strcmp(Type.Extra,'Full Kinematics') == 1 % create average rotations for Extra if desired
    % uses the averages from every gait cycle in the representative trial, just like GAMS did
    [~, Extra_ave_rotations, ~] = GAMS_AnalyTrial(char(ExtraData(1).RepTrial.Name));
    Extra_AveRot = Extra_ave_rotations;
    for i = 1:51
        Extra_AveRot(1).data(i) = Extra_ave_rotations(1).data;
        Extra_AveRot(2).data(i) = Extra_ave_rotations(2).data;
        Extra_AveRot(3).data(i) = Extra_ave_rotations(3).data;
        Extra_AveRot(4).data(i) = Extra_ave_rotations(4).data;
        if length(Extra_ave_rotations) == 12 % if shank rotations are included aka GAIA model
            Extra_AveRot(5).data(i) = Extra_ave_rotations(11).data; % distal shank rotation
            Extra_AveRot(6).data(i) = Extra_ave_rotations(12).data;
        else % not a GAIA model - no shank rotations
            Extra_AveRot(5).data(i) = Extra_ave_rotations(5).data; % Knee
            Extra_AveRot(6).data(i) = Extra_ave_rotations(6).data;
        end
        Extra_AveRot(7).data(i) = Extra_ave_rotations(7).data;
        Extra_AveRot(8).data(i) = Extra_ave_rotations(8).data;
    end
    for i = 1:31
        Extra_AveRot(9).data(i) = Extra_ave_rotations(9).data;
        Extra_AveRot(10).data(i) = Extra_ave_rotations(10).data;
    end
    % export average rotations in GCD format
    CreateGCD('Extra_Ave_Rotations','ave_rot', Extra_AveRot);
end

%% Create Table of Average Gait Performance and Temporal Spatial Variables
if NumConds == 1 % for 1 measure
    T = table;
    T.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}];
    T.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.Avg];
    T.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.Avg];
    T.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.Avg];
    T.L_GDI = [{'N/A'}; 100; KineData(1).GDI(1)];
    T.R_GDI = [{'N/A'}; 100; KineData(1).GDI(2)];
    T.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.Avg];
    T.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.Avg];
    T.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.Avg];
    T.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.Avg + Stance.Current.Right.Avg)/2];
    T.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.Avg + Swing.Current.Right.Avg)/2];
    T.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.Avg + DS1.Current.Right.Avg)/2];
    T.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.Avg + SLS.Current.Right.Avg)/2];
    T.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.Avg + DS2.Current.Right.Avg)/2];
elseif NumConds == 2
    if strcmp(AddAFO,'Yes') == 1
        T = table;
        T.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.AFO)}];
        T.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.Avg; Cadence.AFO.Avg];
        T.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.Avg; StrideLength.AFO.Avg];
        T.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.Avg; GaitSpeed.AFO.Avg];
        T.L_GDI = [{'N/A'}; 100; KineData(1).GDI(1); AFOData(1).GDI(1)];
        T.R_GDI = [{'N/A'}; 100; KineData(1).GDI(2); AFOData(1).GDI(2)];
        T.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.Avg; StepLength.AFO.Left.Avg];
        T.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.Avg; StepLength.AFO.Right.Avg];
        T.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.Avg; StepWidth.AFO.Avg];
        T.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.Avg + Stance.Current.Right.Avg)/2; (Stance.AFO.Left.Avg + Stance.AFO.Right.Avg)/2];
        T.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.Avg + Swing.Current.Right.Avg)/2; (Swing.AFO.Left.Avg + Swing.AFO.Right.Avg)/2];
        T.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.Avg + DS1.Current.Right.Avg)/2; (DS1.AFO.Left.Avg + DS1.AFO.Right.Avg)/2];
        T.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.Avg + SLS.Current.Right.Avg)/2; (SLS.AFO.Left.Avg + SLS.AFO.Right.Avg)/2];
        T.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.Avg + DS2.Current.Right.Avg)/2; (DS2.AFO.Left.Avg + DS2.AFO.Right.Avg)/2];
    elseif strcmp(AddLast,'Yes') == 1
        T = table;
        T.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.Last)}];
        T.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.Avg; Cadence.Last.Avg];
        T.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.Avg; StrideLength.Last.Avg];
        T.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.Avg; GaitSpeed.Last.Avg];
        T.L_GDI = [{'N/A'}; 100; KineData(1).GDI(1); LastData(1).GDI(1)];
        T.R_GDI = [{'N/A'}; 100; KineData(1).GDI(2); LastData(1).GDI(2)];
        T.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.Avg; StepLength.Last.Left.Avg];
        T.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.Avg; StepLength.Last.Right.Avg];
        T.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.Avg; StepWidth.Last.Avg];
        T.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.Avg + Stance.Current.Right.Avg)/2; (Stance.Last.Left.Avg + Stance.Last.Right.Avg)/2];
        T.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.Avg + Swing.Current.Right.Avg)/2; (Swing.Last.Left.Avg + Swing.Last.Right.Avg)/2];
        T.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3);  (DS1.Current.Left.Avg + DS1.Current.Right.Avg)/2; (DS1.Last.Left.Avg + DS1.Last.Right.Avg)/2];
        T.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.Avg + SLS.Current.Right.Avg)/2; (SLS.Last.Left.Avg + SLS.Last.Right.Avg)/2];
        T.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.Avg + DS2.Current.Right.Avg)/2; (DS2.Last.Left.Avg + DS2.Last.Right.Avg)/2];
    elseif strcmp(AddExtra,'Yes') == 1
        T = table;
        T.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.Extra)}];
        T.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.Avg; Cadence.Extra.Avg];
        T.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.Avg; StrideLength.Extra.Avg];
        T.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.Avg; GaitSpeed.Extra.Avg];
        T.L_GDI = [{'N/A'}; 100;  KineData(1).GDI(1); ExtraData(1).GDI(1)];
        T.R_GDI = [{'N/A'}; 100; KineData(1).GDI(2); ExtraData(1).GDI(2)];
        T.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.Avg; StepLength.Extra.Left.Avg];
        T.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.Avg; StepLength.Extra.Right.Avg];
        T.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.Avg; StepWidth.Extra.Avg];
        T.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.Avg + Stance.Current.Right.Avg)/2; (Stance.Extra.Left.Avg + Stance.Extra.Right.Avg)/2];
        T.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.Avg + Swing.Current.Right.Avg)/2; (Swing.Extra.Left.Avg + Swing.Extra.Right.Avg)/2];
        T.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.Avg + DS1.Current.Right.Avg)/2; (DS1.Extra.Left.Avg + DS1.Extra.Right.Avg)/2];
        T.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.Avg + SLS.Current.Right.Avg)/2; (SLS.Extra.Left.Avg + SLS.Extra.Right.Avg)/2];
        T.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.Avg + DS2.Current.Right.Avg)/2; (DS2.Extra.Left.Avg + DS2.Extra.Right.Avg)/2];
    end
elseif NumConds == 3
    if strcmp(AddExtra,'No') == 1
        T = table;
        T.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.AFO)}; {char(Name.Last)}];
        T.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.Avg; Cadence.AFO.Avg; Cadence.Last.Avg];
        T.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.Avg; StrideLength.AFO.Avg; StrideLength.Last.Avg];
        T.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.Avg; GaitSpeed.AFO.Avg; GaitSpeed.Last.Avg];
        T.L_GDI = [{'N/A'}; 100; KineData(1).GDI(1); AFOData(1).GDI(1); LastData(1).GDI(1)];
        T.R_GDI = [{'N/A'}; 100; KineData(1).GDI(2); AFOData(1).GDI(2); LastData(1).GDI(2)];
        T.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.Avg; StepLength.AFO.Left.Avg; StepLength.Last.Left.Avg];
        T.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.Avg; StepLength.AFO.Right.Avg; StepLength.Last.Right.Avg];
        T.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.Avg; StepWidth.AFO.Avg; StepWidth.Last.Avg];
        T.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.Avg + Stance.Current.Right.Avg)/2; (Stance.AFO.Left.Avg + Stance.AFO.Right.Avg)/2; (Stance.Last.Left.Avg + Stance.Last.Right.Avg)/2];
        T.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.Avg + Swing.Current.Right.Avg)/2; (Swing.AFO.Left.Avg + Swing.AFO.Right.Avg)/2; (Swing.Last.Left.Avg + Swing.Last.Right.Avg)/2];
        T.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.Avg + DS1.Current.Right.Avg)/2; (DS1.AFO.Left.Avg + DS1.AFO.Right.Avg)/2; (DS1.Last.Left.Avg + DS1.Last.Right.Avg)/2];
        T.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.Avg + SLS.Current.Right.Avg)/2; (SLS.AFO.Left.Avg + SLS.AFO.Right.Avg)/2; (SLS.Last.Left.Avg + SLS.Last.Right.Avg)/2];
        T.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.Avg + DS2.Current.Right.Avg)/2; (DS2.AFO.Left.Avg + DS2.AFO.Right.Avg)/2; (DS2.Last.Left.Avg + DS2.Last.Right.Avg)/2];
    elseif strcmp(AddLast,'No') == 1
        T = table;
        T.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.AFO)}; {char(Name.Extra)}];
        T.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.Avg; Cadence.AFO.Avg; Cadence.Extra.Avg];
        T.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.Avg; StrideLength.AFO.Avg; StrideLength.Extra.Avg];
        T.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.Avg; GaitSpeed.AFO.Avg; GaitSpeed.Extra.Avg];
        T.L_GDI = [{'N/A'}; 100; KineData(1).GDI(1); AFOData(1).GDI(1); ExtraData(1).GDI(1)];
        T.R_GDI = [{'N/A'}; 100; KineData(1).GDI(2); AFOData(1).GDI(2); ExtraData(1).GDI(2)];
        T.L_StepLength = [{'cm'}; TSN.TimeDist(4);StepLength.Current.Left.Avg; StepLength.AFO.Left.Avg; StepLength.Extra.Left.Avg];
        T.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.Avg; StepLength.AFO.Right.Avg; StepLength.Extra.Right.Avg];
        T.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.Avg; StepWidth.AFO.Avg; StepWidth.Extra.Avg];
        T.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.Avg + Stance.Current.Right.Avg)/2; (Stance.AFO.Left.Avg + Stance.AFO.Right.Avg)/2; (Stance.Extra.Left.Avg + Stance.Extra.Right.Avg)/2];
        T.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.Avg + Swing.Current.Right.Avg)/2; (Swing.AFO.Left.Avg + Swing.AFO.Right.Avg)/2; (Swing.Extra.Left.Avg + Swing.Extra.Right.Avg)/2];
        T.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.Avg + DS1.Current.Right.Avg)/2; (DS1.AFO.Left.Avg + DS1.AFO.Right.Avg)/2; (DS1.Extra.Left.Avg + DS1.Extra.Right.Avg)/2];
        T.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.Avg + SLS.Current.Right.Avg)/2; (SLS.AFO.Left.Avg + SLS.AFO.Right.Avg)/2; (SLS.Extra.Left.Avg + SLS.Extra.Right.Avg)/2];
        T.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.Avg + DS2.Current.Right.Avg)/2; (DS2.AFO.Left.Avg + DS2.AFO.Right.Avg)/2; (DS2.Extra.Left.Avg + DS2.Extra.Right.Avg)/2];
    elseif strcmp(AddAFO,'No') == 1
        T = table;
        T.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.Extra)}; {char(Name.Last)}];
        T.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.Avg; Cadence.Extra.Avg; Cadence.Last.Avg];
        T.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.Avg; StrideLength.Extra.Avg; StrideLength.Last.Avg];
        T.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.Avg; GaitSpeed.Extra.Avg; GaitSpeed.Last.Avg];
        T.L_GDI = [{'N/A'}; 100; KineData(1).GDI(1); ExtraData(1).GDI(1); LastData(1).GDI(1)];
        T.R_GDI = [{'N/A'}; 100; KineData(1).GDI(2); ExtraData(1).GDI(2); LastData(1).GDI(2)];
        T.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.Avg; StepLength.Extra.Left.Avg; StepLength.Last.Left.Avg];
        T.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.Avg; StepLength.Extra.Right.Avg; StepLength.Last.Right.Avg];
        T.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.Avg; StepWidth.Extra.Avg; StepWidth.Last.Avg];
        T.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.Avg + Stance.Current.Right.Avg)/2; (Stance.Extra.Left.Avg + Stance.Extra.Right.Avg)/2; (Stance.Last.Left.Avg + Stance.Last.Right.Avg)/2];
        T.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.Avg + Swing.Current.Right.Avg)/2; (Swing.Extra.Left.Avg + Swing.Extra.Right.Avg)/2; (Swing.Last.Left.Avg + Swing.Last.Right.Avg)/2];
        T.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.Avg + DS1.Current.Right.Avg)/2; (DS1.Extra.Left.Avg + DS1.Extra.Right.Avg)/2; (DS1.Last.Left.Avg + DS1.Last.Right.Avg)/2];
        T.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.Avg + SLS.Current.Right.Avg)/2; (SLS.Extra.Left.Avg + SLS.Extra.Right.Avg)/2; (SLS.Last.Left.Avg + SLS.Last.Right.Avg)/2];
        T.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.Avg + DS2.Current.Right.Avg)/2; (DS2.Extra.Left.Avg + DS2.Extra.Right.Avg)/2; (DS2.Last.Left.Avg + DS2.Last.Right.Avg)/2];
    end
elseif NumConds == 4 % If 4 conditions present
    T = table;
    T.Conditions = [{'Units'};{'Reference'};{char(Name.Curr)}; {char(Name.AFO)}; {char(Name.Last)}; {char(Name.Extra)}];
    T.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.Avg; Cadence.AFO.Avg; Cadence.Last.Avg; Cadence.Extra.Avg];
    T.StepLength = [{'m'};TSN.TimeDist(2); StrideLength.Current.Avg; StrideLength.AFO.Avg; StrideLength.Last.Avg; StrideLength.Extra.Avg];
    T.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.Avg; GaitSpeed.AFO.Avg; GaitSpeed.Last.Avg; GaitSpeed.Extra.Avg];
    T.L_GDI = [{'N/A'}; 100; KineData(1).GDI(1); AFOData(1).GDI(1); LastData(1).GDI(1); ExtraData(1).GDI(1)];
    T.R_GDI = [{'N/A'}; 100; KineData(1).GDI(2); AFOData(1).GDI(2); LastData(1).GDI(2); ExtraData(1).GDI(2)];
    T.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.Avg; StepLength.AFO.Left.Avg; StepLength.Last.Left.Avg; StepLength.Extra.Left.Avg];
    T.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.Avg; StepLength.AFO.Right.Avg; StepLength.Last.Right.Avg; StepLength.Extra.Right.Avg];
    T.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.Avg; StepWidth.AFO.Avg; StepWidth.Last.Avg; StepWidth.Extra.Avg];
    T.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.Avg + Stance.Current.Right.Avg)/2; (Stance.AFO.Left.Avg + Stance.AFO.Right.Avg)/2; (Stance.Last.Left.Avg + Stance.Last.Right.Avg)/2; (Stance.Extra.Left.Avg + Stance.Extra.Right.Avg)/2];
    T.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.Avg + Swing.Current.Right.Avg)/2; (Swing.AFO.Left.Avg + Swing.AFO.Right.Avg)/2; (Swing.Last.Left.Avg + Swing.Last.Right.Avg)/2; (Swing.Extra.Left.Avg + Swing.Extra.Right.Avg)/2];
    T.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.Avg + DS1.Current.Right.Avg)/2; (DS1.AFO.Left.Avg + DS1.AFO.Right.Avg)/2; (DS1.Last.Left.Avg + DS1.Last.Right.Avg)/2; (DS1.Extra.Left.Avg + DS1.Extra.Right.Avg)/2];
    T.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.Avg + SLS.Current.Right.Avg)/2; (SLS.AFO.Left.Avg + SLS.AFO.Right.Avg)/2; (SLS.Last.Left.Avg + SLS.Last.Right.Avg)/2; (SLS.Extra.Left.Avg + SLS.Extra.Right.Avg)/2];
    T.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.Avg + DS2.Current.Right.Avg)/2; (DS2.AFO.Left.Avg + DS2.AFO.Right.Avg)/2; (DS2.Last.Left.Avg + DS2.Last.Right.Avg)/2; (DS2.Extra.Left.Avg + DS2.Extra.Right.Avg)/2];
end

%% Create Table of STD Gait Performance and Temporal Spatial Variables
if NumConds == 1 % for 1 measure
    Tstd = table;
    Tstd.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}];
    Tstd.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.STD];
    Tstd.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.STD];
    Tstd.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.STD];
    Tstd.L_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(1))];
    Tstd.R_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(2))];
    Tstd.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.STD];
    Tstd.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.STD];
    Tstd.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.STD];
    Tstd.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.STD + Stance.Current.Right.STD)/2];
    Tstd.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.STD + Swing.Current.Right.STD)/2];
    Tstd.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.STD + DS1.Current.Right.STD)/2];
    Tstd.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.STD + SLS.Current.Right.STD)/2];
    Tstd.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.STD + DS2.Current.Right.STD)/2];
elseif NumConds == 2
    if strcmp(AddAFO,'Yes') == 1
        Tstd = table;
        Tstd.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.AFO)}];
        Tstd.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.STD; Cadence.AFO.STD];
        Tstd.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.STD; StrideLength.AFO.STD];
        Tstd.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.STD; GaitSpeed.AFO.STD];
        Tstd.L_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(1)); std(AFOData(1).TrialsGDI(1))];
        Tstd.R_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(2)); std(AFOData(1).TrialsGDI(2))];
        Tstd.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.STD; StepLength.AFO.Left.STD];
        Tstd.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.STD; StepLength.AFO.Right.STD];
        Tstd.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.STD; StepWidth.AFO.STD];
        Tstd.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.STD + Stance.Current.Right.STD)/2; (Stance.AFO.Left.STD + Stance.AFO.Right.STD)/2];
        Tstd.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.STD + Swing.Current.Right.STD)/2; (Swing.AFO.Left.STD + Swing.AFO.Right.STD)/2];
        Tstd.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.STD + DS1.Current.Right.STD)/2; (DS1.AFO.Left.STD + DS1.AFO.Right.STD)/2];
        Tstd.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.STD + SLS.Current.Right.STD)/2; (SLS.AFO.Left.STD + SLS.AFO.Right.STD)/2];
        Tstd.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.STD + DS2.Current.Right.STD)/2; (DS2.AFO.Left.STD + DS2.AFO.Right.STD)/2];
    elseif strcmp(AddLast,'Yes') == 1
        Tstd = table;
        Tstd.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.Last)}];
        Tstd.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.STD; Cadence.Last.STD];
        Tstd.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.STD; StrideLength.Last.STD];
        Tstd.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.STD; GaitSpeed.Last.STD];
        Tstd.L_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(1)); std(LastData(1).TrialsGDI(1))];
        Tstd.R_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(2)); std(LastData(1).TrialsGDI(2))];
        Tstd.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.STD; StepLength.Last.Left.STD];
        Tstd.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.STD; StepLength.Last.Right.STD];
        Tstd.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.STD; StepWidth.Last.STD];
        Tstd.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.STD + Stance.Current.Right.STD)/2; (Stance.Last.Left.STD + Stance.Last.Right.STD)/2];
        Tstd.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.STD + Swing.Current.Right.STD)/2; (Swing.Last.Left.STD + Swing.Last.Right.STD)/2];
        Tstd.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3);  (DS1.Current.Left.STD + DS1.Current.Right.STD)/2; (DS1.Last.Left.STD + DS1.Last.Right.STD)/2];
        Tstd.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.STD + SLS.Current.Right.STD)/2; (SLS.Last.Left.STD + SLS.Last.Right.STD)/2];
        Tstd.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.STD + DS2.Current.Right.STD)/2; (DS2.Last.Left.STD + DS2.Last.Right.STD)/2];
    elseif strcmp(AddExtra,'Yes') == 1
        Tstd = table;
        Tstd.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.Extra)}];
        Tstd.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.STD; Cadence.Extra.STD];
        Tstd.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.STD; StrideLength.Extra.STD];
        Tstd.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.STD; GaitSpeed.Extra.STD];
        Tstd.L_GDI = [{'N/A'}; 100;  std(KineData(1).TrialsGDI(1)); std(ExtraData(1).TrialsGDI(1))];
        Tstd.R_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(2)); std(ExtraData(1).TrialsGDI(2))];
        Tstd.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.STD; StepLength.Extra.Left.STD];
        Tstd.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.STD; StepLength.Extra.Right.STD];
        Tstd.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.STD; StepWidth.Extra.STD];
        Tstd.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.STD + Stance.Current.Right.STD)/2; (Stance.Extra.Left.STD + Stance.Extra.Right.STD)/2];
        Tstd.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.STD + Swing.Current.Right.STD)/2; (Swing.Extra.Left.STD + Swing.Extra.Right.STD)/2];
        Tstd.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.STD + DS1.Current.Right.STD)/2; (DS1.Extra.Left.STD + DS1.Extra.Right.STD)/2];
        Tstd.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.STD + SLS.Current.Right.STD)/2; (SLS.Extra.Left.STD + SLS.Extra.Right.STD)/2];
        Tstd.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.STD + DS2.Current.Right.STD)/2; (DS2.Extra.Left.STD + DS2.Extra.Right.STD)/2];
    end
elseif NumConds == 3
    if strcmp(AddExtra,'No') == 1
        Tstd = table;
        Tstd.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.AFO)}; {char(Name.Last)}];
        Tstd.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.STD; Cadence.AFO.STD; Cadence.Last.STD];
        Tstd.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.STD; StrideLength.AFO.STD; StrideLength.Last.STD];
        Tstd.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.STD; GaitSpeed.AFO.STD; GaitSpeed.Last.STD];
        Tstd.L_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(1)); std(AFOData(1).TrialsGDI(1)); std(LastData(1).TrialsGDI(1))];
        Tstd.R_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(2)); std(AFOData(1).TrialsGDI(2)); std(LastData(1).TrialsGDI(2))];
        Tstd.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.STD; StepLength.AFO.Left.STD; StepLength.Last.Left.STD];
        Tstd.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.STD; StepLength.AFO.Right.STD; StepLength.Last.Right.STD];
        Tstd.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.STD; StepWidth.AFO.STD; StepWidth.Last.STD];
        Tstd.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.STD + Stance.Current.Right.STD)/2; (Stance.AFO.Left.STD + Stance.AFO.Right.STD)/2; (Stance.Last.Left.STD + Stance.Last.Right.STD)/2];
        Tstd.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.STD + Swing.Current.Right.STD)/2; (Swing.AFO.Left.STD + Swing.AFO.Right.STD)/2; (Swing.Last.Left.STD + Swing.Last.Right.STD)/2];
        Tstd.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.STD + DS1.Current.Right.STD)/2; (DS1.AFO.Left.STD + DS1.AFO.Right.STD)/2; (DS1.Last.Left.STD + DS1.Last.Right.STD)/2];
        Tstd.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.STD + SLS.Current.Right.STD)/2; (SLS.AFO.Left.STD + SLS.AFO.Right.STD)/2; (SLS.Last.Left.STD + SLS.Last.Right.STD)/2];
        Tstd.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.STD + DS2.Current.Right.STD)/2; (DS2.AFO.Left.STD + DS2.AFO.Right.STD)/2; (DS2.Last.Left.STD + DS2.Last.Right.STD)/2];
    elseif strcmp(AddLast,'No') == 1
        Tstd = table;
        Tstd.Conditions = [{'Units'};{'Reference'}; {char(Name.Curr)}; {char(Name.AFO)}; {char(Name.Extra)}];
        Tstd.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.STD; Cadence.AFO.STD; Cadence.Extra.STD];
        Tstd.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.STD; StrideLength.AFO.STD; StrideLength.Extra.STD];
        Tstd.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.STD; GaitSpeed.AFO.STD; GaitSpeed.Extra.STD];
        Tstd.L_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(1)); std(AFOData(1).TrialsGDI(1)); std(ExtraData(1).TrialsGDI(1))];
        Tstd.R_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(2)); std(AFOData(1).TrialsGDI(2)); std(ExtraData(1).TrialsGDI(2))];
        Tstd.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.STD; StepLength.AFO.Left.STD; StepLength.Extra.Left.STD];
        Tstd.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.STD; StepLength.AFO.Right.STD; StepLength.Extra.Right.STD];
        Tstd.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.STD; StepWidth.AFO.STD; StepWidth.Extra.STD];
        Tstd.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.STD + Stance.Current.Right.STD)/2; (Stance.AFO.Left.STD + Stance.AFO.Right.STD)/2; (Stance.Extra.Left.STD + Stance.Extra.Right.STD)/2];
        Tstd.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.STD + Swing.Current.Right.STD)/2; (Swing.AFO.Left.STD + Swing.AFO.Right.STD)/2; (Swing.Extra.Left.STD + Swing.Extra.Right.STD)/2];
        Tstd.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.STD + DS1.Current.Right.STD)/2; (DS1.AFO.Left.STD + DS1.AFO.Right.STD)/2; (DS1.Extra.Left.STD + DS1.Extra.Right.STD)/2];
        Tstd.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.STD + SLS.Current.Right.STD)/2; (SLS.AFO.Left.STD + SLS.AFO.Right.STD)/2; (SLS.Extra.Left.STD + SLS.Extra.Right.STD)/2];
        Tstd.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.STD + DS2.Current.Right.STD)/2; (DS2.AFO.Left.STD + DS2.AFO.Right.STD)/2; (DS2.Extra.Left.STD + DS2.Extra.Right.STD)/2];
    elseif strcmp(AddAFO,'No') == 1
        Tstd = table;
        Tstd.Conditions = [{'Units'}; {'Reference'}; {char(Name.Curr)}; {char(Name.Extra)}; {char(Name.Last)}];
        Tstd.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.STD; Cadence.Extra.STD; Cadence.Last.STD];
        Tstd.StepLength = [{'m'}; TSN.TimeDist(2); StrideLength.Current.STD; StrideLength.Extra.STD; StrideLength.Last.STD];
        Tstd.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.STD; GaitSpeed.Extra.STD; GaitSpeed.Last.STD];
        Tstd.L_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(1)); std(ExtraData(1).TrialsGDI(1)); std(LastData(1).TrialsGDI(1))];
        Tstd.R_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(2)); std(ExtraData(1).TrialsGDI(2)); std(LastData(1).TrialsGDI(2))];
        Tstd.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.STD; StepLength.Extra.Left.STD; StepLength.Last.Left.STD];
        Tstd.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.STD; StepLength.Extra.Right.STD; StepLength.Last.Right.STD];
        Tstd.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.STD; StepWidth.Extra.STD; StepWidth.Last.STD];
        Tstd.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.STD + Stance.Current.Right.STD)/2; (Stance.Extra.Left.STD + Stance.Extra.Right.STD)/2; (Stance.Last.Left.STD + Stance.Last.Right.STD)/2];
        Tstd.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.STD + Swing.Current.Right.STD)/2; (Swing.Extra.Left.STD + Swing.Extra.Right.STD)/2; (Swing.Last.Left.STD + Swing.Last.Right.STD)/2];
        Tstd.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.STD + DS1.Current.Right.STD)/2; (DS1.Extra.Left.STD + DS1.Extra.Right.STD)/2; (DS1.Last.Left.STD + DS1.Last.Right.STD)/2];
        Tstd.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.STD + SLS.Current.Right.STD)/2; (SLS.Extra.Left.STD + SLS.Extra.Right.STD)/2; (SLS.Last.Left.STD + SLS.Last.Right.STD)/2];
        Tstd.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.STD + DS2.Current.Right.STD)/2; (DS2.Extra.Left.STD + DS2.Extra.Right.STD)/2; (DS2.Last.Left.STD + DS2.Last.Right.STD)/2];
    end
elseif NumConds == 4 % If 4 conditions present
    Tstd = table;
    Tstd.Conditions = [{'Units'};{'Reference'};{char(Name.Curr)}; {char(Name.AFO)}; {char(Name.Last)}; {char(Name.Extra)}];
    Tstd.Cadence = [{'steps/min'}; TSN.TimeDist(1); Cadence.Current.STD; Cadence.AFO.STD; Cadence.Last.STD; Cadence.Extra.STD];
    Tstd.StepLength = [{'m'};TSN.TimeDist(2); StrideLength.Current.STD; StrideLength.AFO.STD; StrideLength.Last.STD; StrideLength.Extra.STD];
    Tstd.GaitSpeed = [{'m/min'}; TSN.TimeDist(3); GaitSpeed.Current.STD; GaitSpeed.AFO.STD; GaitSpeed.Last.STD; GaitSpeed.Extra.STD];
    Tstd.L_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(1)); std(AFOData(1).TrialsGDI(1)); std(LastData(1).TrialsGDI(1)); std(ExtraData(1).TrialsGDI(1))];
    Tstd.R_GDI = [{'N/A'}; 100; std(KineData(1).TrialsGDI(2)); std(AFOData(1).TrialsGDI(2)); std(LastData(1).TrialsGDI(2)); std(ExtraData(1).TrialsGDI(2))];
    Tstd.L_StepLength = [{'cm'}; TSN.TimeDist(4); StepLength.Current.Left.STD; StepLength.AFO.Left.STD; StepLength.Last.Left.STD; StepLength.Extra.Left.STD];
    Tstd.R_StepLength = [{'cm'}; TSN.TimeDist(5); StepLength.Current.Right.STD; StepLength.AFO.Right.STD; StepLength.Last.Right.STD; StepLength.Extra.Right.STD];
    Tstd.StepWidth = [{'cm'}; TSN.TimeDist(6); StepWidth.Current.STD; StepWidth.AFO.STD; StepWidth.Last.STD; StepWidth.Extra.STD];
    Tstd.Stance = [{'% Gait Cycle'}; TSN.StanceSwing(1); (Stance.Current.Left.STD + Stance.Current.Right.STD)/2; (Stance.AFO.Left.STD + Stance.AFO.Right.STD)/2; (Stance.Last.Left.STD + Stance.Last.Right.STD)/2; (Stance.Extra.Left.STD + Stance.Extra.Right.STD)/2];
    Tstd.Swing = [{'% Gait Cycle'}; TSN.StanceSwing(2); (Swing.Current.Left.STD + Swing.Current.Right.STD)/2; (Swing.AFO.Left.STD + Swing.AFO.Right.STD)/2; (Swing.Last.Left.STD + Swing.Last.Right.STD)/2; (Swing.Extra.Left.STD + Swing.Extra.Right.STD)/2];
    Tstd.Initial_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(3); (DS1.Current.Left.STD + DS1.Current.Right.STD)/2; (DS1.AFO.Left.STD + DS1.AFO.Right.STD)/2; (DS1.Last.Left.STD + DS1.Last.Right.STD)/2; (DS1.Extra.Left.STD + DS1.Extra.Right.STD)/2];
    Tstd.Single_Limb_Support = [{'% Gait Cycle'}; TSN.StanceSwing(4); (SLS.Current.Left.STD + SLS.Current.Right.STD)/2; (SLS.AFO.Left.STD + SLS.AFO.Right.STD)/2; (SLS.Last.Left.STD + SLS.Last.Right.STD)/2; (SLS.Extra.Left.STD + SLS.Extra.Right.STD)/2];
    Tstd.Final_Double_Support = [{'% Gait Cycle'}; TSN.StanceSwing(5); (DS2.Current.Left.STD + DS2.Current.Right.STD)/2; (DS2.AFO.Left.STD + DS2.AFO.Right.STD)/2; (DS2.Last.Left.STD + DS2.Last.Right.STD)/2; (DS2.Extra.Left.STD + DS2.Extra.Right.STD)/2];
end

Tstd(2,:) = []; % remove reference row due to lack of STD norms
Tstd(1,:) = []; % remove Units row

%% Calculate #s of trials and cycles
clearvars NumTrials MinCycles h prompt
% Current condition
if strcmp(Type.Current, 'Full Kinematics') == 1
    for i = 1:length(KineData)
        LCycles.Current(i) = KineData(i).NumCycles(1);
        RCycles.Current(i) = KineData(i).NumCycles(2);
    end
    LCycles.Currentsum = sum(LCycles.Current);
    RCycles.Currentsum = sum(RCycles.Current);
    NumTrials.CurrentFullCyc = length(KineData(1).TrialsGDI)  - sum(isnan(KineData(1).TrialsGDI(:,1))); 
else
    for i = 1:length(KineData)
        LCycles.Current(i) = KineData(i).TempSpat(34).data;
        RCycles.Current(i) = KineData(i).TempSpat(35).data;
        Mins.Current(i) = min(LCycles.Current(i), RCycles.Current(i)); 
    end
    LCycles.Currentsum = sum(LCycles.Current);
    RCycles.Currentsum = sum(RCycles.Current);
    NumTrials.CurrentFullCyc = sum(Mins.Current);
end
NumTrials.Current = length(KineData);

% AFO condition
if strcmp(AddAFO,'Yes') == 1
    if strcmp(Type.AFO, 'Full Kinematics') == 1
        for i = 1:length(AFOData)
            LCycles.AFO(i) = AFOData(i).NumCycles(1);
            RCycles.AFO(i) = AFOData(i).NumCycles(2);
        end
        LCycles.AFOsum = sum(LCycles.AFO);
        RCycles.AFOsum = sum(RCycles.AFO);
        NumTrials.AFOFullCyc = length(AFOData(1).TrialsGDI) - sum(isnan(AFOData(1).TrialsGDI(:,1)));
    else
        for i = 1:length(AFOData)
            LCycles.AFO(i) = AFOData(i).TempSpat(34).data;
            RCycles.AFO(i) = AFOData(i).TempSpat(35).data;
            Mins.AFO(i) = min(LCycles.AFO(i), RCycles.AFO(i)); 
        end
        LCycles.AFOsum = sum(LCycles.AFO);
        RCycles.AFOsum = sum(RCycles.AFO);
        NumTrials.AFOFullCyc = sum(Mins.AFO);
    end
    NumTrials.AFO = length(AFOData);
end
% Last condition
if strcmp(AddLast,'Yes') == 1
    if strcmp(Type.Last, 'Full Kinematics') == 1
        for i = 1:length(LastData)
            LCycles.Last(i) = LastData(i).NumCycles(1);
            RCycles.Last(i) = LastData(i).NumCycles(2);
        end
        LCycles.Lastsum = sum(LCycles.Last);
        RCycles.Lastsum = sum(RCycles.Last);
        NumTrials.LastFullCyc = length(LastData(1).TrialsGDI)  - sum(isnan(LastData(1).TrialsGDI(:,1)));
    else
        for i = 1:length(LastData)
            LCycles.Last(i) = LastData(i).TempSpat(34).data;
            RCycles.Last(i) = LastData(i).TempSpat(35).data;
            Mins.Last(i) = min(LCycles.Last(i), RCycles.Last(i)); 
        end
        LCycles.Lastsum = sum(LCycles.Last);
        RCycles.Lastsum = sum(RCycles.Last);
        NumTrials.LastFullCyc = sum(Mins.Last);
    end
    NumTrials.Last = length(LastData);
end
% Extra condition
if strcmp(AddExtra,'Yes') == 1
    if strcmp(Type.Extra, 'Full Kinematics') == 1
        for i = 1:length(ExtraData)
            LCycles.Extra(i) = ExtraData(i).NumCycles(1);
            RCycles.Extra(i) = ExtraData(i).NumCycles(2);
        end
        LCycles.Extrasum = sum(LCycles.Extra);
        RCycles.Extrasum = sum(RCycles.Extra);
        NumTrials.ExtraFullCyc = length(ExtraData(1).TrialsGDI)  - sum(isnan(ExtraData(1).TrialsGDI(:,1)));
    else
        for i = 1:length(ExtraData)
            LCycles.Extra(i) = ExtraData(i).TempSpat(34).data;
            RCycles.Extra(i) = ExtraData(i).TempSpat(35).data;
            Mins.Extra(i) = min(LCycles.Extra(i), RCycles.Extra(i)); 
        end
        LCycles.Extrasum = sum(LCycles.Extra);
        RCycles.Extrasum = sum(RCycles.Extra);
        NumTrials.ExtraFullCyc = sum(Mins.Extra); 
    end
    NumTrials.Extra = length(ExtraData);
end

%% Create table of number of trials and cycles for each condition
    TTrialCyc = table;
if NumConds == 1 % for 1 condition
    TTrialCyc.Conditions = Tstd.Conditions;
    TTrialCyc.NumTrials = NumTrials.Current;
    TTrialCyc.LCyc = LCycles.Currentsum; 
    TTrialCyc.RCyc = RCycles.Currentsum;
    TTrialCyc.FullCyc = NumTrials.CurrentFullCyc;
elseif NumConds == 2 % 2 conditions
    if strcmp(AddAFO, 'Yes') == 1
        TTrialCyc.Conditions = Tstd.Conditions;
        TTrialCyc.NumTrials = [NumTrials.Current; NumTrials.AFO] ;
        TTrialCyc.LCyc = [LCycles.Currentsum; LCycles.AFOsum];
        TTrialCyc.RCyc = [RCycles.Currentsum; RCycles.AFOsum];
        TTrialCyc.FullCyc = [NumTrials.CurrentFullCyc; NumTrials.AFOFullCyc];
    elseif strcmp(AddLast, 'Yes') == 1
        TTrialCyc.Conditions = Tstd.Conditions;
        TTrialCyc.NumTrials = [NumTrials.Current; NumTrials.Last] ;
        TTrialCyc.LCyc = [LCycles.Currentsum; LCycles.Lastsum];
        TTrialCyc.RCyc = [RCycles.Currentsum; RCycles.Lastsum];
        TTrialCyc.FullCyc = [NumTrials.CurrentFullCyc; NumTrials.LastFullCyc];
    elseif strcmp(AddExtra, 'Yes') == 1
        TTrialCyc.Conditions = Tstd.Conditions;
        TTrialCyc.NumTrials = [NumTrials.Current; NumTrials.Extra];
        TTrialCyc.LCyc = [LCycles.Currentsum; LCycles.Extrasum];
        TTrialCyc.RCyc = [RCycles.Currentsum; RCycles.Extrasum];
        TTrialCyc.FullCyc = [NumTrials.CurrentFullCyc; NumTrials.ExtraFullCyc];
    end
elseif NumConds == 3 % 3 conditions
    if strcmp(AddAFO,'No') == 1
        TTrialCyc.Conditions = Tstd.Conditions;
        TTrialCyc.NumTrials = [NumTrials.Current; NumTrials.Last; NumTrials.Extra] ;
        TTrialCyc.LCyc = [LCycles.Currentsum; LCycles.Lastsum; LCycles.Extrasum];
        TTrialCyc.RCyc = [RCycles.Currentsum; RCycles.Lastsum; RCycles.Extrasum];
        TTrialCyc.FullCyc = [NumTrials.CurrentFullCyc; NumTrials.LastFullCyc; NumTrials.ExtraFullCyc];
    elseif strcmp(AddLast,'No') == 1
        TTrialCyc.Conditions = Tstd.Conditions;
        TTrialCyc.NumTrials = [NumTrials.Current; NumTrials.AFO; NumTrials.Extra] ;
        TTrialCyc.LCyc = [LCycles.Currentsum; LCycles.AFOsum; LCycles.Extrasum];
        TTrialCyc.RCyc = [RCycles.Currentsum; RCycles.AFOsum; RCycles.Extrasum];
        TTrialCyc.FullCyc = [NumTrials.CurrentFullCyc; NumTrials.AFOFullCyc; NumTrials.ExtraFullCyc];
    elseif strcmp(AddExtra,'No') == 1
        TTrialCyc.Conditions = Tstd.Conditions;
        TTrialCyc.NumTrials = [NumTrials.Current; NumTrials.AFO; NumTrials.Last] ;
        TTrialCyc.LCyc = [LCycles.Currentsum; LCycles.AFOsum; LCycles.Lastsum];
        TTrialCyc.RCyc = [RCycles.Currentsum; RCycles.AFOsum; RCycles.Lastsum];
        TTrialCyc.FullCyc = [NumTrials.CurrentFullCyc; NumTrials.AFOFullCyc; NumTrials.LastFullCyc];
    end
elseif NumConds == 4 % 4 conditions
    TTrialCyc.Conditions = Tstd.Conditions;
    TTrialCyc.NumTrials = [NumTrials.Current; NumTrials.AFO; NumTrials.Last; NumTrials.Extra] ;
    TTrialCyc.LCyc = [LCycles.Currentsum; LCycles.AFOsum; LCycles.Lastsum; LCycles.Extrasum;]; 
    TTrialCyc.RCyc = [RCycles.Currentsum; RCycles.AFOsum; RCycles.Lastsum; RCycles.Extrasum];
    TTrialCyc.FullCyc = [NumTrials.CurrentFullCyc; NumTrials.AFOFullCyc; NumTrials.LastFullCyc; NumTrials.ExtraFullCyc]; 
end

%% Copy and rename WAAAG template excel sheet for simple renaming
InputFileName = fullfile(pwd, 'WAAAG - Template.xlsx');
OutputFileName = fullfile(pwd, 'WAAAG.xlsx');

copyfile(InputFileName, OutputFileName);

clearvars InputFileName OutputFileName

%% Output to WAAG excel spreadsheet
xlswrite('WAAAG.xlsx',{date},'WAAAG','M2');
xlswrite('WAAAG.xlsx',Tstd.Conditions,'WAAAG','A3');
i = 1;
% current Condition
if strcmp(Type.Current, 'Full Kinematics') == 1 
   RepName{i} = [KineData(1).RepTrial.Name]; 
    RepTrials(i, :) = [KineData(1).RepTrial.RepCycles(1), KineData(1).RepTrial.RepCycles(2)];
else
    RepName{i} = [' '];
    RepTrials(i, :) = [0, 0]; 
end
Types{i} = [Type.Current]; 
% AFO Condition
if strcmp(AddAFO, 'Yes') == 1 
    i = i + 1;
    if strcmp(Type.AFO, 'Full Kinematics') == 1
        RepName{i} = [AFOData(1).RepTrial.Name];
        RepTrials(i, :) = [AFOData(1).RepTrial.RepCycles(1), AFOData(1).RepTrial.RepCycles(2)];
    else
        RepName{i} = [' '];
        RepTrials(i, :) = [0, 0];
    end
    Types{i} = [Type.AFO]; 
end
    % Last Condition
if strcmp(AddLast, 'Yes') == 1 
    i = i + 1;
    if strcmp(Type.Last, 'Full Kinematics') == 1
        RepName{i} = [LastData(1).RepTrial.Name];
        RepTrials(i, :) = [LastData(1).RepTrial.RepCycles(1), LastData(1).RepTrial.RepCycles(2)];
    else
        RepName{i} = [' '];
        RepTrials(i, :) = [0, 0];
    end
    Types{i} = [Type.Last]; 
end
% Extra Condition
if strcmp(AddExtra, 'Yes') == 1 
    i = i + 1;
    if strcmp(Type.Extra, 'Full Kinematics') == 1
        RepName{i} = [ExtraData(1).RepTrial.Name];
        RepTrials(i, :) = [ExtraData(1).RepTrial.RepCycles(1), ExtraData(1).RepTrial.RepCycles(2)];
    else
        RepName{i} = [' '];
        RepTrials(i, :) = [0, 0];
    end
    Types{i} = [Type.Extra]; 
end

xlswrite('WAAAG.xlsx', RepName','WAAAG','B3'); % export Rep trial names
xlswrite('WAAAG.xlsx', RepTrials,'WAAAG','D3'); % export Rep Trial cycles
xlswrite('WAAAG.xlsx', Types','WAAAG','G3'); % export trial type
         

xlswrite('WAAAG.xlsx',age,'WAAAG','J3'); % export subject age

xlswrite('WAAAG.xlsx',table2cell(T),'WAAAG','A9'); % export averages
xlswrite('WAAAG.xlsx',table2cell(Tstd),'WAAAG','A28'); % export STDs
xlswrite('WAAAG.xlsx',table2cell(TTrialCyc),'WAAAG','J18');

%% output average rotations
if strcmp(Type.Current,'Full Kinematics') == 1
    j = 1;
    for i = 1:2:9
        Rots(j,1) = ave_rotations(i).data;
        j = j+1;
    end
    j = 1;
    for i = 2:2:10
        Rots(j,2) = ave_rotations(i).data;
        j = j+1;
    end
    xlswrite('WAAAG.xlsx',Name.Curr,'WAAAG','B17');
    xlswrite('WAAAG.xlsx',Rots,'WAAAG','B19');
end
    
if strcmp(AddAFO,'Yes') == 1 % output AFO average rotations if present
    if strcmp(Type.AFO,'Full Kinematics') == 1
        j = 1;
        for i = 1:2:9
            AFO_Rots(j,1) = AFO_ave_rotations(i).data;
            j = j+1;
        end
        j = 1;
        for i = 2:2:10
            AFO_Rots(j,2) = AFO_ave_rotations(i).data;
            j = j+1;
        end
        xlswrite('WAAAG.xlsx',Name.AFO,'WAAAG','D17');
        xlswrite('WAAAG.xlsx',AFO_Rots,'WAAAG','D19');
    end
end

if strcmp(AddLast,'Yes') == 1 % output Last average rotations if present
    if strcmp(Type.Last,'Full Kinematics') == 1
        j = 1;
        for i = 1:2:9
            Last_Rots(j,1) = Last_ave_rotations(i).data;
            j = j+1;
        end
        j = 1;
        for i = 2:2:10
            Last_Rots(j,2) = Last_ave_rotations(i).data;
            j = j+1;
        end
        xlswrite('WAAAG.xlsx',Name.Last,'WAAAG','F17');
        xlswrite('WAAAG.xlsx',Last_Rots,'WAAAG','F19');
    end
end

if strcmp(AddExtra,'Yes') == 1 % output Extra average rotations if present
    if strcmp(Type.Extra,'Full Kinematics') == 1
        j = 1;
        for i = 1:2:9
            Extra_Rots(j,1) = Extra_ave_rotations(i).data;
            j = j+1;
        end
        j = 1;
        for i = 2:2:10
            Extra_Rots(j,2) = Extra_ave_rotations(i).data;
            j = j+1;
        end
        xlswrite('WAAAG.xlsx',Extra_Rots,'WAAAG','H19');
        xlswrite('WAAAG.xlsx',Name.Extra,'WAAAG','H17');
    end
end

%% Create tables of trial and cycle rankings - Current Condition
if strcmp(Type.Current, 'Full Kinematics') == 1
    % Trial Rankings
    TrialCount = 1:length(KineData);
    for i = 1:length(KineData)
        CurrFileNames{i} = KineData(i).FileName;
    end
    TrialRanks = table(TrialCount', CurrFileNames', KineData(1).RepTrial.TrialVariance.Measures(:,1),...
        KineData(1).RepTrial.TrialVariance.Measures(:,2), KineData(1).RepTrial.TrialVariance.Average, KineData(1).RepTrial.Rank');
    xlswrite('WAAAG.xlsx',table2cell(TrialRanks),'Curr_Trial_Cyc','A5'); %export trial ranks
    
    % Cycle Rankings for Current Condition
    % For 1st ranked trial
    Trial = KineData(1).RepTrial.Num(1);
    LeftLength = length(KineData(1).RepTrial.CycleRanks(Trial).left);
    RightLength = length(KineData(1).RepTrial.CycleRanks(Trial).right);
    while LeftLength < RightLength
        KineData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
        LeftLength = length(KineData(1).RepTrial.CycleRanks(Trial).left);
    end
    while RightLength < LeftLength
        KineData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
        RightLength = length(KineData(1).RepTrial.CycleRanks(Trial).right);
    end
    NumCycles = max(max(KineData(1).RepTrial.CycleRanks(Trial).left, KineData(1).RepTrial.CycleRanks(Trial).right));
    GCs = 1:NumCycles;
    CycleRanks = table(GCs', KineData(1).RepTrial.CycleRanks(Trial).left', KineData(1).RepTrial.CycleRanks(Trial).right');
    xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Curr_Trial_Cyc','G5'); %export trial ranks
    xlswrite('WAAAG.xlsx',{CurrFileNames{Trial}},'Curr_Trial_Cyc','H3'); % export #1 trial name
    clearvars GCs CycleRanks LeftLength RightLength
    
    if length(KineData) > 1
        % For 2nd ranked trial
        Trial = KineData(1).RepTrial.Num(2);
        LeftLength = length(KineData(1).RepTrial.CycleRanks(Trial).left);
        RightLength = length(KineData(1).RepTrial.CycleRanks(Trial).right);
        while LeftLength < RightLength
            KineData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
            LeftLength = length(KineData(1).RepTrial.CycleRanks(Trial).left);
        end
        while RightLength < LeftLength
            KineData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
            RightLength = length(KineData(1).RepTrial.CycleRanks(Trial).right);
        end
        NumCycles = max(max(KineData(1).RepTrial.CycleRanks(Trial).left, KineData(1).RepTrial.CycleRanks(Trial).right));
        GCs = 1:NumCycles;
        
        CycleRanks = table(GCs', KineData(1).RepTrial.CycleRanks(Trial).left', KineData(1).RepTrial.CycleRanks(Trial).right');
        xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Curr_Trial_Cyc','J5'); %export trial ranks
        xlswrite('WAAAG.xlsx',{CurrFileNames{Trial}},'Curr_Trial_Cyc','K3'); % export #2 trial name
        clearvars GCs CycleRanks LeftLength RightLength
        
        if length(KineData) > 2
            % For 3rd ranked trial
            Trial = KineData(1).RepTrial.Num(3);
            LeftLength = length(KineData(1).RepTrial.CycleRanks(Trial).left);
            RightLength = length(KineData(1).RepTrial.CycleRanks(Trial).right);
            while LeftLength < RightLength
                KineData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
                LeftLength = length(KineData(1).RepTrial.CycleRanks(Trial).left);
            end
            while RightLength < LeftLength
                KineData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
                RightLength = length(KineData(1).RepTrial.CycleRanks(Trial).right);
            end
            NumCycles = max(max(KineData(1).RepTrial.CycleRanks(Trial).left, KineData(1).RepTrial.CycleRanks(Trial).right));
            GCs = 1:NumCycles;
            CycleRanks = table(GCs', KineData(1).RepTrial.CycleRanks(Trial).left', KineData(1).RepTrial.CycleRanks(Trial).right');
            xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Curr_Trial_Cyc','M5'); %export trial ranks
            xlswrite('WAAAG.xlsx',{CurrFileNames{Trial}},'Curr_Trial_Cyc','N3'); % export #3 trial name
            clearvars GCs CycleRanks LeftLength RightLength
        end
    end
    
    clearvars TrialCount TrialRanks DotSize
end
%% Create tables of trial and cycle rankings - AFO Condition
if strcmp(AddAFO, 'Yes') == 1 % for AFO condition
    if strcmp(Type.AFO,'Full Kinematics') == 1
        % Trial Rankings
        TrialCount = 1:length(AFOData);
        for i = 1:length(AFOData)
            AFOFileNames{i} = AFOData(i).FileName;
        end
        TrialRanks = table(TrialCount', AFOFileNames', AFOData(1).RepTrial.TrialVariance.Measures(:,1),...
            AFOData(1).RepTrial.TrialVariance.Measures(:,2), AFOData(1).RepTrial.TrialVariance.Average, AFOData(1).RepTrial.Rank');
        xlswrite('WAAAG.xlsx',table2cell(TrialRanks),'AFO_Trial_Cyc','A5'); %export trial ranks
        
        % Cycle Rankings for Current Condition
        % For 1st ranked trial
        Trial = AFOData(1).RepTrial.Num(1);
        LeftLength = length(AFOData(1).RepTrial.CycleRanks(Trial).left);
        RightLength = length(AFOData(1).RepTrial.CycleRanks(Trial).right);
        while LeftLength < RightLength
            AFOData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
            LeftLength = length(AFOData(1).RepTrial.CycleRanks(Trial).left);
        end
        while RightLength < LeftLength
            AFOData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
            RightLength = length(AFOData(1).RepTrial.CycleRanks(Trial).right);
        end
        NumCycles = max(max(AFOData(1).RepTrial.CycleRanks(Trial).left, AFOData(1).RepTrial.CycleRanks(Trial).right));
        GCs = 1:NumCycles;
        CycleRanks = table(GCs', AFOData(1).RepTrial.CycleRanks(Trial).left', AFOData(1).RepTrial.CycleRanks(Trial).right');
        xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'AFO_Trial_Cyc','G5'); %export trial ranks
        xlswrite('WAAAG.xlsx',{AFOFileNames{Trial}},'AFO_Trial_Cyc','H3'); % export #1 trial name
        clearvars GCs CycleRanks LeftLength RightLength
        
        if length(AFOData) > 1
            % For 2nd ranked trial
            Trial = AFOData(1).RepTrial.Num(2);
            LeftLength = length(AFOData(1).RepTrial.CycleRanks(Trial).left);
            RightLength = length(AFOData(1).RepTrial.CycleRanks(Trial).right);
            while LeftLength < RightLength
                AFOData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
                LeftLength = length(AFOData(1).RepTrial.CycleRanks(Trial).left);
            end
            while RightLength < LeftLength
                AFOData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
                RightLength = length(AFOData(1).RepTrial.CycleRanks(Trial).right);
            end
            NumCycles = max(max(AFOData(1).RepTrial.CycleRanks(Trial).left, AFOData(1).RepTrial.CycleRanks(Trial).right));
            GCs = 1:NumCycles;
            
            CycleRanks = table(GCs', AFOData(1).RepTrial.CycleRanks(Trial).left', AFOData(1).RepTrial.CycleRanks(Trial).right');
            xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'AFO_Trial_Cyc','J5'); %export trial ranks
            xlswrite('WAAAG.xlsx',{AFOFileNames{Trial}},'AFO_Trial_Cyc','K3'); % export #2 trial name
            clearvars GCs CycleRanks LeftLength RightLength
            
            if length(AFOData) > 2
                % For 3rd ranked trial
                Trial = AFOData(1).RepTrial.Num(3);
                LeftLength = length(AFOData(1).RepTrial.CycleRanks(Trial).left);
                RightLength = length(AFOData(1).RepTrial.CycleRanks(Trial).right);
                while LeftLength < RightLength
                    AFOData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
                    LeftLength = length(AFOData(1).RepTrial.CycleRanks(Trial).left);
                end
                while RightLength < LeftLength
                    AFOData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
                    RightLength = length(AFOData(1).RepTrial.CycleRanks(Trial).right);
                end
                NumCycles = max(max(AFOData(1).RepTrial.CycleRanks(Trial).left, AFOData(1).RepTrial.CycleRanks(Trial).right));
                GCs = 1:NumCycles;
                CycleRanks = table(GCs', AFOData(1).RepTrial.CycleRanks(Trial).left', AFOData(1).RepTrial.CycleRanks(Trial).right');
                xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'AFO_Trial_Cyc','M5'); %export trial ranks
                xlswrite('WAAAG.xlsx',{AFOFileNames{Trial}},'AFO_Trial_Cyc','N3'); % export #3 trial name
                clearvars GCs CycleRanks LeftLength RightLength
            end
        end
    end
end
clearvars TrialCount TrialRanks

%% Create tables of trial and cycle rankings - Last Condition
if strcmp(AddLast, 'Yes') == 1 % for Last condition
    if strcmp(Type.Last,'Full Kinematics') == 1
        % Trial Rankings
        TrialCount = 1:length(LastData);
        for i = 1:length(LastData)
            LastFileNames{i} = LastData(i).FileName;
        end
        TrialRanks = table(TrialCount', LastFileNames', LastData(1).RepTrial.TrialVariance.Measures(:,1),...
            LastData(1).RepTrial.TrialVariance.Measures(:,2), LastData(1).RepTrial.TrialVariance.Average, LastData(1).RepTrial.Rank');
        xlswrite('WAAAG.xlsx',table2cell(TrialRanks),'Last_Trial_Cyc','A5'); %export trial ranks
        
        % Cycle Rankings for Current Condition
        % For 1st ranked trial
        Trial = LastData(1).RepTrial.Num(1);
        LeftLength = length(LastData(1).RepTrial.CycleRanks(Trial).left);
        RightLength = length(LastData(1).RepTrial.CycleRanks(Trial).right);
        while LeftLength < RightLength
            LastData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
            LeftLength = length(LastData(1).RepTrial.CycleRanks(Trial).left);
        end
        while RightLength < LeftLength
            LastData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
            RightLength = length(LastData(1).RepTrial.CycleRanks(Trial).right);
        end
        NumCycles = max(max(LastData(1).RepTrial.CycleRanks(Trial).left, LastData(1).RepTrial.CycleRanks(Trial).right));
        GCs = 1:NumCycles;
        CycleRanks = table(GCs', LastData(1).RepTrial.CycleRanks(Trial).left', LastData(1).RepTrial.CycleRanks(Trial).right');
        xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Last_Trial_Cyc','G5'); %export trial ranks
        xlswrite('WAAAG.xlsx',{LastFileNames{Trial}},'Last_Trial_Cyc','H3'); % export #1 trial name
        clearvars GCs CycleRanks LeftLength RightLength
        
        if length(LastData) > 1
            % For 2nd ranked trial
            Trial = LastData(1).RepTrial.Num(2);
            LeftLength = length(LastData(1).RepTrial.CycleRanks(Trial).left);
            RightLength = length(LastData(1).RepTrial.CycleRanks(Trial).right);
            while LeftLength < RightLength
                LastData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
                LeftLength = length(LastData(1).RepTrial.CycleRanks(Trial).left);
            end
            while RightLength < LeftLength
                LastData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
                RightLength = length(LastData(1).RepTrial.CycleRanks(Trial).right);
            end
            NumCycles = max(max(LastData(1).RepTrial.CycleRanks(Trial).left, LastData(1).RepTrial.CycleRanks(Trial).right));
            GCs = 1:NumCycles;
            
            CycleRanks = table(GCs', LastData(1).RepTrial.CycleRanks(Trial).left', LastData(1).RepTrial.CycleRanks(Trial).right');
            xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Last_Trial_Cyc','J5'); %export trial ranks
            xlswrite('WAAAG.xlsx',{LastFileNames{Trial}},'Last_Trial_Cyc','K3'); % export #2 trial name
            clearvars GCs CycleRanks LeftLength RightLength
            
            if length(LastData) > 2
                % For 3rd ranked trial
                Trial = LastData(1).RepTrial.Num(3);
                LeftLength = length(LastData(1).RepTrial.CycleRanks(Trial).left);
                RightLength = length(LastData(1).RepTrial.CycleRanks(Trial).right);
                while LeftLength < RightLength
                    LastData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
                    LeftLength = length(LastData(1).RepTrial.CycleRanks(Trial).left);
                end
                while RightLength < LeftLength
                    LastData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
                    RightLength = length(LastData(1).RepTrial.CycleRanks(Trial).right);
                end
                NumCycles = max(max(LastData(1).RepTrial.CycleRanks(Trial).left, LastData(1).RepTrial.CycleRanks(Trial).right));
                GCs = 1:NumCycles;
                CycleRanks = table(GCs', LastData(1).RepTrial.CycleRanks(Trial).left', LastData(1).RepTrial.CycleRanks(Trial).right');
                xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Last_Trial_Cyc','M5'); %export trial ranks
                xlswrite('WAAAG.xlsx',{LastFileNames{Trial}},'Last_Trial_Cyc','N3'); % export #3 trial name
                clearvars GCs CycleRanks LeftLength RightLength
            end
        end
    end
end
clearvars TrialCount TrialRanks

%% Create tables of trial and cycle rankings - Extra Condition
if strcmp(AddExtra, 'Yes') == 1 % for Extra condition
    if strcmp(Type.Extra,'Full Kinematics') == 1
        % Trial Rankings
        TrialCount = 1:length(ExtraData);
        for i = 1:length(ExtraData)
            ExtraFileNames{i} = ExtraData(i).FileName;
        end
        TrialRanks = table(TrialCount', ExtraFileNames', ExtraData(1).RepTrial.TrialVariance.Measures(:,1),...
            ExtraData(1).RepTrial.TrialVariance.Measures(:,2), ExtraData(1).RepTrial.TrialVariance.Average, ExtraData(1).RepTrial.Rank');
        xlswrite('WAAAG.xlsx',table2cell(TrialRanks),'Extra_Trial_Cyc','A5'); %export trial ranks
        
        % Cycle Rankings for Current Condition
        % For 1st ranked trial
        Trial = ExtraData(1).RepTrial.Num(1);
        LeftLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).left);
        RightLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).right);
        while LeftLength < RightLength
            ExtraData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
            LeftLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).left);
        end
        while RightLength < LeftLength
            ExtraData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
            RightLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).right);
        end
        NumCycles = max(max(ExtraData(1).RepTrial.CycleRanks(Trial).left, ExtraData(1).RepTrial.CycleRanks(Trial).right));
        GCs = 1:NumCycles;
        CycleRanks = table(GCs', ExtraData(1).RepTrial.CycleRanks(Trial).left', ExtraData(1).RepTrial.CycleRanks(Trial).right');
        xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Extra_Trial_Cyc','G5'); %export trial ranks
        xlswrite('WAAAG.xlsx',{ExtraFileNames{Trial}},'Extra_Trial_Cyc','H3'); % export #1 trial name
        clearvars GCs CycleRanks LeftLength RightLength
        
        if length(ExtraData) > 1
            % For 2nd ranked trial
            Trial = ExtraData(1).RepTrial.Num(2);
            LeftLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).left);
            RightLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).right);
            while LeftLength < RightLength
                ExtraData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
                LeftLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).left);
            end
            while RightLength < LeftLength
                ExtraData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
                RightLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).right);
            end
            NumCycles = max(max(ExtraData(1).RepTrial.CycleRanks(Trial).left, ExtraData(1).RepTrial.CycleRanks(Trial).right));
            GCs = 1:NumCycles;
            
            CycleRanks = table(GCs', ExtraData(1).RepTrial.CycleRanks(Trial).left', ExtraData(1).RepTrial.CycleRanks(Trial).right');
            xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Extra_Trial_Cyc','J5'); %export trial ranks
            xlswrite('WAAAG.xlsx',{ExtraFileNames{Trial}},'Extra_Trial_Cyc','K3'); % export #2 trial name
            clearvars GCs CycleRanks LeftLength RightLength
            
            if length(ExtraData) > 2
                % For 3rd ranked trial
                Trial = ExtraData(1).RepTrial.Num(3);
                LeftLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).left);
                RightLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).right);
                while LeftLength < RightLength
                    ExtraData(1).RepTrial.CycleRanks(Trial).left(LeftLength + 1) = 0;
                    LeftLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).left);
                end
                while RightLength < LeftLength
                    ExtraData(1).RepTrial.CycleRanks(Trial).right(RightLength + 1) = 0;
                    RightLength = length(ExtraData(1).RepTrial.CycleRanks(Trial).right);
                end
                NumCycles = max(max(ExtraData(1).RepTrial.CycleRanks(Trial).left, ExtraData(1).RepTrial.CycleRanks(Trial).right));
                GCs = 1:NumCycles;
                CycleRanks = table(GCs', ExtraData(1).RepTrial.CycleRanks(Trial).left', ExtraData(1).RepTrial.CycleRanks(Trial).right');
                xlswrite('WAAAG.xlsx',table2cell(CycleRanks),'Extra_Trial_Cyc','M5'); %export trial ranks
                xlswrite('WAAAG.xlsx',{ExtraFileNames{Trial}},'Extra_Trial_Cyc','N3'); % export #3 trial name
                clearvars GCs CycleRanks LeftLength RightLength
            end
        end
    end
end
clearvars TrialCount TrialRanks

%% Clear erroneous variables to clean up workspace
clearvars i Blue RedGreen prompt X1 X2 Version ax Columns XTicks
save('WAAAG_data.mat');
close all;
