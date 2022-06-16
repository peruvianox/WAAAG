function[RepTrial, GDI] =  DetermineRepTrialCyc (files2get)
% % -------------------------------------------------------------------------
% Determines Representative Trial and Cycles
% -------------------------------------------------------------------------
% Author: Ricky Pimentel
% Date: January 5, 2017
%
% Description: This function will determine the representative trial and
% cycle from a single trial or numerous trials. Must select C3D files as inputs.
%
% Input:       files2get         name of file(s) desired to select the representative trial and cycle
%                                           * is called within, thus does not need to be specified
%
% Output:     RepTrial          Representative Trial Name and Number (of selected trials)
%                   RepCycles       The representative cycle for the left and right sides
%
% % DEBUG LOOP %%%%%%%%%
% close all
% clear
% clc
% %%% END DEBUG LOOP %%%%%

warning off; 

%% Select data for IMPORT
if exist('files2get','var') == 0  % test to see if data exists, if not -> load with uigetfile below
    [files2get, ~] = uigetfile('*.c3d', 'Please select C3D file(s) to be analyzed', 'MultiSelect', 'on');
end

%% Batch Processes All Files
% files2get = files2get';
if iscell(files2get) == 1
    % if multiple files are selected
    for n = 1:size(files2get, 2)
        % converts name of each file from cell structure to char
        FullFileName = char(cell2mat( files2get(n)));
        % Analyzes current C3D file
        [variance_ratios, residuals, ~, ~, ~, ~, ~] = GDI_AnalyTrial(FullFileName);
        % builds structure of current trial's name
        c3d_files(n).name = FullFileName;
        % structure of each trial's temporal spatial variables
        %c3d_files(n).tempspat = TempSpat;
        % structure of trial's variance ratios (12x1)
        c3d_files(n).variance_ratios = variance_ratios.trial;
        % ordered structure of trial's representative gait cycles for all left & right gait cycles
        c3d_files(n).residuals = residuals;
    end
else % if only one file is selected
    FullFileName = files2get;
    % Analyzes current C3D file
    [variance_ratios, residuals, ~, ~, ~, ~, ~] = GDI_AnalyTrial(FullFileName);
    % builds structure of current trial's name
    c3d_files.name = FullFileName;
    % structure of each trial's temporal spatial variables
    %c3d_files.tempspat = TempSpat;
    % structure of trial's variance ratios (12x1)
    c3d_files.variance_ratios = variance_ratios.trial;
    % ordered structure of trial's representative gait cycles for all left & right gait cycles
    c3d_files.residuals = residuals;
end

%% Determines Representative Trial
% Determines the temporal spatial residuals
%[sum_resids, c3d_files] = GAMS_TempSpatResids(c3d_files);
for t = 1:length(c3d_files)
    % sorts each trial's variance ratios
    trial_variance(:,t) = c3d_files(t).variance_ratios;
end
% ranks the representative trails, in ascending order, left = 1, right = 2
for i = 1:2
    if i == 1
        for j = 1:length(trial_variance)
            t_v(j) = trial_variance(j).left;
        end
        [~, trial] = sort(t_v);
        Representative.trial_order.left = trial;
    else
        for j = 1:length(trial_variance)
            t_v(j) = trial_variance(j).right;
        end
        [~, trial] = sort(t_v);
        Representative.trial_order.right = trial;
    end
end

% trial ranking by averaging residuals of left and right sides 
for i = 1:length(trial_variance)
TrialVariance.Measures(i,1) = trial_variance(i).left;
TrialVariance.Measures(i,2) = trial_variance(i).right;
end
TrialVariance.Average = mean(TrialVariance.Measures,2); 
% [~, ind] = min(TrialVariance.Average); 
if length(c3d_files) == 1
    ind = 1; 
    Rank = 1; 
else 
    [~, ind] = sort(TrialVariance.Average);
    for i = 1:length(ind)
        a = ind(i);
        Rank(a) = i;
    end
end

% % replicates matrix of trial names from structures format
% rep_trials = (repmat({c3d_files.name}, 1,1));
% 
% % sorts names of trials by rank order
% for t = 1:length(rep_trials)
%     ordered_rep_trials(t) = rep_trials(trial(t));
% end

%% Save rep trial for output
% RepTrial.Name = ordered_rep_trials(1);
% RepTrialNum = [Representative.trial_order.left(1), Representative.trial_order.right(1)];
if length(c3d_files) == 1
    RepTrial.Name = FullFileName; 
else
    RepTrial.Name = files2get{ind(1)};
end

RepTrial.Num = ind;
RepTrial.Rank = Rank; 
RepTrial.TrialVariance = TrialVariance;

%% Define start times for L and R Reps
% determine minimum # of cycles for R and L sides
% Num_GCs = min(length(residuals.all_GCs.left), length(residuals.all_GCs.right)); 
% if Num_GCs > 3 % if more than 3 cycles for each side, only use the top 3 choices to speed up code
%     cyc = 3;
% else
%     cyc = Num_GCs;
% end
% L_Rep = residuals.all_GCs.left(1:cyc);
% R_Rep = residuals.all_GCs.right(1:cyc);
RepCycles(1) = c3d_files(RepTrial.Num(1)).residuals.all_GCs.left(1); 
RepCycles(2) = c3d_files(RepTrial.Num(1)).residuals.all_GCs.right(1); 

% % Determine gait cycle times
% [~, ~, ~, ~, ~, ~, GaitEvents] = GDI_AnalyTrial(RepTrial.Name);
% [~, LGaitCycles, RGaitCycles] = GaitEventTimes(GaitEvents);
% LCycleStart = LGaitCycles(1,1,1);
% RCycleStart = RGaitCycles(1,1,1);

for i = 1:length(c3d_files) % define the rankings for all the gait cycles
RepTrial.CycleRanks(i).left = c3d_files(i).residuals.all_GCs.left;
RepTrial.CycleRanks(i).right = c3d_files(i).residuals.all_GCs.right;
end

% %% Determine L & R Rep Cycles
% if exist('Hemi','var') == 0 % if Hemi is not defined, prompt input
%     prompt = ('Is the subject hemiplegic?');
%     Hemi = questdlg(prompt, 'Hemi','Left','Right','No','No');
% end
% 
% if strcmp(Hemi, 'No') == 1 % If no Hemi find which side has the lowest GDI
%     GDI = GDI_Calculator_Fast(L_Rep(1), R_Rep(1), RepTrial.Name,0);
% elseif strcmp(Hemi, 'Left') == 1 % if left Hemi, temporarily define GDI with left side lower
%     GDI = [1 2];
% elseif strcmp(Hemi, 'Right') == 1 % if right Hemi, temporarily define GDI with right side lower
%     GDI = [2 1];
% end
% 
% if GDI(1) < GDI(2) % If left GDI is less than the right use left representative
%     RepCycles(1) = L_Rep(1); % cycle with the least amount of residual is chosen as the left rep
%     if LCycleStart < RCycleStart % If the first cycle occurs on the left foot
%         if L_Rep(1) == 1 % if rep cycle is #1
%             RepCycles(2) = 1; % Right cycle must be rep bc there are no other adjacent cycles
%         elseif L_Rep(1) == length(L_Rep) % if the left cycle is the last rep on that side
%             if length(R_Rep) < length(L_Rep) % if the right side is shorter than the right
%                 RepCycles(2) = L_Rep(1) - 1; % choose the right cycle directly before
%             else % if there is a right cycle after the final left
%                 ind1 = find(R_Rep == L_Rep(1));  % compare the reps before and after
%                 ind2 = find(R_Rep == L_Rep(1) - 1);
%                 RepCycles(2) = R_Rep(min(ind1, ind2)); % if there is no left cycle after the one on the left, the last left cycle is the only adjacent cycle
%             end
%         else % if the rep cycle is not the first or the last
%             ind1 = find(R_Rep == L_Rep(1));
%             ind2 = find(R_Rep == L_Rep(1) - 1);
%             RepCycles(2) = R_Rep(min(ind1, ind2)); % choose the left cycle ranked the lowest
%         end
%     else %If the first cycle occurs on the right foot
%         if L_Rep(1) == 1 % if left rep cycle is #1
%             ind1 = find(R_Rep == L_Rep(1));  % compare the reps before and after
%             ind2 = find(R_Rep == L_Rep(1) + 1);
%             RepCycles(2) = R_Rep(min(ind1, ind2));
%         elseif L_Rep(1) == length(L_Rep) % if the left cycle is the last rep on that side
%             if length(R_Rep) > length(L_Rep) % if there is a cycle on the right after the left
%                 ind1 = find(R_Rep == L_Rep(1));  % compare the reps before and after
%                 ind2 = find(R_Rep == L_Rep(1) + 1);
%                 RepCycles(2) = R_Rep(min(ind1, ind2)); % choose the left cycle ranked the lowest
%             else
%                 RepCycles(2) = L_Rep(1); % if there is no left cycle after the one on the left
%             end
%         else % if the rep cycle is not the first or the last
%             ind1 = find(R_Rep == L_Rep(1));
%             ind2 = find(R_Rep == L_Rep(1) + 1);
%             RepCycles(2) = R_Rep(min(ind1, ind2)); % choose the left cycle ranked the lowest
%         end
%     end
% else %If right GDI is lower than the left
%     RepCycles(2) = R_Rep(1); % cycle with the least amount of residual is chosen as the right rep
%     if LCycleStart < RCycleStart % If the first cycle occurs on the left foot
%         if R_Rep(1) == 1 % if rep cycle is #1
%             ind1 = find(L_Rep == R_Rep(1));
%             ind2 = find(L_Rep == R_Rep(1) +1);
%             RepCycles(1) = L_Rep(min(ind1, ind2)); % choose the left cycle ranked the lowest
%         elseif R_Rep(1) == length(R_Rep) % if the right cycle is the last rep on that side
%             if length(L_Rep) > length(R_Rep) % if there is a cycle on the left after the right
%                 ind1 = find(L_Rep == R_Rep(1));  % compare the reps before and after
%                 ind2 = find(L_Rep == R_Rep(1) + 1);
%                 RepCycles(1) = L_Rep(min(ind1, ind2)); % choose the left cycle ranked the lowest
%             else
%                 RepCycles(1) = R_Rep(1); % if there is no left cycle after the one on the left, the last left cycle is the only adjacent cycle
%             end
%         else % if the rep cycle is not the first or the last
%             ind1 = find(L_Rep == R_Rep(1));
%             ind2 = find(L_Rep == R_Rep(1) + 1);
%             RepCycles(1) = L_Rep(min(ind1, ind2)); % choose the left cycle ranked the lowest
%         end
%     else %If the first cycle occurs on the right foot
%         if R_Rep(1) == 1 % if rep cycle is #1
%             RepCycles(1) = 1; % the only adjacent left cycle is directly after
%         elseif R_Rep(1) == length(R_Rep) % if the right cycle is the last rep on that side
%             if length(L_Rep) == length(R_Rep) % if there is a cycle on the left after the right
%                 ind1 = find(L_Rep == R_Rep(1));  % compare the reps before and after
%                 ind2 = find(L_Rep == R_Rep(1) -1);
%                 RepCycles(1) = L_Rep(min(ind1, ind2)); % choose the left cycle ranked the lowest
%             else
%                 RepCycles(1) = R_Rep(1) - 1; % if there is no left cycle after the one on the left
%             end
%         else % if the rep cycle is not the first or the last
%             ind1 = find(L_Rep == R_Rep(1));
%             ind2 = find(L_Rep == R_Rep(1) -1);
%             RepCycles(1) = L_Rep(min(ind1, ind2)); % choose the left cycle ranked the lowest
%             if ind2 == 0
%                 RepCycles(1) = ind1;
%             end
%         end
%     end
% end
GDI = GDI_Calculator_Fast(RepCycles(1), RepCycles(2), FullFileName, 0); % re-calculate GDI with specific cycles

%%  save the representative cycles
RepTrial.RepCycles = RepCycles;

end