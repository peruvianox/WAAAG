function [sum_resids, c3d_files] = GAMS_TempSpatResids(c3d_files)
% -------------------------------------------------------------------------
% TEMPORAL SPATIAL RESIDUALS CALCULATION
% -------------------------------------------------------------------------
% Author: Kate Worster
% Date: December 23, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This function calculates and sums the residuals for all temporal 
%               spatial variables.
%               
%
% Input:        c3d_files           structure of each c3d file's data
%               
%
% Output:       sum_resids          matrix of summed residuals for TempSpat
%               c3d_files           added sum_resids to structure
%
%


%% Variance ratio for Temporal Spatial Structure

for f = 1:length(c3d_files)

    % matrix of all Cadence values for all trials
    all_Cadence(f) = c3d_files(f).tempspat(1).data;
    all_WalkSpeed(f) = c3d_files(f).tempspat(2).data;
    all_StrideLength(f) = c3d_files(f).tempspat(3).data;
    
    all_LStride_time(f) = c3d_files(f).tempspat(4).data;
    all_LStep_time(f) = c3d_files(f).tempspat(5).data;
    all_LStep_length(f) = c3d_files(f).tempspat(6).data;
    all_LOppFOff(f) = c3d_files(f).tempspat(7).data;
    all_LOppFOn(f) = c3d_files(f).tempspat(8).data;
    all_LFootOff(f) = c3d_files(f).tempspat(9).data;
    all_LIDS(f) = c3d_files(f).tempspat(10).data;
    all_LSS(f) = c3d_files(f).tempspat(11).data;
    all_LFDS(f) = c3d_files(f).tempspat(12).data;
    
    all_RStride_time(f) = c3d_files(f).tempspat(13).data;
    all_RStep_time(f) = c3d_files(f).tempspat(14).data;
    all_RStep_length(f) = c3d_files(f).tempspat(15).data;
    all_ROppFOff(f) = c3d_files(f).tempspat(16).data;
    all_ROppFOn(f) = c3d_files(f).tempspat(17).data;
    all_RFootOff(f) = c3d_files(f).tempspat(18).data;
    all_RIDS(f) = c3d_files(f).tempspat(19).data;
    all_RSS(f) = c3d_files(f).tempspat(20).data;
    all_RFDS(f) = c3d_files(f).tempspat(21).data;
end


mean_all_Cadence = mean(all_Cadence);
variable = all_Cadence;  mean_var = mean_all_Cadence;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_Cadence(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_WalkSpeed = mean(all_WalkSpeed);
variable = all_WalkSpeed;  mean_var = mean_all_WalkSpeed;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_WalkSpeed(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_StrideLength = mean(all_StrideLength);
variable = all_StrideLength;  mean_var = mean_all_StrideLength;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_StrideLength(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LStridet = mean(all_LStride_time);
variable = all_LStride_time;  mean_var = mean_all_LStridet;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LStridet(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LStept = mean(all_LStep_time);
variable = all_LStep_time;  mean_var = mean_all_LStept;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LStept(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LStepL = mean(all_LStep_length);
variable = all_LStep_length;  mean_var = mean_all_LStepL;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LStepL(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LOFOff = mean(all_LOppFOff);
variable = all_LOppFOff;  mean_var = mean_all_LOFOff;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LOFOff(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LOFOn = mean(all_LOppFOn);
variable = all_LOppFOn;  mean_var = mean_all_LOFOn;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LOFOn(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LFOff = mean(all_LFootOff);
variable = all_LFootOff;  mean_var = mean_all_LFOff;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LFOff(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LIDS = mean(all_LIDS);
variable = all_LIDS;  mean_var = mean_all_LIDS;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LIDS(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LSS = mean(all_LSS);
variable = all_LSS;  mean_var = mean_all_LSS;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LSS(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_LFDS = mean(all_LFDS);
variable = all_LFDS;  mean_var = mean_all_LFDS;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_LFDS(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_RStridet = mean(all_RStride_time);
variable = all_RStride_time;  mean_var = mean_all_RStridet;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_RStridet(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_RStept = mean(all_RStep_time);
variable = all_RStep_time;  mean_var = mean_all_RStept;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_RStept(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_RStepL = mean(all_RStep_length);
variable = all_RStep_length;  mean_var = mean_all_RStepL;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_RStepL(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_ROFOff = mean(all_ROppFOff);
variable = all_ROppFOff;  mean_var = mean_all_ROFOff;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_ROFOff(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_ROFOn = mean(all_ROppFOn);
variable = all_ROppFOn;  mean_var = mean_all_ROFOn;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_ROFOn(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_RFOff = mean(all_RFootOff);
variable = all_RFootOff;  mean_var = mean_all_RFOff;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_RFOff(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_RIDS = mean(all_RIDS);
variable = all_RIDS;  mean_var = mean_all_RIDS;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_RIDS(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_RSS = mean(all_RSS);
variable = all_RSS;  mean_var = mean_all_RSS;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_RSS(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

mean_all_RFDS = mean(all_RFDS);
variable = all_RFDS;  mean_var = mean_all_RFDS;
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        dev_RFDS(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end


%% Variance Ratios & Coeff. of Variation for All Kinematic & Temporal Spatial Data
% % Composite matrix of all residuals for the C3D files
residuals = [dev_Cadence; dev_WalkSpeed; dev_StrideLength;
             dev_LStridet; dev_LStept; dev_LStepL;
             dev_LOFOff; dev_LOFOn; dev_LFOff;
             dev_LIDS; dev_LSS; dev_LFDS;
             dev_RStridet; dev_RStept; dev_RStepL;
             dev_ROFOff; dev_ROFOn; dev_RFOff;
             dev_RIDS; dev_RSS; dev_RFDS];

% sum of residuals for each trial
for s = 1:length(residuals(1,:))
    sum_resids(1,s) = sum(residuals(:,s));
    c3d_files(s).residTS = sum_resids(1,s);
end



