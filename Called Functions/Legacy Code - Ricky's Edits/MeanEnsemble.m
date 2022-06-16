function [mean_var, std_dev] = MeanEnsemble(variable)

% -------------------------------------------------------------------------
% MEAN ENSEMBLE
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: July 20, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This is function calculates the mean and standard deviation
%               for the input variable, for all gait cycles.
%
% Input:        variable        data from which mean ensemble is calculated
%
% OUTPUTS:  
%               mean_var        mean of input variable
%               std_dev         standard deviation of variable
%
% UPDATE: Corrected equation 5/1/12

%% Calculates the mean for input variable

for t = 1:size(variable,1)
    % mean of variable for all gait cycles and both sides
    mean_var(t,1) = sum(variable(t,:))/length(variable(1,:));
end

%% Calculates the standard deviation for variable
for t = 1:size(variable,1)
    for numGC = 1:length(variable(1,:))
        % stdev = square root of (1/n-1)*sum(var-mean)^2
        std_dev(t,numGC) = sqrt( (1/(length(variable)-1))*(sum(variable(t,numGC)-mean_var(t,1)).^2) );

    end
end
