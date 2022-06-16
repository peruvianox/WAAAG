function [resid] = Residual(mean_var, variable)

% -------------------------------------------------------------------------
% RESIDUAL
% -------------------------------------------------------------------------

% Author: Kate Worster
% Date: December 14, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This is function calculates residual between the mean and data.
%
% Input:        mean_var                 mean of variable
%               variable                 variable to compare to mean
%
% OUTPUTS:  
%               resid                    residual for each gait cycle
%                                        (column)
%
%


%% Residual
  
for numGC = 1:length(variable(1,:))  % number of gait cycles
    for t = 1:length(variable(:,1))  % number of time points
        % delta = mean-variable for each gait cycle
        delta(t,numGC) = abs(mean_var(t,1) - variable(t,numGC));
    end
end

% Residual for each gait cycle
resid = sum(delta)/(length(mean_var(:,1))-1);
    

