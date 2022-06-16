function [resamp_var] = Resample(var_in)

% -------------------------------------------------------------------------
% GDI RESAMPLE
% -------------------------------------------------------------------------

% *Slightly edited by Ricky Pimentel, 11,21,16 in order to interpolate
% variables down to 101 points from any length of input

% Original Author: Kate Worster 
% Date: September 11, 2009
% For use at the Center for Gait and Movement Laboratory at the Children's
% Hospital in Denver, CO.
%
% Description:	This is function up-samples the input variable (51pts) to
%               101 data points. 
%
% Input:        var_in           variable to be resampled (nX1 array)
%               new_samp_size    size of data var_in is to be resampled
%
%
% OUTPUTS:  
%               resamp_var       resampled variable
%
%% DEBUG LOOP
% 
% close all
% clear all
% clc
% 
% x = 1:51;
% var_in(:,1) = sin(x);


%% Resample Data
    
N = length(var_in);
var_x = 1:N; % number of data points to interpolate from
var_y = var_in(:,1); % data points to interpolate from
var_xinterp = 1:N/101:N; % new sample rate (interpolated curve's x points)
resamp_var(:,1) = interp1(var_x, var_y, var_xinterp, 'pchip'); % interpolated curve's y points


%% DEBUG PLOTS
% subplot(3,1,1)
% t = 1:length(var_in);
% plot(t, var_in(t),'k')
% title('Original Data')
% 
% subplot(3,1,2)
% t2 = 1:length(resamp_var);
% plot(t2, resamp_var(t2),'b')
% title('Interpolated')


