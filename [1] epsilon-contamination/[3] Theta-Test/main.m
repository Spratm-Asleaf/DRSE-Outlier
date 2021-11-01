%{
Online supplementary materials of the paper titled 
"Distributionally Robust State Estimation for Linear Systems Subject to Uncertainty and Outlier"
Authored By: Shixiong Wang, and Zhisheng Ye
From the Department of Industrial Systems Engineering and Management, National University of Singapore

@Author: Shixiong Wang
@Date: 10 May 2021, 8 Sep 2021
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE-Outlier

Acknowledgements:

"WKF.m", "FrankWolfe.m" and "tauKF.m" are adapted from sources at
    "https://github.com/sorooshafiee/WKF".

For "tauKF.m", we also acknowledge M. Zorz.
    "https://mathworks.com/matlabcentral/fileexchange/54308-robust-kalman-filtering-package?focused=5751770&tab=function"
%}

clc;
clear all;
close all;


% Problem setting
% n = 2; % dimension of state vector
% m = 1; % dimension of measurement
sys.F = [0.9802, 0.0196; 0, 0.9802];
sys.G = [1, 0; 0, 1];
sys.H = [1, -1];
sys.D = 1;
sys.Q = [1.9608, 0.0195; 0.0195, 1.9605];
sys.R = 1;

% Simulate 1000 steps, length of an episode
T = 10000;

% RMSE v.s. theta
expected_theta = 0:0.01:0.1;
theta_len = length(expected_theta);

RMSE_MKF_Array = zeros(theta_len, 1);

for theta_index = 1:theta_len
    disp(['Loading: ' num2str(theta_index*100/theta_len) '%']);
    
    real_theta = expected_theta(theta_index);
    
    % Init conditions
    x_0 = randn(2,1);
    V_0 = eye(2);

    % Parameters
    isTimeVariant = true;
    coeff_alpha = 1;
    epsilon = 0.05;
    
    %% Simulate data
    [x, y, real_y, y0, real_F] = GenerateData(sys, x_0, T, coeff_alpha, epsilon, isTimeVariant);
    
    %% Moment filter
    tic
    theta = real_theta;   % 0.025
    theta_1 = 1 - theta;
    theta_2 = 1 + theta;
    k = 1.4;      % use 2 or 1.4, make no big difference; try 1.4 if you do not believe; 

    xhat_MKF = MKF(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0,k,epsilon,theta_2);
    time_MKF = toc/T;
    % disp(['Avg Time of Moment: ' num2str(time_MKF)]);
    err_MKF = sum((x-xhat_MKF).^2,1);
    RMSE_MKF = sqrt(mean(err_MKF));
    % disp(['---------------------RMSE of MKF: ' num2str(RMSE_MKF)]);
    RMSE_MKF_Array(theta_index) = RMSE_MKF;
end

% theta_2 := theta + 1
plot(expected_theta + 1, RMSE_MKF_Array, 'r',...
    'linewidth', 2);
set(gca, 'FontSize', 18);
ylabel('RMSE','FontSize', 20, 'Interpreter', 'latex');
xlabel('$\theta_2$','FontSize', 20, 'Interpreter', 'latex');

