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

% RMSE v.s. epsilon
variable_epsilon = 0:0.01:0.5;
epsilon_len = length(variable_epsilon);
RMSE_TMKF_Array = zeros(epsilon_len,1);
RMSE_KF_Array = zeros(epsilon_len,1);
RMSE_Huber_Array = zeros(epsilon_len,1);
RMSE_tauKF_Array = zeros(epsilon_len,1);
RMSE_WKF_Array = zeros(epsilon_len,1);
RMSE_MKF_Array = zeros(epsilon_len,1);

for eps_index = 1:epsilon_len
    disp(['Loading: ' num2str(eps_index*100/epsilon_len) '%']);
    
    real_epsilon = variable_epsilon(eps_index);
    
    % Init conditions
    x_0 = randn(2,1);
    V_0 = eye(2);

    % Parameters
    isTimeVariant = true;
    coeff_alpha = 1;
    epsilon = 0.05;
    
    %% Simulate data
    [x, y, real_y, y0, real_F] = GenerateData(sys, x_0, T, coeff_alpha, real_epsilon, isTimeVariant);
    
    %% TMKF: Optimal Kalman filter (i.e., KF with the true model)
    % suffix "TM" is for "True Model"
    x_0 = [0;0];
    tic
    xhat_TMKF = TMKF(real_y,real_F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
    time_TMKF = toc/T;
%     disp(['Avg Time of TMKF: ' num2str(time_TMKF)]);
    err_TMKF = sum((x - xhat_TMKF).^2,1);
    RMSE_TMKF = sqrt(mean(err_TMKF));
%     disp(['---------------------RMSE of TMKF: ' num2str(RMSE_TMKF)]);
    RMSE_TMKF_Array(eps_index) = RMSE_TMKF;


    %% KF: Standard Kalman filter
    % The immediate line below is to test the validity of WKF.
    % Note that when the radius of uncertainty set equals to 0, ...
    % WKF degrades into the standard KF.
    % xhat = WKF(sys, 0, y, x_0, V_0); err_KF_0 = sum((x-xhat).^2,1);
    x_0 = [0;0];
    tic
    xhat_KF = KF(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
    time_KF = toc/T;
    % disp(['Avg Time of KF: ' num2str(time_KF)]);
    err_KF = sum((x - xhat_KF).^2,1);
    RSME_KF = sqrt(mean(err_KF));
    % disp(['---------------------RMSE of KF: ' num2str(RSME_KF)]);
    RMSE_KF_Array(eps_index) = RSME_KF;


    %% Huber Kalman filter
    x_0 = [0;0];
    k = 1.4;      % use 2 or 1.4, does not matter; try 1.4 if you do not believe; 
    tic
    xhat_Huber = Huber(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0,k,epsilon);
    time_Huber = toc/T;
    % disp(['Avg Time of Huber: ' num2str(time_Huber)]);
    err_Huber = sum((x - xhat_Huber).^2,1);
    RSME_Huber = sqrt(mean(err_Huber));
    % disp(['---------------------RMSE of Huber: ' num2str(RSME_Huber)]);
    RMSE_Huber_Array(eps_index) = RSME_Huber;
    
    %% tau filter (Note: when tau-divergence filter takes tau = 0, it gives the Kullback-Leibler divergence filter)
    tau = 0;
    c = 1.5e-4;
    y_delay = [y0, y(:,1:end-1)]';
    x_0 = [0;0];
    tic
    xhat_tauKF = tauKF(sys, c, tau, y_delay, x_0, V_0);
    time_tauKF = toc/T;
%     disp(['Avg Time of tauKF: ' num2str(time_tauKF)]);
    err_tauKF = sum((x - xhat_tauKF').^2,1);
    RMSE_tauKF = sqrt(mean(err_tauKF));
%     disp(['---------------------RMSE of tauKF: ' num2str(RMSE_tauKF)]);
    RMSE_tauKF_Array(eps_index) = RMSE_tauKF;


    %% Wasserstein filter
    x_0 = [0;0];
    tic
    rho = 0.1;
    [xhat_WKF, V_WKF, G_WKF, S_WKF] = WKF(sys, rho, y, x_0, V_0);
    time_WKF = toc/T;
%     disp(['Avg Time of WKF: ' num2str(time_WKF)]);
    err_WKF = sum((x - xhat_WKF).^2,1);
    RMSE_WKF = sqrt(mean(err_WKF));
%     disp(['---------------------RMSE of WKF: ' num2str(RMSE_WKF)]);
    RMSE_WKF_Array(eps_index) = RMSE_WKF;


    %% Moment filter
    tic
    theta = 0.02;   % 0.025
    theta_1 = 1 - theta;
    theta_2 = 1 + theta;
    xhat_MKF = MKF(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0,k,epsilon,theta_2);
    time_MKF = toc/T;
    % disp(['Avg Time of Moment: ' num2str(time_MKF)]);
    err_MKF = sum((x-xhat_MKF).^2,1);
    RMSE_MKF = sqrt(mean(err_MKF));
    % disp(['---------------------RMSE of MKF: ' num2str(RMSE_MKF)]);
    RMSE_MKF_Array(eps_index) = RMSE_MKF;
end

plot(variable_epsilon, RMSE_TMKF_Array, 'k',...
    variable_epsilon, RMSE_KF_Array, 'r',...
    variable_epsilon, RMSE_Huber_Array, 'm',...
    variable_epsilon, RMSE_tauKF_Array, 'g',...
    variable_epsilon, RMSE_WKF_Array, 'c',...
    variable_epsilon, RMSE_MKF_Array, 'b',...
    'linewidth', 2);
set(gca, 'FontSize', 18);
ylabel('RMSE','FontSize', 20, 'Interpreter', 'latex');
xlabel('$\epsilon_{real}$','FontSize', 20, 'Interpreter', 'latex');
leg1 = legend({'TMKF', 'KF', 'Huber', '$\tau$-KF', 'WKF', 'MKF'}, 'Interpreter', 'latex', 'Location', 'northwest');

set(leg1,'FontSize',15);

%% Save data
%save ONLY_UNCERTAINTY RMSE_TMKF_Array RMSE_KF_Array RMSE_Huber_Array RMSE_tauKF_Array RMSE_WKF_Array RMSE_MKF_Array variable_epsilon
