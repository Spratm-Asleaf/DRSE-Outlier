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


%% Problem setting
% n = 2; % dimension of state vector
% m = 1; % dimension of measurement
sys.F = [0.9802, 0.0196; 0, 0.9802];
sys.G = [1, 0; 0, 1];
sys.H = [1, -1];
sys.D = 1;
sys.Q = [1.9608, 0.0195; 0.0195, 1.9605];
sys.R = 1;

%% Simulate 10000 steps, length of an episode
T = 10000;
% Init conditions
x_0 = randn(2,1);
V_0 = eye(2);

% Parameters
isTimeVariant = true;
coeff_alpha = 1*0;      % has no model uncertainty

% Simulate data
real_epsilon = 0.3;
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


%% KF: Standard Kalman filter
x_0 = [0;0];
tic
xhat_KF = KF(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0);
time_KF = toc/T;
%     disp(['Avg Time of KF: ' num2str(time_KF)]);
err_KF = sum((x - xhat_KF).^2,1);
RMSE_KF = sqrt(mean(err_KF));
%     disp(['---------------------RMSE of KF: ' num2str(RMSE_KF)]);

%% Breakdown test
LoadData;

[epsilon_len_contamination, ~] = size(Epsilon_Contamination_Parameter);

% RMSE v.s. rparameter epsilon
RMSE_Huber_Epsilon_Contamination_Array = zeros(epsilon_len_contamination,1);
for epsilon_index = 1:epsilon_len_contamination
	%% Note that when theta_2 := 0 and alpha := 0 (i.e., no parameter uncertainty), MKF reduces to HKF
    %%Huber Kalman filter (epsilon Contamination)
    x_0 = [0;0];
    epsilon = Epsilon_Contamination_Parameter(epsilon_index,1);
    k = Epsilon_Contamination_Parameter(epsilon_index,2);
    tic
    xhat_Huber_Epsilon_Contamination = Huber_Epsilon_Contamination(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0,k,epsilon);
    time_Huber_Epsilon_Contamination = toc/T;
%     disp(['Avg Time of Huber (Contamination): ' num2str(time_Huber_Epsilon_Contamination)]);
    err_Huber_Epsilon_Contamination = sum((x - xhat_Huber_Epsilon_Contamination).^2,1);
    RMSE_Huber_Epsilon_Contamination = sqrt(mean(err_Huber_Epsilon_Contamination));
%     disp(['---------------------RMSE of Huber (Contamination): ' num2str(RMSE_Huber_Epsilon_Contamination)]);
    RMSE_Huber_Epsilon_Contamination_Array(epsilon_index) = RMSE_Huber_Epsilon_Contamination;
end

[epsilon_len_normal, ~] = size(Epsilon_Normal_Parameter);

% RMSE v.s. rparameter epsilon
RMSE_Huber_Epsilon_Normal_Array = zeros(epsilon_len_normal,1);
for epsilon_index = 1:epsilon_len_normal
	%% Note that when theta_2 := 0 and alpha := 0 (i.e., no parameter uncertainty), MKF reduces to HKF
    % Huber Kalman filter (epsilon Normal)
    x_0 = [0;0];
    a = Epsilon_Normal_Parameter(epsilon_index,2);
    b = Epsilon_Normal_Parameter(epsilon_index,4);
    c = Epsilon_Normal_Parameter(epsilon_index,3);
    tic
    xhat_Huber_Epsilon_Normal = Huber_Epsilon_Normal(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0,a,b,c);
    time_Huber_Epsilon_Normal = toc/T;
%     disp(['Avg Time of Huber (Normal): ' num2str(time_Huber_Normal)]);
    err_Huber_Epsilon_Normal = sum((x - xhat_Huber_Epsilon_Normal).^2,1);
    RMSE_Huber_Epsilon_Normal = sqrt(mean(err_Huber_Epsilon_Normal));
%     disp(['---------------------RMSE of Huber (Normal): ' num2str(RMSE_Huber_Normal)]);
    RMSE_Huber_Epsilon_Normal_Array(epsilon_index) = RMSE_Huber_Epsilon_Normal;
end

%% Plot results
yline(RMSE_TMKF, 'k', 'linewidth', 2); hold on;
yline(RMSE_KF, 'r', 'linewidth', 2); hold on;
plot(Epsilon_Contamination_Parameter(:,1), RMSE_Huber_Epsilon_Contamination_Array, 'b',...
     Epsilon_Normal_Parameter(:,1), RMSE_Huber_Epsilon_Normal_Array, 'm',...
     'linewidth', 2); hold on;
xline(real_epsilon, 'g--', 'linewidth', 2);
leg1 = legend({'TMKF','KF','$\epsilon$-Contamination','$\epsilon$-Normal','$\epsilon_{real}$'}, 'Interpreter', 'latex', 'Location', 'northwest');

set(gca, 'FontSize', 18);
ylabel('RMSE','FontSize', 20, 'Interpreter', 'latex');
xlabel('$\epsilon$','FontSize', 20, 'Interpreter', 'latex');

set(leg1,'FontSize',15);

axis([0 0.5 RMSE_TMKF-0.2 RMSE_KF+0.2]);

box on;
