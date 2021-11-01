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


%% Tables below are from: 
% Huber (2009). Robust Statistics (2nd edition). Wiley.

%% Page 85, Exhibit 4.3
Epsilon_Contamination_Parameter = [
% Epsilon       K
0               inf
0.001           2.630 
0.002           2.435 
0.005           2.160 
0.01            1.945 
0.02            1.717 
0.05            1.399 
0.10            1.140 
0.15            0.980 
0.20            0.862 
0.25            0.766 
0.3             0.685 
0.4             0.550 
0.5             0.436
];

%% Page 90, Exhibit 4.6
Epsilon_Normal_Parameter = [
% Epsilon   a       c       b
0           0       1.4142  inf
0.001       0.6533  1.3658  2.4364
0.002       0.7534  1.3507  2.2317
0.005       0.9118  1.3234  1.9483
0.01        1.0564  1.2953  1.7241
0.02        1.2288  1.2587  1.4921
0.03033     1.3496  1.2316  1.3496
];



