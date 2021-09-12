%{
Online supplementary materials of the paper titled 
"Distributionally Robust State Estimation for Linear Systems Subject to Uncertainty and Outlier"
Authored By: Shixiong Wang, and Zhisheng Ye
From the Department of Industrial Systems Engineering and Management, National University of Singapore

@Author: Shixiong Wang
@Date: 10 May 2021
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE-Outlier

Acknowledgements:

"WKF.m", "FrankWolfe.m" and "tauKF.m" are adapted from sources at
    "https://github.com/sorooshafiee/WKF".

For "tauKF.m", we also acknowledge M. Zorz.
    "https://mathworks.com/matlabcentral/fileexchange/54308-robust-kalman-filtering-package?focused=5751770&tab=function"
%}

function [hat_X,P] = TMKF(Y,real_F,G,H,Q,R,Pi_0,hat_x_0)
[~,N] = size(Y);
P = Pi_0;
X = hat_x_0;
n = length(hat_x_0);
hat_X = [];

for i=1:N
    F = real_F{i};
    X = F*X;
    Z_ = H*X;

    P = F*P*F' + G*Q*G';

    K = P*H'*(H*P*H' + R)^-1;
    X = X + K*(Y(:,i) - Z_);
    
    hat_X = [hat_X X];

    P = (eye(n,n) - K*H)*P*(eye(n,n) - K*H)' + K*R*K';
end
