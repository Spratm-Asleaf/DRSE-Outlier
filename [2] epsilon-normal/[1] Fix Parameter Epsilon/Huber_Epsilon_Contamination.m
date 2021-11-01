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

function [hat_X,P] = Huber_Epsilon_Contamination(Y,F,G,H,Q,R,Pi_0,hat_x_0,k,epsilon)
[~,N] = size(Y);
P = Pi_0;
X = hat_x_0;
n = length(hat_x_0);
hat_X = [];

for i=1:N
    
    X = F*X;
    Z_ = H*X;

    P = F*P*F' + G*Q*G';
    
    S = H*P*H' + R;
    S_square = chol(S^(-1));
    
    X = X + P*H'*S_square*phi(S_square*(Y(i) - Z_), k);
    
    hat_X = [hat_X X];

    P = P - P*H'*S^(-1)*H*P*(1-epsilon)*(1-2*Phi(-k));
end

end

function value = phi(normalized_innovation, k)
    m = length(normalized_innovation);
    value = normalized_innovation;
    
    for j = 1:m
        if normalized_innovation(j) < -k
            %warning('here 1')
            value(j) = -k;
        elseif normalized_innovation(j) > k
            value(j) = k;
            %warning('here 2')
        end
    end
end

function val = Phi(x)
    val = (1/2)*(1+erf(x/sqrt(2)));
end