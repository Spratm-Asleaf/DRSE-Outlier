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

function [x,y,real_y,y0,true_F] = GenerateData(sys, x_0, T, coeff, is_TV)
    [m,n] = size(sys.H);
    
    dof = 3;

    y0 = sys.H * x_0 + sys.D * chol(sys.R) * randn(m,1);
    
    x = zeros(n,T);
    y = zeros(m,T);
    real_y = zeros(m,T);
    true_F = cell(1,T);
    Delta = 2 * rand - 1;   % for changable
                            % Delta = 1, for fixed 
    A_purt = sys.F + [0, coeff * Delta; 0, 0];
    
    x_prev = x_0;
    for t = 1 : T
        %A_purt
        x(:,t) = A_purt * x_prev + sys.G * chol(sys.Q) * randn(n,1);
        y(:,t) = sys.H * x(:,t) + sys.D * chol(sys.R) * trandn(m,1,dof);       % t-distribution measurement noise
        
        x_prev = x(:,t);
        
        true_F{t} = A_purt;
        
        if is_TV
            Delta = 2 * rand - 1;
            A_purt = sys.F + [0, coeff * Delta; 0, 0];
        end

        real_y(:,t) = y(:,t);
        
        %Add outliers
        %None
    end
end

% student-t distributed noise
function samples = trandn(m,n,dof)
    variance = dof/(dof-2);
    coef = 1/sqrt(variance);
    
    samples = coef * trnd(dof,m,n);     % "samples" has zero mean and unit variance
end