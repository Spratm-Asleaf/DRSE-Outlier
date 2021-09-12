expected_gamma = 0:0.001:0.1;
gamma_len = length(expected_gamma);
RMSE_array = zeros(gamma_len,1);
for gamma_index = 1:gamma_len
gamma = expected_gamma(gamma_index);
gamma_1 = 1 - gamma;
gamma_2 = 1 + gamma;
epsilon = 0.05;
k = 1.4;      % use 1.4 or 2, make no big difference
xhat_MKF = MKF(y,sys.F,sys.G,sys.H,sys.Q,sys.R,V_0,x_0,k,epsilon,gamma_2);
err_MKF = sum((x-xhat_MKF).^2,1);
RMSE_MKF = sqrt(mean(err_MKF));
RMSE_array(gamma_index) = RMSE_MKF;
end
plot(expected_gamma,RMSE_array)