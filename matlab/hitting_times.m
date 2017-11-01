function [HIT_TIME UPPERLOWER_FLAG] = hitting_times(S_CURRENT, S_UPPER, S_LOWER, R, SIGMA)

%% Function to compute random hitting times for a GEOMETRIC BM with two barrrier levels (upper and lower)

%% Inputs and Outputs

% Inputs
% S_CURRENT - Current value of the spot price (scalar)
% S_UPPER - Current level of the upper barrier (scalar)
% S_LOWER - Current level of the lower barrier (scalar)
% R - Risk free interest rate (drift for the risk neutral probability measure)
% SIGMA - Volatility (difusion term)

% outputs
% HIT_TIME - hitting time for either upper or lower levels
% UPPERLOWER_FLAG - flag indicating if the hit happened in the upper (1) or lower (-1) levels

%% Algorithm

year2hour = (12*21*24); % rescale the time dimension to hours

R = R/year2hour; % rescale R to hours
SIGMA = SIGMA/sqrt(year2hour); % rescale SIGMA to hours

alpha = (1/SIGMA)*log(S_UPPER/S_CURRENT); % defines the equivalent level of the upper barrier of the GBM to the standard BM
beta = (1/SIGMA)*log(S_LOWER/S_CURRENT); % defines the equivalent level of the upper barrier of the GBM to the standard BM
nu = (R/SIGMA) - (SIGMA/2); % defines parameter nu
U = rand(1,2); % generates two samples of the uniform distribution - the first for the upper hit and the second for the lower hit

 % solves the equations by using optimization routines
t_upper = fminsearch(@(t) equation1(t, alpha, nu, U(1)), 0.5*exp(alpha*nu-abs(alpha)*sqrt(nu^2+1)));
t_lower = fminsearch(@(t) equation2(t, beta, nu, U(2)), 0.5*exp(beta*nu-abs(beta)*sqrt(nu^2+1)));

% defines the hitting time as the minimum between the stopping times of the upper and lower levels
HIT_TIME = min(t_upper,t_lower);
% defines the flag variable with 1 (if the first hit happened for the upper level) or -1 (if the first hit happened for the lower level)
UPPERLOWER_FLAG = 1*(t_upper>=t_lower) - 1*(t_upper<t_lower);


end

% function to use in a optimization routine to find the first passage time for the upper barrier
function eq1 = equation1(t, alpha, nu, u)

eq1 = (u - normcdf((-alpha + nu*t)/sqrt(t)) - exp(2*nu*alpha)*normcdf((-nu*t - alpha)/sqrt(t)))^2;


end

% function to use in a optimization routine to find the first passage time for the lower barrier
function eq2 = equation2(t, beta, nu, u)

eq2 = (u - normcdf((beta - nu*t)/sqrt(t)) - exp(2*nu*beta)*normcdf((nu*t + beta)/sqrt(t)))^2;

end