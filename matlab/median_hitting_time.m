function [MEDIAN_HIT_TIME_UPPER MEDIAN_HIT_TIME_LOWER] = median_hitting_time(S_CURRENT, S_UPPER, S_LOWER, R, SIGMA)

%% Function to compute random hitting times for a GEOMETRIC BM with two barrrier levels (upper and lower)
%% Inputs and Outputs

% Inputs
% S_CURRENT - Current value of the spot price (scalar)
% S_UPPER - Current level of the upper barrier (scalar)
% S_LOWER - Current level of the lower barrier (scalar)
% R - Risk free interest rate (drift for the risk neutral probability measure)
% SIGMA - Volatility (difusion term)

% outputs
% MEDIAN_HIT_TIME_UPPER - median of the hitting time at the upper barrier
% MEDIAN_HIT_TIME_LOWER - median of the hitting time at the lower barrier

%% Algorithm

year2hour = (12*21*24); % rescale the time dimension to hours

R = R/year2hour; % rescale R to hours
SIGMA = SIGMA/sqrt(year2hour); % rescale SIGMA to hours

alpha = (1/SIGMA)*log(S_UPPER/S_CURRENT); % defines the equivalent level of the upper barrier of the GBM to the standard BM
beta = (1/SIGMA)*log(S_LOWER/S_CURRENT); % defines the equivalent level of the upper barrier of the GBM to the standard BM
nu = (R/SIGMA) - (SIGMA/2); % defines parameter nu
U = [0.5 0.5]; % generates 0.5 and 0.5 as the uniform "random numbers" to ensure we can compute the median of both 
% upper hitting times and lower hitting times

MEDIAN_HIT_TIME_UPPER = fzero(@(t) equation1(t, alpha, nu, U(1)), [0 1000000000]);
MEDIAN_HIT_TIME_LOWER = fzero(@(t) equation2(t, beta, nu, U(2)), [0 1000000000]);


end

%% function to use in the bisection method for the upper barrier
function eq1 = equation1(t, alpha, nu, u)

eq1 = (u - normcdf((-alpha + nu*t)/sqrt(t)) - exp(2*nu*alpha)*normcdf((-nu*t - alpha)/sqrt(t)));


end

%% function to use in the bisection method for the lower barrier
function eq2 = equation2(t, beta, nu, u)

eq2 = (u - normcdf((beta - nu*t)/sqrt(t)) - exp(2*nu*beta)*normcdf((nu*t + beta)/sqrt(t)));

end
