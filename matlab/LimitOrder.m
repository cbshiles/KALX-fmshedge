function [S_CURRENT S_UPPER S_LOWER, DELTA_UPPER, DELTA_LOWER] = LimitOrder(UPPERLOWER, S_UPPER_LAST, S_LOWER_LAST, K, T, R, SIGMA, dSmin)

%% Function to generate the process from one stopping time to another.

%% Inputs and Outputs

% Inputs
% UPPERLOWER - V indicating whether the stock price has reached the upper or lower limit. 
%                 1 - upper, 0 - never hit, -1 - lower
% HIT_TIME - Time that the price reached the limit. Also another output from hitting_times function
% R - Risk free interest rate
% SIGMA - Volatility (diffusion term)
% dSmin - Price change before we rehedge.(Difference between the limit
%        order and the spot price)

% Outputs
% S_CURRENT - spot price
% S_UPPER - price of the new upper limit order
% S_LOWER - price of the new lower limit order
% DELTA_UPPER - size of the new upper limit order
% DELTA_LOWER - size of the new lowerlimit order


%% Algorithm

%This is an expressino will compute the current spot price.
% S_CURRENT = S_UPPER_LAST if UPPERLOWER = 1
% S_CURRENT = 0            if UPPERLOWER = 0
% S_CURRENT = S_LOWER_LAST if UPPERLOWER = -1
S_CURRENT = 1/2 .* abs(UPPERLOWER) * ((UPPERLOWER+1) .* S_UPPER_LAST + (1-UPPERLOWER) .* S_LOWER_LAST);

% Set next limit prices
S_UPPER = S_CURRENT + dSmin;
S_LOWER = S_CURRENT - dSmin;


% Next six lines to estimate the expected hitting time
% Probably wise to revise this part of the code to speed up the process
iter = 1000;
hit_times = zeros(1,1000);
for i = 1:iter
    hit_times(i) = hitting_times_bisection(S_CURRENT, S_UPPER, S_LOWER, R, SIGMA);
end
EXP_HIT_TIME = mode(hit_times);


%If the EXPECTED HITTING TIME is greater than the TIME TO EXPIRATION, we
%will esitmate the delta assuming we hit the limit orders half way.
if EXP_HIT_TIME > T
    EXP_HIT_TIME = T./2;
end

% Computes the size of the limit order
year2hour = (12*21*24);
DELTA_UPPER = blsdelta(S_UPPER, K, R, (T-EXP_HIT_TIME)/year2hour, SIGMA);
DELTA_LOWER = blsdelta(S_LOWER, K, R, (T-EXP_HIT_TIME)/year2hour, SIGMA);

