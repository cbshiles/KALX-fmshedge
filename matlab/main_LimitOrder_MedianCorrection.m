
function [S_END C_END Q_S_END Q_B_END PNL_END N_UPPER_HITS N_LOWER_HITS] = main_LimitOrder_MedianCorrection(S0, K, R, T, SIGMA, dSmin)

% this function simulates a sigle path of the limit order hedging strategy

%% Inputs / Outputs

% Inputs
% - S0 - Initial value of the spot price 
% - K - Strike price
% - R - Risk free (constant) interest rate
% - T - Option큦 maturity (given in hours)
% - SIGMA - Black and Scholes flat volatility
% - dSmin - Symmetric distance between one limit order and the current value of the spot

% Outputs
% - S_END - terminal value of the spot price
% - C_END - terminal value of the call option
% - Q_S_END - terminal quantity of stocks held
% - Q_B_END - terminal quantity of units of the money account held
% - PNL_END - terminal value of the portfolio큦 P&L (portfolio is short in a call and long in the replicating portfolio)
% - N_UPPER_HITS - number of times that we hit the upper barrier
% - N_LOWER_HITS - number of times that we hit the lower barrier

%% Initial definitions

year2hour = (12*21*24); % Time conversion factor (years - hours)

TtoExp = T; % initial value of the time to expiration
S_CURRENT = S0; % initial value of the spot
C_CURRENT = blsprice(S_CURRENT, K, R, (TtoExp)/year2hour, SIGMA); % initial value of the call option
% initial levels of the upper / lower limit orders
S_UPPER_CURRENT = S_CURRENT +dSmin;
S_LOWER_CURRENT = S_CURRENT -dSmin;
% initial quantity of stocks held. At time zero the problem is perfectly delta hedged and the quantity is equal to the B&S delta
Q_S_CURRENT = blsdelta(S_CURRENT, K, R, (TtoExp)/year2hour, SIGMA);

B_CURRENT = 1; % initial value of the money account
Q_B_CURRENT = (C_CURRENT - S_CURRENT.*Q_S_CURRENT)/B_CURRENT; % initial units held in the money account

PNL_CURRENT = 0; % starting value of the hedged portfolio (zero by construction)
% initializes the number of upper and lower hits
N_UPPER_HITS = 0; 
N_LOWER_HITS = 0;

%a=0;b=0;

%% looping for the hitting times and updated positions before the option큦 expiration

while (TtoExp>0) % condition that the option did not expire yet
    % computes the median of the hitting times for both upper and lower levels
    [MEDIAN_HIT_TIME_UPPER MEDIAN_HIT_TIME_LOWER] = median_hitting_time(S_CURRENT, S_UPPER_CURRENT, S_LOWER_CURRENT, R, SIGMA);
    
    % computes the time to maturity to be used in the delta expression (for both upper and lower barriers)
    % the 0.01 correction is just to ensure this time wont be equal to zero when the option has expired
    T_DELTA_UPPER = (TtoExp-MEDIAN_HIT_TIME_UPPER)*(TtoExp>MEDIAN_HIT_TIME_UPPER) + 0.01*(TtoExp<MEDIAN_HIT_TIME_UPPER);
    T_DELTA_LOWER = (TtoExp-MEDIAN_HIT_TIME_LOWER)*(TtoExp>MEDIAN_HIT_TIME_LOWER) + 0.01*(TtoExp<MEDIAN_HIT_TIME_LOWER);
    
    % computes the expected delta at the next hit (still unknown when it will hit) based on the barrier levels and the medians of the hitting times
    DELTA_UPPER_CURRENT = blsdelta(S_UPPER_CURRENT, K, R, (T_DELTA_UPPER)/year2hour, SIGMA);
    DELTA_LOWER_CURRENT = blsdelta(S_LOWER_CURRENT, K, R, (T_DELTA_LOWER)/year2hour, SIGMA);
    
    % computes the sizes of the limit orders based on the final expected deltas in both conditions and the current number of stocks held
    ORDER_UPPER_CURRENT = DELTA_UPPER_CURRENT - Q_S_CURRENT;
    ORDER_LOWER_CURRENT = DELTA_LOWER_CURRENT - Q_S_CURRENT;
    
    % calls the function that generates the next hitting time and what order was hit/executed (upper or lower)
    [HIT_TIME UPPERLOWER_FLAG] = hitting_times_bisection(S_CURRENT, S_UPPER_CURRENT, S_LOWER_CURRENT, R, SIGMA);
    
    %a = a + (HIT_TIME >= 0.5*(MEDIAN_HIT_TIME_UPPER+MEDIAN_HIT_TIME_LOWER))
    %b = b + (HIT_TIME < 0.5*(MEDIAN_HIT_TIME_UPPER+MEDIAN_HIT_TIME_LOWER))
    
    N_UPPER_HITS = N_UPPER_HITS + (UPPERLOWER_FLAG==1);  % updates the number of upper hits
    N_LOWER_HITS = N_LOWER_HITS + (UPPERLOWER_FLAG==-1); % updates the number of lower hits
    
    TtoExp = TtoExp - HIT_TIME; % Updates the time to expiration
    
    OptionNotExpired = (TtoExp>0); % flag variable that will be equal to 1 if the option did not expire yet and zero otherwise
    % this flag is used to determine when we apply the updates in the following variables or not
    
    % updates the value of the money market account
    B_CURRENT = (exp(R*(T-TtoExp)/year2hour))*(OptionNotExpired) + (exp(R*(T)/year2hour))*(1-OptionNotExpired);
    
    % updates the value of the spot based on what barrier was hit (if it did not expire)
    S_CURRENT = (0.5*(UPPERLOWER_FLAG + 1)*S_UPPER_CURRENT + 0.5*(1 - UPPERLOWER_FLAG)*S_LOWER_CURRENT)*(OptionNotExpired) + S_CURRENT*(1-OptionNotExpired);
    
    % updates the value of the call option based on the new spot and time (if it did not expire)
    C_CURRENT = blsprice(S_CURRENT, K, R, (TtoExp)/year2hour*(OptionNotExpired), SIGMA)*(OptionNotExpired);
    
    % assigns the previous positions (units held)
    Q_S_PAST = Q_S_CURRENT;
    Q_B_PAST = Q_B_CURRENT;
    
    % computes the new/updated quantity of stocks held - based on what order was executed (if it did not expire yet)
    Q_S_CURRENT = (Q_S_PAST + (0.5*(UPPERLOWER_FLAG + 1)*ORDER_UPPER_CURRENT + 0.5*(1 - UPPERLOWER_FLAG)*ORDER_LOWER_CURRENT))*(OptionNotExpired) + Q_S_PAST*(1-OptionNotExpired);
    
    % computes the new/updated quantity of bond account held based on the self financing condition
    Q_B_CURRENT = (Q_B_PAST*B_CURRENT + (Q_S_PAST - Q_S_CURRENT)*S_CURRENT)/B_CURRENT;
    
    % computes the new P&L of the portfolio (short call and long replicating portfolio)
    PNL_CURRENT = ((Q_S_CURRENT*S_CURRENT + Q_B_CURRENT*B_CURRENT) - C_CURRENT)*(OptionNotExpired);
    
    % computes the new levels of the limit orders to be placed (both upper and lower)
    S_UPPER_CURRENT = (S_CURRENT +dSmin)*(OptionNotExpired) + S_UPPER_CURRENT*(1-OptionNotExpired);
    S_LOWER_CURRENT = (S_CURRENT -dSmin)*(OptionNotExpired) + S_LOWER_CURRENT*(1-OptionNotExpired);
    
    
    
end

% if the hitting time generated close to the option expiration leads to a time that is after T, it means that the option 
% expired before hitting the current barrier levels

% in this case the terminal value of the spot will surely fall between the current barrier levels. We used a uniform approx 
% for the terminal value of S, since the barriers are supposed to be close

%% Generates the market values for the condition of the option큦 expiration

if (OptionNotExpired==0) % condition of option being expired
    
    S_CURRENT = S_LOWER_CURRENT + rand(1,1)*(S_UPPER_CURRENT - S_LOWER_CURRENT); % terminal value of the spot price 
    
    C_CURRENT = (S_CURRENT - K)*(S_CURRENT>K); % terminal value of the call (option큦 payoff)
    
    PNL_CURRENT = ((Q_S_CURRENT*S_CURRENT + Q_B_CURRENT*B_CURRENT) - C_CURRENT); % terminal value of the P&L
    
end

%% Final Outputs

% definines the final outputs as the terminal conditions of the hedging
% strategy for the single path we simulated
S_END = S_CURRENT;
C_END = C_CURRENT;
Q_S_END = Q_S_CURRENT;
Q_B_END = Q_B_CURRENT;
PNL_END = PNL_CURRENT;



end