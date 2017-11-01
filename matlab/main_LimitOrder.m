%% Main flie to run a full path of the hedging strategy

%% SET PARAMETERS HERE
S_CURRENT = 100;
S_UPPER = 101;
S_LOWER = 99;
dSmin = 1;

R = 0;
SIGMA = 0.7;
K = 100;
T = 1/12; %Initial Time to Expiration in Years

%%
EXP_HIT_TIME = 0;%Felipe's function goes here
numSHARE = blsdelta(S_CURRENT, K, R, T, SIGMA); %Initially You will hold this many shares
CASH = -numSHARE;
TtoExp = T; % Initially Time to Expiration is set to T
year2hour = (12*21*24);
DELTA_UPPER = blsdelta(S_UPPER, K, R, (T-EXP_HIT_TIME)/year2hour, SIGMA);
DELTA_LOWER = blsdelta(S_LOWER, K, R, (T-EXP_HIT_TIME)/year2hour, SIGMA);

[HIT_TIME UPPERLOWER] = hitting_times_bisection(S_CURRENT, S_UPPER, S_LOWER, R, SIGMA);%generate the next hitting time and the upperlower flag
while UPPERLOWER ~= 0
    newSHARE = 1/2.*((UPPERLOWER+1).*DELTA_UPPER + (1-UPPERLOWER).*DELTA_LOWER);  
    BOT = newSHARE - numSHARE;
    CASH = CASH - BOT.*S_CURRENT;
    TtoExp = TtoExp - HIT_TIME; %Time to Expiration has now decreased by the HIT_TIME amount
    S_UPPER_LAST = S_UPPER;
    S_LOWER_LAST = S_LOWER;
    [S_CURRENT S_UPPER S_LOWER, DELTA_UPPER, DELTA_LOWER] = LimitOrder(UPPERLOWER, S_UPPER_LAST, S_LOWER_LAST, K, TtoExp, R, SIGMA, dSmin);
    numSHARE = newSHARE;
    [HIT_TIME UPPERLOWER] = hitting_times_bisection(S_CURRENT, S_UPPER, S_LOWER, R, SIGMA);%generate the next hitting time and the upperlower flag
end

%Once we go past the Expiration, we randomly generate the PRICE AT EXP and
%then compute the TERMINAL PNL

S_EXP = rand * (S_UPPER - S_LOWER) + S_LOWER;
PNL = CASH + numSHARE.*S_EXP;