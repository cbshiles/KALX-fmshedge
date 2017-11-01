function [S buysell Cost_Hedging mean_cost var_cost] = simple_hedging(S0, K, dt, T, R, SIGMA, npaths)

% This code tries to simulates the strategy in which we buy the call option
% when the underlying price goes over the stirke price K and sells the
% call option when the price goes below.
% This simulation does NOT take the transaction cost into the account.

%% Initial Definitions

nsteps = T/dt; % defines the number of time steps
SIGMA_step = SIGMA.*(sqrt(dt)); % defines the volatility related to a single time step dt

Z = normrnd(0,1,[nsteps npaths]); % Simulation of normally distributed draws to be used in the Monte Carlo simulation for S

S = S0.*(ones(nsteps+1, npaths)); % Initialize the matrix S with the value (constant) of S0
dS = zeros(nsteps, npaths); % Initialize the matrix with the differences dS with zeros
buysell = zeros(nsteps+1, npaths);

C0 = blsprice(S0,K,R,T,SIGMA); % Call price at time zero
CALL_PRICE = C0.*(ones(nsteps+1, npaths)); % Initializes the matrix of call prices with the same value of C0

Time_vector = 0:dt:T; % Vector with the time evolution
Maturity_vector = T - Time_vector; % Vector with the maturity (T - time) evolution

for i=1:1:nsteps 

    dS(i,:) = (R*dt).*(S(i,:)) + (Z(i,:).*(S(i,:))).*(SIGMA_step); % simulates the (finite) variations of S on time (for different paths) 
    S(i+1,:) = S(i,:) + dS(i,:); % Computes the value of S on the next time step (for different paths)
end


for j = 1:1:npaths
    if S(1,j) > K
        buysell(1,j) = 1;
    end
        
    for i = 1:1:nsteps
        if (S(i,j)<K && S(i+1,j) >= K)
            buysell(i+1,j) = 1;
        elseif (S(i,j)>= K && S(i+1,j) < K)
            buysell(i+1,j) = -1;
        end
    end
end

PNL = S .* buysell;
final_pos = sum(buysell, 1);
Cost_Hedging = sum(PNL, 1)-final_pos.*S(end,:);
mean_cost = mean(Cost_Hedging);
var_cost = var(Cost_Hedging);
hist(Cost_Hedging,100)
title('Stop Loss Start Gain Cost')
xlabel('Cost of Hedging')
ylabel('Frequency')

