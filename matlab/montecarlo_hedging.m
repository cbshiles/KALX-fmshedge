function [S CALL_PRICE UPPER_ORDER LOWER_ORDER UPPER_REHEDGE LOWER_REHEDGE QUANTITY_S QUANTITY_BOND MTM_S MTM_BOND PNL E_FINAL_PNL STD_FINAL_PNL CONDITIONAL_E_FINAL_LOSS] = montecarlo_hedging(S0, K, dt, T, R, SIGMA, npaths, dSmin, flag_plot)


%% Variable Definitions

% INPUTS
% - S0 - Scalar representing the Initial price of the underlying instrument S
% - K - Scalar representing the strike of the option one will sell (write)
% - dt - Scalar representing the (finite) time step when one will place the new upper and lower limit orders to re-hedge. It should be informed in
%   years (e.g: dt = 1/(12*21*24) for 1 hour, since one year = 12 months *21 trading days * 24 hours)
% - T - Scalar representing the Total maturity of the option. Should also be informed in years. E.g: T = 1/12 = 1 month
% - SIGMA - Scalar representing the volatility of the underlying (e.g: SIGMA = 0.3 for 30% of annual vol)
% - npaths - Scalar represending the number of paths to be considered for the Monte Carlo simulation  
% - dSmin -  Scalar. It defines the minimum variation (on each dt) that the
% stock should move (upwards and downwards) in order to be re-hedged by the limit orders
% - flag_plot - Flag variable that indicates if one wants to plot the final results (if equal to 1) or not (otherwise)

% OUTPUTS (All Matrices show the time evolution in rows and the different paths organized in columns)
% - S - Matrix (time in rows, paths in columns) with the evolution of the underlying S 
% - CALL_PRICE - Matrix (time in rows, paths in columns) with the evolution of the price of the call option  
% - UPPER_ORDER - Matrix (time in rows, paths in columns) showing what should be the upper order (on that time and for that path) to be placed
% - LOWER_ORDER - Matrix (time in rows, paths in columns) showing what should be the lower order (on that time and for that path) to be placed
% - UPPER_REHEDGE - Matrix (time in rows, paths in columns) with the indicator variables related to the upper re-hedge. If it is equal to 1 it means
%   that the re-hedge was performed for that time step and for that path. Zero otherwise. 
% - LOWER_REHEDGE - Matrix (time in rows, paths in columns) with the indicator variables related to the lower re-hedge. If it is equal to 1 it means
%   that the re-hedge was performed for that time step and for that path. Zero otherwise.
% - QUANTITY_S - Matrix (time in rows, paths in columns) with the actual quantities (positions) on the risky asset for each time and for each path
% - QUANTITY_BOND - Matrix (time in rows, paths in columns) with the actual quantities (positions) on the money account for each time and for each path
% - MTM_S - Matrix (time in rows, paths in columns) with the MtM exposure (or simply P&L) for the position on S. It is the component-wise 
%   product of QUANTITY_S and S
% - MTM_BOND - Matrix (time in rows, paths in columns) with the MtM exposure (or simply P&L) for the position on the money account. It is the 
%   component-wise product of QUANTITY_BOND and exp(R*Time)
% - PNL - Matrix (time in rows, paths in columns)
% - E_FINAL_PNL - Scalar. It represents the expected value (mean across diff paths) of the terminal P&L of the total portfolio

%% Scale Conversions

Conv_Year2Hour = 24*12*21;
Conv_Year2Minute = 24*12*21*60;


%% Initial Definitions

nsteps = T/dt; % defines the number of time steps
SIGMA_step = SIGMA.*(sqrt(dt)); % defines the volatility related to a single time step dt

Z = normrnd(0,1,[nsteps npaths]); % Simulation of normally distributed draws to be used in the Monte Carlo simulation for S

S = S0.*(ones(nsteps+1, npaths)); % Initialize the matrix S with the value (constant) of S0
dS = zeros(nsteps, npaths); % Initialize the matrix with the differences dS with zeros

C0 = blsprice(S0,K,R,T,SIGMA); % Call price at time zero
CALL_PRICE = C0.*(ones(nsteps+1, npaths)); % Initializes the matrix of call prices with the same value of C0
DELTA0 = blsdelta(S0,K,R,T,SIGMA); % Initial delta hedge
ETA0 = (C0 - DELTA0*S0);

QUANTITY_S = DELTA0.*ones(nsteps+1, npaths); % Initializes the matrix of quantities on S with the starting value DELTA0
ETA = ETA0.*ones(nsteps+1, npaths); % Initializes the matrix of positions in the money account with the starting value ETA0


DELTA_UPPER = blsdelta(S0+dSmin,K,R,T,SIGMA).*ones(nsteps+1, npaths); % Initializes the matrix of upper deltas with the first value (first time step dt)
DELTA_LOWER = blsdelta(S0-dSmin,K,R,T,SIGMA).*ones(nsteps+1, npaths); % Initializes the matrix of lower deltas with the first value (first time step dt)

Time_vector = 0:dt:T; % Vector with the time evolution
Maturity_vector = T - Time_vector; % Vector with the maturity (T - time) evolution
Bond_vector = exp(R.*Time_vector); % Vector with the money market price exp(R*Time)

Maturity_vector(end) = 10^-30; % corrects the last maturity value (in order to use function blsdelta)

%% Monte Carlo Simulation for the prices of S, the Call option and the corresponding re-hedging strategies

for i=1:1:nsteps 

    dS(i,:) = (R*dt).*(S(i,:)) + (Z(i,:).*(S(i,:))).*(SIGMA_step); % simulates the (finite) variations of S on time (for different paths) 
    S(i+1,:) = S(i,:) + dS(i,:); % Computes the value of S on the next time step (for different paths)
    
    CALL_PRICE(i+1,:) = blsprice(S(i+1,:),K,R,Maturity_vector(i+1),SIGMA); % Computes the corresponding market value of the call option at time i+1 (for different paths)
    
    if (i<nsteps) % Computes what should be the next delta hedge in both the upper and lower conditions for S on the next time step
        DELTA_UPPER(i+1,:) = blsdelta(S(i+1,:)+dSmin,K,R,Maturity_vector(i+2),SIGMA);
        DELTA_LOWER(i+1,:) = blsdelta(S(i+1,:)-dSmin,K,R,Maturity_vector(i+2),SIGMA);
    end
    
end

UPPER_REHEDGE = [zeros(1,npaths); (dS>=dSmin)]; % Defines the dummy variable for the upper re-hedge as a function of dS
LOWER_REHEDGE = [zeros(1,npaths); (dS<=-dSmin)]; % Defines the dummy variable for the lower re-hedge as a function of dS

%% Computes the new quantity of stocks for each simulation and the corresponding quantity of money mkt account (self financing condition)
for i=1:1:nsteps
   
 QUANTITY_S(i+1,:) = QUANTITY_S(i,:).*(ones(1,npaths) - UPPER_REHEDGE(i+1,:) - LOWER_REHEDGE(i+1,:)) + DELTA_UPPER(i,:).*(UPPER_REHEDGE(i+1,:)) + DELTA_LOWER(i,:).*(LOWER_REHEDGE(i+1,:));    
 ETA(i+1,:) = (S(i+1,:).*(QUANTITY_S(i,:)-QUANTITY_S(i+1,:))+Bond_vector(i+1).*(ETA(i,:)))./(Bond_vector(i+1));
 
end

%% Defining the Total MtM (or P&L) for the Stock, the Money account, the Option and the total portfolio

QUANTITY_BOND = ETA; % Just another name

MTM_BOND = QUANTITY_BOND.*(repmat(Bond_vector',[1 npaths])); % Exposure (or MtM) on S
MTM_S = QUANTITY_S.*S; % Exposure (or MtM) on B (money account)

PNL = -CALL_PRICE + MTM_S + MTM_BOND; % Portfolio total value (or P&L evolution)
E_FINAL_PNL = mean(PNL(end,:)); % Expected value of the terminal P&L

STD_FINAL_PNL = std(PNL(end,:)); % Standard deviation of the terminal P&L
CONDITIONAL_E_FINAL_LOSS = -sum((PNL(end,:)<0).*PNL(end,:))/(length((PNL(end,:)<0))); % Conditional (expected) terminal loss (negative P&L)

UPPER_ORDER = DELTA_UPPER - QUANTITY_S; % net value of the order placed for the next time step if the stock reachs the upper condition
LOWER_ORDER = DELTA_LOWER - QUANTITY_S; % net value of the order placed for the next time step if the stock reachs the upper condition


%% Final Corrections
% Corrects the matrices to take away the last row (because at the final time (maturity) there is no need for defining a hedge and placing orders
DELTA_UPPER = DELTA_UPPER((1:(end-1)),:);
DELTA_LOWER = DELTA_LOWER((1:(end-1)),:);
UPPER_ORDER = UPPER_ORDER((1:(end-1)),:);
LOWER_ORDER = LOWER_ORDER((1:(end-1)),:);

nplot1 = 10;
nplot2 = 50;

%% Plotting final results
if (flag_plot==1)
    % The first figure plots the indicator variables related to the upper re-hedge and the lower re-hedge
    figure(1)
    subplot(2,1,1)
    plot(Time_vector,UPPER_REHEDGE(:,1:nplot1));
    axis([Time_vector(1) Time_vector(end) -0.5 1.5]);
    title('Upper re-hedge (1-yes, 0-no) - First 10 paths')
    xlabel('Time (years)');
    subplot(2,1,2)
    plot(Time_vector,LOWER_REHEDGE(:,1:nplot1));
    axis([Time_vector(1) Time_vector(end) -0.5 1.5]);
    title('Lower re-hedge (1-yes, 0-no) - First 10 paths')
    xlabel('Time (years)');
    % The second figure plots for each path the value of the upper and the lower limit orders that were placed on each time
    figure(2)
    subplot(2,1,1)
    plot(Time_vector(1:(end-1)),UPPER_ORDER(:,1:nplot2));
    title('Upper order placed - First 50 paths')
    xlabel('Time (years)'); ylabel('Quantity of S');
    subplot(2,1,2)
    plot(Time_vector(1:(end-1)),LOWER_ORDER(:,1:nplot2));
    title('Lower order placed - First 50 paths')
    xlabel('Time (years)'); ylabel('Quantity of S');
    % This third figure plots the time evolution of the Underlying and the money market account as their corresponding positions (or number of shares)
    figure(3)
    subplot(2,2,1)
    plot(Time_vector,S(:,1:nplot2))
    title('Underlying Evolution - First 50 paths')
    xlabel('Time (years)'); ylabel('$');
    axis([Time_vector(1) Time_vector(end) 0.9*min(min(S)) 1.1*max(max(S))]);
    subplot(2,2,2)
    plot(Time_vector,Bond_vector)
    title('Money Account Evolution - First 50 paths')
    xlabel('Time (years)'); ylabel('$');
    axis([Time_vector(1) Time_vector(end) 0.9999*min(Bond_vector) 1.000001*max(Bond_vector)]);
    subplot(2,2,3)
    plot(Time_vector,QUANTITY_S(:,1:nplot2))
    title('Underlying Position - First 50 paths')
    xlabel('Time (years)'); ylabel('Units of S');
    axis([Time_vector(1) Time_vector(end) 0.9*min(min(QUANTITY_S)) 1.1*max(max(QUANTITY_S))]);
    subplot(2,2,4)
    plot(Time_vector,QUANTITY_BOND(:,1:nplot2))
    title('Money Market Position - First 50 paths')
    xlabel('Time (years)'); ylabel('Units of Money Mkt');
    axis([Time_vector(1) Time_vector(end) 0.9*min(min(QUANTITY_BOND)) 1.1*max(max(QUANTITY_BOND))]);
    % The fourth figure plots
    figure(4)
    subplot(2,2,1)
    plot(Time_vector,MTM_S(:,1:nplot2));
    title('Underlying P&L (MtM) - First 50 paths')
    xlabel('Time (years)'); ylabel('$');
    axis([Time_vector(1) Time_vector(end) 0.9*min(min(MTM_S)) 1.1*max(max(MTM_S))]);
    subplot(2,2,2)
    plot(Time_vector,MTM_BOND(:,1:nplot2));
    title('Money Account P&L (MtM) - First 50 paths')
    xlabel('Time (years)'); ylabel('$');
    axis([Time_vector(1) Time_vector(end) 0.9*min(min(MTM_BOND)) 1.1*max(max(MTM_BOND))]);
    subplot(2,2,3)
    plot(Time_vector,-CALL_PRICE(:,1:nplot2));
    title('Call Option P&L (MtM) - First 50 paths')
    xlabel('Time (years)'); ylabel('$');
    axis([Time_vector(1) Time_vector(end) 1.1*min(min(-CALL_PRICE)) 0.9*max(max(-CALL_PRICE))]);
    subplot(2,2,4)
    plot(Time_vector,PNL(:,1:nplot2));
    title('Total P&L (MtM) - First 50 paths')
    xlabel('Time (years)'); ylabel('$');
    axis([Time_vector(1) Time_vector(end) 1*min(min(PNL)) 1.1*max(max(PNL))]);
    
    figure(5)
    hist(PNL(end,:),100);
    title('Histogram of terminal P&L'), xlabel('P&L ($)');
    
end




end