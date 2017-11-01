%Simulation method for option Pricing
function [pnl,standard_error]= delta_hedge(Rebalancing_time)
% Define the parameters for the option
K = 105;
T = 22/256;
n = 10000;
N = Rebalancing_time;
S0 = 100;
sigma = 0.1;
r = 0.05;
c = blsprice(S0,K,r,T,sigma,0.0);
% Define what are needed for calculations
Z = randn(n,N);
S = zeros(n,N);
t_increment = T/Rebalancing_time;
timeToMaturity = T-t_increment;
S(:,1)=S0*exp((r-0.5*sigma*sigma)*T/N+sigma*Z(:,1)*sqrt(T/N));
Delta_prev = calculate_callDelta(S(:,1),K,r,sigma,timeToMaturity);
Stock_Position_P = Delta_prev.*S(:,1);
Bond_Position = c-Stock_Position_P;
for i = 2:N
Bond_Position = Bond_Position*exp(r*t_increment); %Bond has grown to this value
S(:,i)=S(:,(i-1)).*exp((r-0.5*sigma*sigma)*T/N+sigma*Z(:,i)*sqrt(T/N)); % New stock process
timeToMaturity = timeToMaturity-t_increment; % time has passed
if i Delta_new = calculate_callDelta(S(:,i),K,r,sigma,timeToMaturity);
Stock_adjustment = (Delta_new-Delta_prev).*S(:,i);
Bond_Position = Bond_Position - Stock_adjustment;
Stock_Position_P = Delta_new.*S(:,i);
else
Stock_Position_P=Delta_new.*S(:,N);
end
Delta_prev = Delta_new;
end
Final_value = Stock_Position_P+Bond_Position;
call = max(S(:,N)-K,0);
pnl = Final_value-call;
hist(pnl,500);
title(strcat('Difference for n =',num2str(Rebalancing_time)));
standard_error=std(pnl)/sqrt(n);
pnl = mean(pnl);
end
function delta = calculate_callDelta(S0,K,R,Sigma,T)
delta = (log(S0/K)+(R+0.5*Sigma^2)*T)/(Sigma*sqrt(T));
delta = normcdf(delta);
end