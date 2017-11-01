% Stop Loss Start Gain Simulator

% 25  days till the expiration date.
% Daily Hedging at the end of the day.
tic
cost = zeros(1,100);
vars = zeros(1,100);

for i = 1:100       
[S buysell Cost meancost var] = simple_hedging(99, 100, 1/(365*10), 25/365, 0, .5, 100000);
cost(i) = meancost;
vars(i) = var;
end
toc
% mean(cost) = -0.0161;
% mean(vars) = 81.6618;


% check the price at 10 times a day. then decide whether to buy the stock
% or not

% for i = 1:100       
% [S buysell Cost meancost var] = simple_hedging(99, 100, 1/(365*10), 25/365, 0, .5, 10000);
% cost(i) = meancost;
% vars(i) = var;
% end

% mean(cost) = 0.0085
% mean(vars) = 83.4

