function [E_FINAL_PNL_MATRIX STD_FINAL_PNL_MATRIX CONDITIONAL_E_FINAL_LOSS_MATRIX E_INTER_HEDGING_MATRIX EXEC_RATIO_MATRIX] = hedging_parameter_space_v2(S0, K, T, R, SIGMA, npaths)

tic;

%% Scale Conversions

Conv_Year2Hour = 24*12*21;
Conv_Year2Minute = 24*12*21*60;

%% Parameter Space Plotting

dt = 1/(24*12*21*10); % 6 min for the Brownian motion
dSmin_min = 0.05;

%[1 2 3 4 6 8 12]

t_factor_vector = [5 10 15 20 30 40 60];
dtau_vector = t_factor_vector.*(dt);
dtau_vector_hour = dtau_vector.*Conv_Year2Hour;

dSmin_vector = dSmin_min:dSmin_min:20*dSmin_min;

E_FINAL_PNL_MATRIX = zeros(length(t_factor_vector),length(dSmin_vector));
STD_FINAL_PNL_MATRIX = zeros(length(t_factor_vector),length(dSmin_vector));
CONDITIONAL_E_FINAL_LOSS_MATRIX = zeros(length(t_factor_vector),length(dSmin_vector));
E_INTER_HEDGING_MATRIX = zeros(length(t_factor_vector),length(dSmin_vector));
EXEC_RATIO_MATRIX = zeros(length(t_factor_vector),length(dSmin_vector));

    for (j=1:1:length(dSmin_vector)) 
        for (i=1:1:length(t_factor_vector))
                   
            [S CALL_PRICE UPPER_ORDER LOWER_ORDER UPPER_REHEDGE LOWER_REHEDGE QUANTITY_S QUANTITY_BOND MTM_S MTM_BOND PNL E_FINAL_PNL_MATRIX(i,j) STD_FINAL_PNL_MATRIX(i,j) CONDITIONAL_E_FINAL_LOSS_MATRIX(i,j) E_INTER_HEDGING_MATRIX(i,j) EXEC_RATIO_MATRIX(i,j)] = montecarlo_hedging_v2(S0, K, dt, t_factor_vector(i), T, R, SIGMA, npaths, dSmin_vector(j), 2);
        
        end
    
    
    end

E_INTER_HEDGING_MATRIX_hours = E_INTER_HEDGING_MATRIX.*Conv_Year2Hour;

STD_FINAL_PNL_MATRIX_COLUMN = [];
EXEC_RATIO_MATRIX_COLUMN = [];

for (j=1:1:length(dSmin_vector)) 
    
    STD_FINAL_PNL_MATRIX_COLUMN = [STD_FINAL_PNL_MATRIX_COLUMN; STD_FINAL_PNL_MATRIX(:,j)];
    EXEC_RATIO_MATRIX_COLUMN = [EXEC_RATIO_MATRIX_COLUMN; EXEC_RATIO_MATRIX(:,j)];

end


figure(1)
surf(dSmin_vector, dtau_vector_hour, E_FINAL_PNL_MATRIX)
axis([dSmin_vector(1) dSmin_vector(end) dtau_vector_hour(1) dtau_vector_hour(end) 0.5*min(min(E_FINAL_PNL_MATRIX)) 1.5*max(max(E_FINAL_PNL_MATRIX))])
str1 = ['Exp. Final P&L for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' S_0 = ' num2str(S0) ' dt = ' num2str(dt) ];
title(str1);
xlabel('dS_m_i_n ($)'); ylabel('\Deltat_o_r_d_e_r (hour)');
    
figure(2)
surf(dSmin_vector, dtau_vector_hour, STD_FINAL_PNL_MATRIX)
axis([dSmin_vector(1) dSmin_vector(end) dtau_vector_hour(1) dtau_vector_hour(end) 0.9*min(min(STD_FINAL_PNL_MATRIX)) 1.1*max(max(STD_FINAL_PNL_MATRIX))])
str1 = ['Standard Dev of Final P&L for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' S_0 = ' num2str(S0) ' dt = ' num2str(dt) ];
title(str1);
xlabel('dS_m_i_n ($)'); ylabel('\Deltat_o_r_d_e_r (hour)')
    
figure(3)
surf(dSmin_vector, dtau_vector_hour, CONDITIONAL_E_FINAL_LOSS_MATRIX)
axis([dSmin_vector(1) dSmin_vector(end) dtau_vector_hour(1) dtau_vector_hour(end) 0.9*min(min(CONDITIONAL_E_FINAL_LOSS_MATRIX)) 1.1*max(max(CONDITIONAL_E_FINAL_LOSS_MATRIX))])
str1 = ['Conditional Final Loss for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' S_0 = ' num2str(S0) ' dt = ' num2str(dt) ];
title(str1);
xlabel('dS_m_i_n ($)'); ylabel('\Deltat_o_r_d_e_r (hour)')

figure(4)
surf(dSmin_vector, dtau_vector_hour, E_INTER_HEDGING_MATRIX_hours)
axis([dSmin_vector(1) dSmin_vector(end) dtau_vector_hour(1) dtau_vector_hour(end) 0.9*min(min(E_INTER_HEDGING_MATRIX_hours)) 1.1*max(max(E_INTER_HEDGING_MATRIX_hours))])
str1 = ['Exp. Inter-hedging time (hour) for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' S_0 = ' num2str(S0) ' dt = ' num2str(dt) ];
title(str1);
xlabel('dS_m_i_n ($)'); ylabel('\Deltat_o_r_d_e_r (hour)')

figure(5)
surf(dSmin_vector, dtau_vector_hour, EXEC_RATIO_MATRIX)
axis([dSmin_vector(1) dSmin_vector(end) dtau_vector_hour(1) dtau_vector_hour(end) 0.9*min(min(EXEC_RATIO_MATRIX)) 1.1*max(max(EXEC_RATIO_MATRIX))])
str1 = ['Re-hedging execution rate for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' S_0 = ' num2str(S0) ' dt = ' num2str(dt) ];
title(str1);
xlabel('dS_m_i_n ($)'); ylabel('\Deltat_o_r_d_e_r (hour)')

figure(6)
plot(EXEC_RATIO_MATRIX_COLUMN, STD_FINAL_PNL_MATRIX_COLUMN, '.r', 'MarkerSize',15); grid on;
str1 = ['P&L Standard Dev versus re-hedging rate for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' S_0 = ' num2str(S0) ' dt = ' num2str(dt) ];
title(str1); xlabel('Re-hedging rate'); ylabel('Terminal P&L Std Deviation');

toc;

end