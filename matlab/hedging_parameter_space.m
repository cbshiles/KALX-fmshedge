function [E_FINAL_PNL_MATRIX STD_FINAL_PNL_MATRIX CONDITIONAL_E_FINAL_LOSS_MATRIX] = hedging_parameter_space(S0, K, T, R, SIGMA, npaths)

%% Scale Conversions

Conv_Year2Hour = 24*12*21;
Conv_Year2Minute = 24*12*21*60;

%% Parameter Space Plotting

dt_min = 0.5/(24*12*21);
dSmin_min = 0.05;

dt_vector = dt_min*[1 2 3 4 6 8 12];
dSmin_vector = dSmin_min:dSmin_min:20*dSmin_min;

E_FINAL_PNL_MATRIX = zeros(length(dt_vector),length(dSmin_vector));
STD_FINAL_PNL_MATRIX = zeros(length(dt_vector),length(dSmin_vector));
CONDITIONAL_E_FINAL_LOSS_MATRIX = zeros(length(dt_vector),length(dSmin_vector));


    for (j=1:1:length(dSmin_vector)) 
        for (i=1:1:length(dt_vector))
            
            [S CALL_PRICE UPPER_ORDER LOWER_ORDER UPPER_REHEDGE LOWER_REHEDGE QUANTITY_S QUANTITY_BOND MTM_S MTM_BOND PNL E_FINAL_PNL_MATRIX(i,j) STD_FINAL_PNL_MATRIX(i,j) CONDITIONAL_E_FINAL_LOSS_MATRIX(i,j)] = montecarlo_hedging(S0, K, dt_vector(i), T, R, SIGMA, npaths, dSmin_vector(j), 2);
       
        end
    
    
    end

dt_vector_hour = dt_vector.*(24*12*21);    
        
figure(1)
surf(dSmin_vector, dt_vector_hour, E_FINAL_PNL_MATRIX)
axis([dSmin_vector(1) dSmin_vector(end) dt_vector_hour(1) dt_vector_hour(end) 0.5*min(min(E_FINAL_PNL_MATRIX)) 1.5*max(max(E_FINAL_PNL_MATRIX))])
str1 = ['Exp. Final P&L for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' and S_0 = ' num2str(S0) ];
title(str1);
xlabel('dS_m_i_n ($)'); ylabel('dt (hour)');
    
figure(2)
surf(dSmin_vector, dt_vector_hour, STD_FINAL_PNL_MATRIX)
axis([dSmin_vector(1) dSmin_vector(end) dt_vector_hour(1) dt_vector_hour(end) 0.9*min(min(STD_FINAL_PNL_MATRIX)) 1.1*max(max(STD_FINAL_PNL_MATRIX))])
str1 = ['Range of  Final P&L for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' and S_0 = ' num2str(S0) ];
title(str1);
xlabel('dS_m_i_n ($)'); ylabel('dt (hour)')
    
figure(3)
surf(dSmin_vector, dt_vector_hour, CONDITIONAL_E_FINAL_LOSS_MATRIX)
axis([dSmin_vector(1) dSmin_vector(end) dt_vector_hour(1) dt_vector_hour(end) 0.9*min(min(CONDITIONAL_E_FINAL_LOSS_MATRIX)) 1.1*max(max(CONDITIONAL_E_FINAL_LOSS_MATRIX))])
str1 = ['Conditional Final Loss for \sigma = ' num2str(SIGMA) ' K = ' num2str(K) ' and S_0 = ' num2str(S0) ];
title(str1);
xlabel('dS_m_i_n ($)'); ylabel('dt (hour)')
    

end