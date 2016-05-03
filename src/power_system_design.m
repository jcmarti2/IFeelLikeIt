%{      
 
this script computes the answers for the questions of part 4 of the
power system design model

to use the script, fill in the user input
    - set true/false for answering each question in the variables
      q_2_#, where # is the question number

authors:
        Juan Carlos Martinez
        Paul Gharzouzi
        Jimmy Chang
%}
close all
clear all

% =================================
%             USER INPUT
% =================================

q_4_28 = false;
q_4_29 = false;
q_4_30 = false;
q_4_31 = false;
q_4_32 = false;
q_4_33 = false;
q_4_34 = false;
q_4_35 = true;

% =================================
%   NO USER INPUT AFTER THIS LINE
% =================================

% Given parameters
x_i_max = [240 210 180 400 500 500 320 400 200]';
x_i_min = [70 0 0 100 0 0 0 0 0]';
o_i_u = [0.05 0.04 0.035 0.05 0.05 0.05 0.035 0.03 0.10]';
o_i_p = [0.15 0.15 0.50 0.10 0.15 0.13 0.70 0.8 0.15]';
g_i = [1020 766 435 16 392 305 15 45 80]';

c_i_v = [25.8 29.1 79.8 11.5 36.9 44.3 0.95 1.0 1.5]';
sigma_i_v = [3.0 3.0 3.0 4.0 2.0 2.0 0.1 0.5 0.25]';
c_i_c = [0 0 0 0 93.4 17.8 83.4 145.0 0]';
sigma_i_c = [0 0 0 0 9 4 18 20 0]';

n_t = [106 310 670 1250 2780 3650];
l = [1390 1005 930 845 455 260];

s_kt_max = [ 80  35 20 15   5 0;
             70  30 15 10.5 5 0;
             100 0  0  0    0 0];
c_k_d = [55 70 100]';
sigma_k_d = [15 4 20]';

% class instance
power_system = power_system(x_i_max, x_i_min, o_i_u, o_i_p, g_i, c_i_v, ...
    sigma_i_v, c_i_c, sigma_i_c, n_t, l, s_kt_max, c_k_d, sigma_k_d);

% 4.28
if q_4_28
    [x, cost, emissions, variance, slack, binding, ...
                shadow_cost, shadow_emi, shadow_var, exitflag, output, ...
                lambda] = power_system.minimize_cost();
    fprintf('COST MINIMIZATION \n')
    fprintf('Cost [$]: %i.\n',cost)
    fprintf('Emissions [tons]: %i.\n',emissions)
    fprintf('Cost standard deviation [$]: %i.\n',sqrt(variance))
    fprintf('\n')
end

% 4.29
if q_4_29
    [x, cost, emissions, variance, slack, binding, ...
                shadow_cost, shadow_emi, shadow_var, exitflag, output, ...
                lambda] = power_system.minimize_emissions();
    fprintf('EMISSIONS MINIMIZATION \n')
    fprintf('Cost [$]: %i.\n',cost)
    fprintf('Emissions [tons]: %i.\n',emissions)
    fprintf('Cost standard deviation [$]: %i.\n',sqrt(variance))
    fprintf('\n')
end

% 4.30
if q_4_30
    [costs, emissions, decisions] = power_system.minimize_cost_and_emissions();
    plot(costs,emissions)
    title('Cost-Emissions Pareto Optimality','Interpreter','latex')
    xlabel('Cost [\$]','Interpreter','latex')
    ylabel('Emissions [tons $CO_{2}$]','Interpreter','latex')
    grid on
end

% 4.31
if q_4_31
    [costs, caps] = power_system.minimize_cost_with_emissions_cap();
    plot(costs,caps)
    xlabel('Emissions cap [tons CO2]')
    ylabel('Costs [$]')
    grid on
end

% 4.32
if q_4_32
    [rates,cost,emissions] = power_system.tax_emissions();
    figure
    plot(cost,emissions)
    title('Cost-Emissions Curve after Tax Implementation','Interpreter','latex')
    xlabel('Costs [\$]','Interpreter','latex')
    ylabel('Emissions [tons $CO_{2}$]','Interpreter','latex')
    grid on
    figure
    hold on
    [ax, h1, h2] = plotyy(rates, cost, rates, emissions, 'plot');
    title('Cost and Emissions with Tax Implementation','Interpreter','latex')
    xlabel('Tax rate [\%]','Interpreter','latex')
    set(get(ax(1), 'Ylabel'), 'String', 'Cost [\$]','Interpreter','latex');
    set(get(ax(2), 'Ylabel'), 'String', 'Emissions [tons $CO_{2}$]','Interpreter','latex');
    legend('Costs','Emissions')
    grid on
end

% 4.33
if q_4_33
    
    [x, cost, emissions, variance, slack, binding, ...
                shadow_cost, shadow_emi, shadow_var, exitflag, output, ...
                lambda] = power_system.minimize_variance();
    fprintf('VARIANCE MINIMIZATION \n')
    fprintf('Cost [$]: %i.\n',cost)
    fprintf('Emissions [tons]: %i.\n',emissions)
    fprintf('Cost standard deviation [$]: %i.\n',sqrt(variance))
    fprintf('\n')
end

% 4.34
if q_4_34
    [costs, variances, decisions] = power_system.minimize_cost_and_variance();
    plot(costs,variances)
    title('Cost-Variance Pareto Optimality','Interpreter','latex')
    xlabel('Cost [\$]','Interpreter','latex')
    ylabel('Variance [$\$^2$]','Interpreter','latex')
    grid on
end

% 4.35
if q_4_35
    [costs, emissions, variances, decisions] = ...
                power_system.minimize_cost_emissions_and_variance();
    scatter3(costs, emissions, variances,'.')
    title('Cost-Emissions-Variance Pareto Optimality','Interpreter','latex')
    xlabel('Cost [\$]','Interpreter','latex')
    ylabel('Emissions [tons $CO_{2}$]','Interpreter','latex')
    zlabel('Variance [$\$^2$]','Interpreter','latex')
end
