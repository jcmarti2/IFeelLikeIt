%{      
 
this script computes the answers for the questions of part 2 of the
inoperability input-output model

to use the script, fill in the user input
    - set the file address (location of matrix A)
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

file = '../data/A.csv';
len = 10;
q_2_9 = false;
q_2_10 = false;
q_2_12 = true;
q_2_13 = false;
q_2_15 = false;
q_2_16 = false;
q_2_17 = false; 
q_2_18 = false; 
q_2_19 = false; % TODO (ALSO IN IIM CLASS)

% WARNING: 
% Do not run q_2_12 in combination with other questions.
% Matrix A is changed, so it will change the results of the
% instance. Run separately.

% =================================
%   NO USER INPUT AFTER THIS LINE
% =================================

IIM = IIM(file);

% 2.9
if q_2_9
    gamma = zeros(len,1);
    for i = 1:len
        gamma(i) = IIM.get_gamma_i(i);
    end
    delta = zeros(len,1);
    for j = 1:len
        delta(j) = IIM.get_delta_j(j);
    end
    figure()
    bar(gamma)
    title('Infrastructure Dependency Indices ( $\gamma_{i}$ )','Interpreter','latex')
    ylabel('Index ( $\gamma_{i}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($i$)','Interpreter','latex')
    grid on
    figure()
    bar(delta)
    title('Infrastructure Influence Indices ( $\delta_{j}$ )','Interpreter','latex')
    ylabel('Index ( $\delta_{j}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($j$)','Interpreter','latex')
    grid on
end

% 2.10
if q_2_10
    gamma_bar = zeros(len,1);
    for i = 1:len
        gamma_bar(i) = IIM.get_gamma_bar_i(i);
    end
    delta_bar = zeros(len,1);
    for j = 1:len
        delta_bar(j) = IIM.get_delta_bar_j(j);
    end
    figure()
    bar(gamma_bar)
    title('Overall Infrastructure Dependency Indices ( $\bar{\gamma}$ $_{i}$ )','Interpreter','latex')
    ylabel('Overall Index ( $\gamma_{i}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($i$)','Interpreter','latex')
    grid on
    figure()
    bar(delta_bar)
    title('Overall Infrastructure Influence Indices ( $\bar{\delta}$ $_{j}$ )','Interpreter','latex')
    ylabel('Overall Index ( $\delta_{j}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($j$)','Interpreter','latex')
    grid on
end


% 2.12
% WARNING: Do not run in combination with other questions.
% Matrix A is changed, so it will change the results of the
% instance. Run separately.
if q_2_12
    changes = zeros(len,len);
    changes(:,2) = [0 0 10 10 0 0 0 -10 -10 10]' ;
    
    gamma_bar_old = zeros(len,1);
    for i = 1:len
        gamma_bar_old(i) = IIM.get_gamma_bar_i(i);
    end
    delta_bar_old = zeros(len,1);
    for j = 1:len
        delta_bar_old(j) = IIM.get_delta_bar_j(j);
    end
    
    IIM.modify_A(changes);
    gamma_bar_new = zeros(len,1);
    for i = 1:len
        gamma_bar_new(i) = IIM.get_gamma_bar_i(i);
    end
    delta_bar_new = zeros(len,1);
    for j = 1:len
        delta_bar_new(j) = IIM.get_delta_bar_j(j);
    end
    
    figure()
    bar(gamma_bar_new)
    title('Overall Infrastructure Dependency Indices ( $\bar{\gamma}$ $_{i}$ )','Interpreter','latex')
    ylabel('Overall Index ( $\gamma_{i}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($i$)','Interpreter','latex')
    grid on
    figure()
    bar(delta_bar_new)
    title('Overall Infrastructure Influence Indices ( $\bar{\delta}$ $_{j}$ )','Interpreter','latex')
    ylabel('Overall Index ( $\delta_{j}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($j$)','Interpreter','latex')
    grid on
    
    change_gamma = 100*(gamma_bar_new - gamma_bar_old)./gamma_bar_old;
    change_delta = 100*(delta_bar_new - delta_bar_old)./delta_bar_old;
    figure()
    bar(change_gamma)
    title('Percent Change in Overall Infrastructure Dependency Indices ( $\bar{\gamma}$ $_{i}$ )','Interpreter','latex')
    ylabel('Percent Change in Overall Index ( $\gamma_{i}$ ) [\%]','Interpreter','latex')
    xlabel('Infrastructure Sector ($i$)','Interpreter','latex')
    grid on
    figure()
    bar(change_delta)
    title('Percent Change in Overall Infrastructure Influence Indices ( $\bar{\delta}$ $_{j}$ )','Interpreter','latex')
    ylabel('Percent Change in Overall Index ( $\delta_{j}$ ) [\%]','Interpreter','latex')
    xlabel('Infrastructure Sector ($j$)','Interpreter','latex')
    grid on
end

% 2.13
if q_2_13
    f_213 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    x_213 = IIM.get_x(f_213);
    bar(x_213)
    title('Snowstorm Effects','Interpreter','latex')
    ylabel('Overall Damage','Interpreter','latex')
    xlabel('Infrastructure Sector','Interpreter','latex')
    grid on
end

%2.15
if q_2_15
    len = 10;
    f_215 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    k_215 = 10;
    x_k_215 = IIM.get_damage_propagation(f_215,k_215);
    diff_x_k_215 = diff(x_k_215,1,2);
    
    % generate plot for damage after kth tier
    num_fig = 0;
    num_per_fig = 5;
    for i = 1:len
        if mod(i,num_per_fig) == 1
            figure()
            num_fig = num_fig + 1;
        end
        subplot(num_per_fig,1,i - num_per_fig*(num_fig - 1))
        plot(x_k_215(i,:))
        if mod(i,num_per_fig) == 1
            title('Propagation of the Initial Damage Through the $k^{th}$ Tier','Interpreter','latex')
        end
        if mod(i,num_per_fig) == 0 || i == len 
            xlabel('$k^{th}$ tier','Interpreter','latex')
        end
        YLim = get(gca,'YLim');
        ypos = sum(YLim)/2;
        pos = [0.6 ypos 0];
        ylabel(['Sector ',num2str(i)],'Interpreter','latex', 'position',pos)
        grid on
    end
    
    % generate plot for incremental damage after each tier
    num_fig = 0;
    num_per_fig = 5;
    for i = 1:len
        if mod(i,num_per_fig) == 1
            figure()
            num_fig = num_fig + 1;
        end
        subplot(num_per_fig,1,i - num_per_fig*(num_fig - 1))
        plot(diff_x_k_215(i,:))
        if mod(i,num_per_fig) == 1
            title('Incremental Propagation of Damage Through the $k^{th}$ Tier','Interpreter','latex')
        end
        if mod(i,num_per_fig) == 0 || i == len 
            xlabel('$k^{th}$ tier','Interpreter','latex')
        end
        YLim = get(gca,'YLim');
        ypos = sum(YLim)/2;
        pos = [0.6 ypos 0];
        ylabel(['Sector ',num2str(i)],'Interpreter','latex', 'position',pos)
        grid on
    end
end

%2.16
if q_2_16
    f_216 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    k_216 = 10000;
    low_216 = -15;
    up_216 = 15;
    mc_k = IIM.simulate_monte_carlo_A(f_216, k_216, low_216, up_216);
    
    % generate plot for histograms of each sector
    num_bins_216 = 100;
    num_fig_216 = 0;
    num_col_per_fig_216 = 2;
    num_per_fig_216 = 10; % must be an even number
    for i = 1:len
        if mod(i,num_per_fig_216) == 1
            figure()
            num_fig_216 = num_fig_216 + 1;
        end
        subplot(floor(num_per_fig_216/num_col_per_fig_216), ... 
            num_col_per_fig_216,i - num_per_fig_216*(num_fig_216 - 1))
        hist(mc_k(i,:),num_bins_216) 
        xlabel(['Sector ',num2str(i)],'Interpreter','latex')
        ylabel('Frequency','Interpreter','latex')
    end
    annotation('textbox', [0 0.9 1 0.1], ... 
               'String', 'Monte Carlo Simulation with Uncertainity in A', ...
               'EdgeColor', 'none', 'fontsize', 14, ...
               'HorizontalAlignment', 'center', 'interpreter', 'latex')
    figure()
    boxplot(mc_k')
    title('Monte Carlo Simulation with Uncertainty in A','interpreter','latex')
    xlabel('Infrastructure Sector','interpreter', 'latex')
    ylabel('Overall Damage','interpreter', 'latex')
    ylim([0 1])
    grid on
end

%2.17
if q_2_17
    f_217 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    k_217 = 10000;
    low_217 = -15;
    up_217 = 15;
    mc_k = IIM.simulate_monte_carlo_f(f_217, k_217, low_217, up_217);
    
    % generate plot for histograms of each sector
    num_bins_217 = 100;
    num_fig_217 = 0;
    num_col_per_fig_217 = 2;
    num_per_fig_217 = 10; % must be an even number
    for i = 1:len
        if mod(i,num_per_fig_217) == 1
            figure()
            num_fig_217 = num_fig_217 + 1;
        end
        subplot(floor(num_per_fig_217/num_col_per_fig_217), ... 
            num_col_per_fig_217,i - num_per_fig_217*(num_fig_217 - 1))
        hist(mc_k(i,:),num_bins_217)
        xlabel(['Sector ',num2str(i)],'Interpreter','latex')
        ylabel('Frequency','Interpreter','latex')
    end
    annotation('textbox', [0 0.9 1 0.1], ... 
               'String', 'Monte Carlo Simulation with Uncertainity in f', ...
               'EdgeColor', 'none', 'fontsize', 14, ...
               'HorizontalAlignment', 'center', 'interpreter', 'latex')
    figure()
    boxplot(mc_k')
    title('Monte Carlo Simulation with Uncertainty in f','interpreter','latex')
    xlabel('Infrastructure Sector','interpreter', 'latex')
    ylabel('Overall Damage','interpreter', 'latex')
    ylim([0 1])
    grid on
end

% 2.18
if q_2_18
    f_218 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    k_218 = 1;
    low_218 = -15;
    up_218 = 15;
    x_orig = IIM.get_x(f_218);
    mc_k_A_low = IIM.simulate_monte_carlo_A(f_218, k_218, low_218, low_218);
    mc_k_f_low = IIM.simulate_monte_carlo_f(f_218, k_218, low_218, low_218);
    mc_k_A_up = IIM.simulate_monte_carlo_A(f_218, k_218, up_218, up_218);
    mc_k_f_up = IIM.simulate_monte_carlo_f(f_218, k_218, up_218, up_218);
    figure()
    bar([mc_k_A_low x_orig mc_k_A_up])
    title('Damage Bounds with Uncertainty in A','Interpreter','latex')
    legend('Lower Bound (-15%)', 'Original Damage', 'Upper Bound (15%)')
    ylabel('Damage Coefficient','Interpreter','latex')
    xlabel('Infrastructure Sector','Interpreter','latex')
    grid on
    figure()
    bar([mc_k_f_low x_orig mc_k_f_up])
    title('Damage Bounds with Uncertainty in f','Interpreter','latex')
    legend('Lower Bound (-15%)', 'Original Damage', 'Upper Bound (15%)')
    ylabel('Damage Coefficient','Interpreter','latex')
    xlabel('Infrastructure Sector','Interpreter','latex')
    grid on
end