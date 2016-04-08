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
q_2_9 = false;
q_2_10 = false;
q_2_12 = false; % TODO (ALSO IN IIM CLASS)
q_2_13 = false; % TODO PLOTS AND ANALYSIS
q_2_15 = false; % TODO FIX PLOTS AND ANALYSIS
q_2_16 = false; % TODO PLOTS AND ANALYSIS
q_2_17 = false; % TODO PLOTS AND ANALYSIS
q_2_18 = false; % TODO PLOTS AND ANALYSIS
q_2_19 = false; % TODO (ALSO IN IIM CLASS)


% =================================
%   NO USER INPUT AFTER THIS LINE
% =================================

IIM = IIM(file);

% 2.9
if q_2_9
    len = 10;
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
    figure()
    bar(delta)
    title('Infrastructure Influence Indices ( $\delta_{j}$ )','Interpreter','latex')
    ylabel('Index ( $\delta_{j}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($j$)','Interpreter','latex')
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
    figure()
    bar(delta_bar)
    title('Overall Infrastructure Influence Indices ( $\bar{\delta}$ $_{j}$ )','Interpreter','latex')
    ylabel('Overall Index ( $\delta_{j}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($j$)','Interpreter','latex')
end


% 2.12
if q_2_12
    % TODO ADD FUNCTION IN CLASS THAT CHANGES ELEMENTS OF A AND COMPUTES 
    % NEW S
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
    figure()
    bar(delta_bar)
    title('Overall Infrastructure Influence Indices ( $\bar{\delta}$ $_{j}$ )','Interpreter','latex')
    ylabel('Overall Index ( $\delta_{j}$ )','Interpreter','latex')
    xlabel('Infrastructure Sector ($j$)','Interpreter','latex')
end

% 2.13
if q_2_13
    f_213 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    x_213 = IIM.get_x(f_213);
end

%2.15
if q_2_15
    f_215 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    k_215 = 10;
    x_k_215 = IIM.get_damage_propagation(f_215,k_215);
    num_fig = 0;
    num_k_per_fig = 5;
    for i = 1:len
        if mod(i,num_k_per_fig) == 1
            figure()
            num_fig = num_fig + 1;
        end
        subplot(num_k_per_fig,1,i - num_k_per_fig*(num_fig - 1))
        plot(x_k_215(i,:))
    end
end

%2.16
if q_2_16
    f_216 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    k_216 = 10;
    low_216 = -15;
    up_216 = 15;
    mc_k = IIM.simulate_monte_carlo_A(f_216, k_216, low_216, up_216);
end

%2.17
if q_2_17
    f_217 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    k_217 = 10;
    low_217 = -15;
    up_217 = 15;
    mc_k = IIM.simulate_monte_carlo_f(f_217, k_217, low_217, up_217);
end

% 2.18
if q_2_18
    f_218 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
    k_218 = 10;
    low_218 = -15;
    up_218 = 15;
    mc_k_A_low = IIM.simulate_monte_carlo_A(f_218, k_218, low_218, low_218);
    mc_k_f_low = IIM.simulate_monte_carlo_f(f_218, k_218, low_218, low_218);
    mc_k_A_up = IIM.simulate_monte_carlo_A(f_218, k_218, up_218, up_218);
    mc_k_f_up = IIM.simulate_monte_carlo_f(f_218, k_218, up_218, up_218);
end