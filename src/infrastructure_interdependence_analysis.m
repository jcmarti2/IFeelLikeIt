clear all
close all

len = 10;
file = 'A.csv';

IIM = IIM(file);

% 2.9
gamma = zeros(len,1);
for i = 1:len
    gamma(i) = IIM.get_gamma_i(i);
end
delta = zeros(len,1);
for j = 1:len
    delta(j) = IIM.get_delta_j(j);
end

% 2.10
gamma_bar = zeros(len,1);
for i = 1:len
    gamma_bar(i) = IIM.get_gamma_bar_i(i);
end
delta_bar = zeros(len,1);
for j = 1:len
    delta_bar(j) = IIM.get_delta_bar_j(j);
end

% 2.12

% 2.13
f_213 = [0.35 0.1 0 0 0 0.5 0 0 0 0]';
x_213 = IIM.get_x(f_213);

% 2.18
IIM.simulate_monte_carlo_A(f, 10, -15, -15)
IIM.simulate_monte_carlo_f(f, 10, -15, -15)

IIM.simulate_monte_carlo_A(f, 10, 15, 15)
IIM.simulate_monte_carlo_f(f, 10, 15, 15)