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


% Space to store A and b 
A = zeros(157,61);
b = NaN(157,1);

% Matrix to store x_it indices for each row of A
x_it = zeros(6,9); % CHANGE NAME LATER TO x_ti
for i = 1:numel(x_it)
    x_it(i) = i; 
end

% Vector to store y_i indices for each row of A
y_i = NaN(1,9)';
y_i(5:8) = [55 56 57 58]';

% Vector to store z_k indices for each row of A
z_k = [59 60 61];


% Building matrix A and vector b
idx = 1;

% Load constraint (6 rows)
% sum_i=1^9{x_it} + sum_k=1^3{s_kt_max} >= l_t for t={1,2,3,4,5,6}
for t = 1:6
    A(idx,x_it(t,:)) = -1;
    A(idx,z_k) = -s_kt_max(:,t)';
    b(idx) = -l(t);
    idx = idx + 1;
end

% Instantaneous capacity constraint (54 rows)
% x_it <= (1-o_i_u)x_i_max for i={1,2,3,4,5,6,7,8,9}, t={1,2,3,4,5,6}
for i = 1:9
    for t = 1:6
        A(idx,x_it(t,i)) = 1;
        b(idx) = (1-o_i_u(i))*x_i_max(i);
        idx = idx + 1; 
    end
end 

% Annual energy constraint for exisiting plants (5 rows)
% sum_t=1^6{n_tx_it} - (1-o_i_p)8766x_i_max <= for i={1,2,3,4,9}
for i = 1:4
    A(idx,x_it(:,i)') = n_t;
    b(idx) = (1-o_i_p(i))*x_i_max(i)*8766;
    idx = idx + 1;
end
A(idx,x_it(:,9)') = n_t;
b(idx) = (1-o_i_p(9))*x_i_max(9)*8766;
idx = idx + 1;

% Annual energy constraint for new plants (4 rows)
% sum_t=1^6{n_tx_it} - (1-o_i_p)8766y_i <= for i={5,6,7,8}
for i = 5:8
    A(idx,x_it(:,i)') = n_t;
    A(idx,y_i(i)) = -8766*(1-o_i_p(i));
    b(idx) = 0;
    idx = idx + 1;
end

% Minimum generation constraint (54 rows)
% x_it >= x_i_min for i={1,2,3,4,5,6,7,8,9}, t={1,2,3,4,5,6}
for i = 1:9
    for t = 1:6
        A(idx,x_it(t,i)) = -1;
        b(idx) = -x_i_min(i);
        idx = idx + 1;
    end
end 

% Lower bound on DSM program (3 rows)
% z_k >= 0 for k={1,2,3}
for k = 1:3
    A(idx,z_k(k)) = -1;
    b(idx) = 0;
    idx = idx + 1;
end

% Upper bound on DSM program (3 rows)
% z_k <= 1 for k={1,2,3}
for k = 1:3
    A(idx,z_k(k)) = 1;
    b(idx) = 1;
    idx = idx + 1;
end

% New generation bounds (24 rows)
% x_it <= (1-o_i_u)y_i for i={5,6,7,8}, t={1,2,3,4,5,6}
for i = 5:8
    for t = 1:6
        A(idx,x_it(t,i)) = 1;
        A(idx,y_i(i)) = -(1-o_i_u(i));
        b(idx) = 0;
        idx = idx + 1;
    end
end 

% New generation bounds (4 rows)
% y_i <= x_i_max for i={5,6,7,8}
for i = 5:8
    A(idx,y_i(i)) = 1;
    b(idx) = x_i_max(i);
    idx = idx + 1;
end




c_prime = [c_i_v' c_i_c(5:8)' c_k_d']';

H = zeros(16,61);
for i = 1:9
    H(i,1+6*(i-1):6+6*(i-1)) = n_t;
end
H(10:13,55:58) = 1000*eye(4,4);
var = 0;
for t = 1:6
    var = var + s_kt_max(1,t)*n_t(t);
end
H(14,59) = var;
var = 0;
for t = 1:6
    var = var + s_kt_max(2,t)*n_t(t);
end
H(15,60) = var;
var = 0;
for t = 1:6
    var = var + s_kt_max(3,t)*n_t(t);
end
H(16,61) = var;

f = c_prime'*H;
[X, FVAL, EXITFLAG] = linprog(-f,A,b);

%{

classdef power_system < handle
    
    % power_system 
    % description
    %{
        
        
        authors:
            Juan Carlos Martinez
            Paul Gharzouzi
            Jimmy Chang
    %}
    
    % methods
    %{   
        PUBLIC:
    %}
    
    properties (Access = private)
    end
    
    methods(Access = public)
        
        
       
        
        % optimize for expected cost
        % get cost given f
        % get emissions given f
        % perform cost sensitivity at f
        
        % optimize for emissions
        % get cost given f
        % get emissions given f
        % perform emissions sensitivity at f
        
        % perform multiobjective on expected cost and emissions
        
        % how to handle emissions cap? tax?
        
        % optimize for variance of cost
        % perform multiobjective on expected cost and variance of cost
        
        % optimize for variance of emissions
        % perform multiobjective on expected emissions and variance of
        % emissions
        
        
    end
        
    methods(Access = private)    
    end
    
end

%}