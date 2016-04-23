



% decision variables
% x = [x11, x12, ..., x16, x21, ..., x96 | y5, y6, y7, y8 | z1, z2, z3]'

% initialize variables for Ax <= b
A = zeros(157,61);
b = NaN(157,1);

% construction matrix for x_it location indices
x_it = zeros(6,9);
for i = 1:numel(x_it)
    x_it(i) = i; 
end

% construct vector for y_i location indices
y_i = NaN(1,9);
y_i(5:8) = [55 56 57 58];

% construct vector for z_k location indices
z_k = [59 60 61];

% given
l = [1390 1005 930 845 455 260]; % CHECK IF TRANSPOSE
s_kt_max = [ 80  35  20  15   5 0;
             70  30  15  10.5 5 0;
             100 0   0   0    0 0];
o_i_u = [0.05 0.04 0.035 0.05 0.05 0.05 0.035 0.03 0.10]';
o_i_p = [0.15 0.15 0.50 0.10 0.15 0.13 0.70 0.8 0.15]';
g_i = [1020 766 435 16 392 305 15 45 80]';
x_i_max = [240 210 180 400 500 500 320 400 200]';
x_i_min = [70 0 0 100 0 0 0 0 0]';
n_t = [106 310 670 1250 2780 3650];
c_i_v = [25.8 29.1 79.8 11.5 36.9 44.3 0.95 1.0 1.5]';
sigma_i_v = [3.0 3.0 3.0 4.0 2.0 2.0 0.1 0.5 0.25]';
c_i_c = [0 0 0 0 93.4 17.8 83.4 145.0 0]'; % CHECK IF 0 OR NAN
sigma_i_c = [0 0 0 0 9 4 18 20 0]'; % CHECK IF 0 OR NAN
c_k_d = [55 70 100]';
sigma_k_d = [15 4 20]';

         
% building matrix A
% load constraints
 
i=1;
for ct = 1:6
    % load constraints
     for t = 1:6
         A(ct,x_it(i,t))=-1;
         for k = 1:3
            A(ct, z_k(k)) = -s_kt_max(k,ct);
         end
         b(ct)= -l(ct);
     end
     i=i+1;
end

% instantaneous capacity constraint
ct = 7;
A(7:60,1:54)=eye(54);
for i = 1:9
    b(ct+6*(i-1):ct+6*i)=(1-o_i_u(i))*x_i_max(i); 
end

% annual energy constraints
ct = 60;
for i = 1:4
    A(ct+i,x_it(:,i)') = n_t;
    b(ct+i) = (1-o_i_p(i))*x_i_max(i)*8766;
end
% index 5 refers to row location in A, 9 refers to power plant 9 in x_it
% TODO: FIX
A(ct+5,x_it(:,9)') = n_t;
b(ct+5) = (1-o_i_p(9))*x_i_max(9)*8766;

% constraints on annual energy for new plants
ct = 65;
for i = 5:8
    A(ct+i-4,x_it(:,i)') = n_t;
    A(ct+i-4,y_i(i)) = -(1-o_i_p(i));
    b(ct+i-4) = 0;
end

% minimum generation constraint
ct = 69;
A(ct+1:ct+54,1:54) = -eye(54,54);
for i = 1:9
    b(ct+1+6*(i-1):ct+6*i)=-x_i_min(i); 
end

% constraints on dsm lower bound
ct = 123;
A(ct+1:ct+3,z_k) = -eye(3,3)';
for i = 1:3
    b(ct+i) = 0;
end


% constraints on dsm upper bound
ct = 126;
A(ct+1:ct+3,z_k) = eye(3,3)';
for i = 1:3
    b(ct+i) = 1;
end

% constraints on new generation bounds
ct = 129;
% TODO: FIX, 24 is (5,6,7,8)x(1:6)
A(ct+1:ct+24,x_it(1,5):x_it(6,8)) = eye(24,24);
% 6 corresponds to load blocks
for i = 5:8
    A(ct+1:ct+6,y_i(i)) = -(1-o_i_u(i));
    b(ct+1:ct+6) = 0;
    ct=ct+6;
end


ct = 153;
A(ct+1:ct+4,55:58) = eye(4,4);
for i = 5:8
    b(ct+i-4)=x_i_max(i); 
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
[X, FVAL, EXITFLAG] = linprog(f,A,b);

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