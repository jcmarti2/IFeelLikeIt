classdef IIM < handle
    
    % IIM 
    % description
    %{
        this class is an implementation of the Inoperability Input-Output
        model framework for dependency analysis
        
        authors:
            Juan Carlos Martinez
            Paul Gharzouzi
            Jimmy Chang
    %}
    
    % methods
    %{   
        PUBLIC:
            IIM(loc_A): constructor
            get_A(): return matrix A
            get_S(): return matrix S
            get_x(f): return vector x
            get_gamma_i(i): return gamma for row i
            get_delta_j(j): return delta for column j
            get_gamma_bar_i(i): return gamma bar for row i
            get_delta_bar_j(j): return delta bar for column j
            get_damage_propagation(f,k): return array holding damage
                                         propagations for k steps
            simulate_monte_carlo_A(f,k,low,up)
            simulate_monte_carlo_f(f,k,low,up)
            get_fmax_l1(b)
            get_fmax_l2(b)
    %}
    
    properties (Access = private)
        A
        S
    end
    
    methods(Access = public)
        
        function obj = IIM(loc_A)
            %{
                class constructor
                :param loc_A: address of .csv file holding A
                :return:
            %}
           obj.A = csvread(loc_A);
           obj.S = inv(eye(size(obj.A)) - obj.A);
        end
        
        function A = get_A(obj)
            %{
                get matrix A of IIM class instance
                :param:
                :return A: matrix A of IIM class instance
            %}
            A = obj.A;
        end
        
        function S = get_S(obj)
            %{
                get matrix S of IIM class instance
                :param:
                :return S: matrix S of IIM class instance
            %}
            S = obj.S;
        end
        
        function x = get_x(obj, f)
            %{
                get vector x given vector f
                :param f: external disturbance vector
                :return x: degradation vector
            %}
            x = obj.compute_x(f);
        end
        
        function gamma_i = get_gamma_i(obj, i)
            %{
                get gamma_i from matrix A at row index i
                :param i: row index for gamma computation
                :return gamma_i: gamma from row index i
            %}
            gamma_i = obj.compute_row_avg(obj.A,i);
        end
        
        function delta_j = get_delta_j(obj, j)
            %{
                get delta_j from matrix A at column index j
                :param j: column index for delta computation
                :return delta_j: delta from column index j
            %}
            delta_j = obj.compute_col_avg(obj.A,j);
        end
        
        function gamma_bar_i = get_gamma_bar_i(obj, i)
            %{
                get gamma_i from matrix S at row index i
                :param i: row index for gamma computation
                :return gamma_bar_i: gamma from row index i
            %}
            gamma_bar_i = obj.compute_row_avg(obj.S,i);
        end
        
        function delta_bar_j = get_delta_bar_j(obj, j)
            %{
                get delta_j from matrix S at column index j
                :param j: column index for delta computation
                :return delta_bar_j: delta from column index j
            %}
            delta_bar_j = obj.compute_col_avg(obj.S,j);
        end
        
        function modify_A(obj, changes)
            %{
                modify obj.A according to specified parameters and 
                compute new obj.S
                :param changes: matrix holding changes from original A in
                                percent change per element
                :return:
            %}
            obj.reconstruct_IIM(changes);
        end
        
        function x_k = get_damage_propagation(obj, f, k)
            %{
                get damage propagation x_k for k propagation steps
                :param f: external disturbance vector
                :param k: number of propagation steps
                :return x_k: array holding damage propagations
            %}
            x_k = obj.compute_damage_propagation(f, k);
        end
        
        function mc_k = simulate_monte_carlo_A(obj, f, k, low, up)
            %{
                perform monte carlo simulation with uncertainity on A
                :param f: external disturbance vector
                :param k: number of simulation replications
                :param low: lower percentage bound (from a given element of
                            A) for the uniform distriubtion of noise 
                :param up: upper percentage bound (from a given element of
                           A) for the uniform distriubtion of noise 
                :return mc_k: array holding simulation results
            %}
            mc_k = obj.replicate_monte_carlo_A(f, k, low, up);
        end
              
        function mc_k = simulate_monte_carlo_f(obj, f, k, low, up)
            %{
                perform monte carlo simulation with uncertainity on f
                :param f: external disturbance vector
                :param k: number of simulation replications
                :param low: lower percentage bound (from a given element of
                            f) for the uniform distriubtion of noise 
                :param up: upper percentage bound (from a given element of
                           f) for the uniform distriubtion of noise 
                :return mc_k: array holding simulation results
            %}
            mc_k = obj.replicate_monte_carlo_f(f, k, low, up);
        end
        
        function [fmax,l,exitflag,output,lambda] = get_fmax_l1(obj,b)
            %{
                solve optimization problem solving for vector fmax that has
                the maximum l1 norm
                :param b: constraint vector on x(in the form Sf <= b)
                          this parameter is used for function quadprog
                :return fmax: fmax vector (max overall damage vector)
                :return l: l1 norm of fmax
                :return exitflag: exit flag
                :return output: output information from solver
                :return lambda: lagrange multipliers
            %}
            [fmax,l,exitflag,output,lambda] = obj.compute_fmax(b, 1);
        end
        
        function [fmax,l,exitflag,output,lambda] = get_fmax_l2(obj,b)
            %{
                solve optimization problem solving for vector fmax that has
                the maximum l2 norm
                :param b: constraint vector on x(in the form Sf <= b)
                          this parameter is used for function quadprog
                :return fmax: fmax vector (max overall damage vector)
                :return l: l2 norm of fmax
                :return exitflag: exit flag
                :return output: output information from solver
                :return lambda: lagrange multipliers
            %}
            [fmax,l,exitflag,output,lambda] = obj.compute_fmax(b, 2);
        end
        
    end
        
    methods(Access = private)
        
        function x = compute_x(obj, f)
            %{
                get vector x given vector f
                :param f: external disturbance vector
                :return x: degradation vector
            %}
            x = obj.S*f;
        end
        
        function row_avg = compute_row_avg(~, M, i)
            %{
                compute the average of row i of matrix M
                :param M: matrix from which row average will be computed
                :param i: index for row average
                :return row_avg: row average from matrix M at row i
            %}
            vec = [M(i,1:i-1), M(i,i+1:end)];
            row_avg = sum(vec)/length(vec);
        end
        
        function col_avg = compute_col_avg(~, M, j)
            %{
                compute the average of column j of matrix M
                :param M: matrix from which column average will be computed
                :param j: index for column average
                :return col_avg: column average from matrix M at column j
            %}
            vec = [M(1:j-1,j); M(j+1:end,j)];
            col_avg = sum(vec)/length(vec);
        end
        
        function reconstruct_IIM(obj, changes)
            %{
                change obj.A and obj.S by adding specified changes to obj.A
                :param changes: matrix holding changes from original A in
                                percent change per element
                :return:
            %}
           changes = changes./100;
           obj.A = obj.A + obj.A.*changes;
           obj.S = inv(eye(size(obj.A)) - obj.A);
        end
        
        function x_k = compute_damage_propagation(obj, f, k)
            %{
                get damage propagation x_k for k propagation steps
                :param f: external disturbance vector
                :param k: number of propagation steps
                :return x_k: array holding damage propagations
            %}
            % set initial step and initialize storage structure
            i = 1;
            x_k = zeros(length(obj.A),k+1);
            % call recursive helper function
            x_k = propagate_damage(obj, x_k, i+1, f, k);
        end
        
        function x_k = propagate_damage(obj, x_k, i, f, k)
            %{
                recursive function to propagate damage to the kth step
                :param x_k: structure were damage propagation is stored
                :param i: current propagation step
                :param f: external disturbance vector
                :param k: number of propagation steps
                :return x_k: structure were damage propagation is stored
            %}
            % base case
            x_k(:,i) = obj.A*x_k(:,i-1) + f;
            % recursive step
            if i <= k
                x_k = propagate_damage(obj, x_k, i+1, f, k);
            end
        end
        
        function mc_k = replicate_monte_carlo_A(obj, f, k, low, up)
            %{
                perform monte carlo replications with uncertainity on A
                :param f: external disturbance vector
                :param k: number of simulation replications
                :param low: lower percentage bound (from a given element of
                            A) for the uniform distriubtion of noise 
                :param up: upper percentage bound (from a given element of
                           A) for the uniform distriubtion of noise 
                :return mc_k: array holding simulation results
            %}
            mc_k = zeros(length(obj.A),k);
            for i = 1:k
                Ai = obj.add_noise(obj.A,low,up);
                Si = inv(eye(size(Ai)) - Ai);
                mc_k(:,i) = Si*f;   
            end
        end
        
        function mc_k = replicate_monte_carlo_f(obj, f, k, low, up)
            %{
                perform monte carlo replications with uncertainity on f
                :param f: external disturbance vector
                :param k: number of simulation replications
                :param low: lower percentage bound (from a given element of
                            f) for the uniform distriubtion of noise 
                :param up: upper percentage bound (from a given element of
                           f) for the uniform distriubtion of noise 
                :return mc_k: array holding simulation results
            %}
            mc_k = zeros(length(obj.A),k);
            for i = 1:k
                fi = obj.add_noise(f,low,up);
                mc_k(:,i) = obj.S*fi;   
            end
        end
        
        function M = add_noise(~, M, low, up)
            %{
                add uniformly distributed noise to matrix M
                :param M: matrix which will have noise added
                :param low: lower percentage bound (from a given element of
                            M) for the uniform distriubtion of noise
                :param up: upper percentage bound (from a given element of
                            M) for the uniform distribution of noise
                :return M: matrix with added uniformly distributed rnoise
            %}
            low = low/100;
            up = up/100;
            M = (up - low).*M.*rand(size(M)) + (1+low).*M;  
        end
        
        function [fmax,l,exitflag,output,lambda] = compute_fmax(obj, b, type)
            %{
                solve optimization problem solving for vector fmax that has
                the maximum l1 norm
                :param b: constraint vector on x(in the form Sf <= b)
                          this parameter is used for function quadprog
                :param type: type of norm needed (1 for l1, 2 for l2)
                :return fmax: fmax vector (max overall damage vector)
                :return l: norm of fmax
                :return exitflag: exit flag
                :return output: output information from solver
                :return lambda: lagrange multipliers
            %}
            len = length(obj.A);
            % check for norm type and use matlab built in solver
            if type == 1
                % compute linprog parameters and call function
                % note that the variables f and A are NOT the properties 
                % of the class, but rather the input parameter for helper
                f = -ones(len,1);
                A = obj.S;
                Aeq = [];
                beq = [];
                lb = zeros(len,1);
                ub = ones(len,1);
                [fmax,l,exitflag,output,lambda] = linprog(f, A, b, Aeq, ...
                    beq, lb, ub);
            elseif type == 2
                % compute quadprog parameters and call function
                % note that the variables f and A are NOT the properties 
                % of the class, but rather the input parameter for helper
                H = -1*eye(len,len);
                f = zeros(len,1);
                A = obj.S;
                Aeq = [];
                beq = [];
                lb = zeros(len,1);
                ub = ones(len,1);
                [fmax,l,exitflag,output,lambda] = quadprog(H, f, A, b, ...
                    Aeq, beq, lb, ub);
            end
        end
        
    end
    
end

