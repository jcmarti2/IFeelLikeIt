
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
    
    properties (Access = public)
        
        x_i_max
        x_i_min
        o_i_u
        o_i_p
        g_i
        c_i_v
        sigma_i_v
        c_i_c
        sigma_i_c
        n_t
        l
        s_kt_max
        c_k_d
        sigma_k_d
        
        A
        b
        
        x_ti
        y_i
        z_k
        
        cost_c_o
        cost_H
        cost_min_f
        cost_max_f
        cost_SIGMA
        
        emi_H
        emi_min_f
        emi_max_f
    end
    
    methods(Access = public)
        
        function obj = power_system(x_i_max, x_i_min, o_i_u, o_i_p, ...
                g_i, c_i_v, sigma_i_v, c_i_c, sigma_i_c, n_t, l, ... 
                s_kt_max, c_k_d, sigma_k_d)
            
            %{
                class constructor
                :param x_i_max: vector of plant maximum capacity
                :param x_i_min: vector of plant minimum capacity
                :param o_i_u: vector of unplanned outage rate
                :param o_i_p: vector of planned annual outage rate
                :param g_i: vector of global warming potential
                :param c_i_v: vector variable cost of plant operation
                :param sigma_i_v: vector of standard deviation of variable
                                  cost of plant operation
                :param c_i_c: vector of annaulized capital cost of plant
                :param sigma_i_c: vector of standard deviation of 
                                  annualized capital cost of plant
                :param n_t: array of number of hours in load block
                :param l: array of power load in load block
                :param s_kt_max: matrix of maximum energy saving for each
                                 load block and demand side managemenr
                                 program
                :param c_k_d: vector of cost of implementation of demand
                              side management program
                :param sigma_k_d: vector of cost standard deviation of 
                                  implementation of demand side management
                                  program
                :return:
            %}

            % Set given parameters
            obj.x_i_max = x_i_max;
            obj.x_i_min = x_i_min;
            obj.o_i_u = o_i_u;
            obj.o_i_p = o_i_p;
            obj.g_i = g_i;
            obj.c_i_v = c_i_v;
            obj.sigma_i_v = sigma_i_v;
            obj.c_i_c = c_i_c;
            obj.sigma_i_c = sigma_i_c;
            obj.n_t = n_t;
            obj.l = l;
            obj.s_kt_max = s_kt_max;
            obj.c_k_d = c_k_d;
            obj.sigma_k_d = sigma_k_d;
            
            % construct matrix A and vector b, corresponding to
            % constraints Ax<=b
            obj.construct_A_b()
            
            % construct cost objective functions (min and max)
            obj.gen_cost_obj_funs()
            
            % construct emissions objective functions (min and max)
            obj.gen_emi_obj_funs()
        end
           
        function [X, FVAL, VAR, SLACK, BINDING, SHADOW,  EXITFLAG, ... 
                OUTPUT, LAMBDA] = minimize_cost(obj)
            %{
                This function returns the linprog output for cost 
                minimization, as well as the cost variance and slack var
                :param:
                :return VAR: variance of the minimum cost
                :return SLACKX:
                :return X:
                :return FVAL:
                :return EXITFLAG:
                :return OUTPUT:
                :return LAMBDA:
            %}
                [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog( ...
                    obj.cost_min_f, obj.A, obj.b);
                SLACK = find(~(obj.A*X == obj.b));
                BINDING = find(obj.A*X == obj.b);
                VAR = X'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*X;
                
                SHADOW = zeros(numel(BINDING),1);
                for binding_idx = 1:numel(BINDING)
                    new_b = obj.b;
                    new_b(BINDING(binding_idx)) = ... 
                        new_b(BINDING(binding_idx)) + 1;
                    [~, new_fval] = linprog(obj.cost_min_f, obj.A, new_b);
                    SHADOW(binding_idx) = new_fval - FVAL;
                end
        end
        
        function [X, FVAL, VAR, SLACK, BINDING, EXITFLAG, OUTPUT, ...
                LAMBDA] = maximize_cost(obj)
            %{
                This function returns the linprog output for cost 
                minimization, as well as the cost variance and slack var
                :param:
                :return VAR: variance of the maximum cost
                :return SLACKX:
                :return X:
                :return FVAL:
                :return EXITFLAG:
                :return OUTPUT:
                :return LAMBDA:
            %}
                [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog( ...
                    obj.cost_max_f, obj.A, obj.b);
                SLACK = find(~(obj.A*X == obj.b));
                BINDING = find(obj.A*X == obj.b);
                VAR = X'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*X;     
        end
        
        function [X, FVAL, VAR, SLACK, BINDING, SHADOW, EXITFLAG, OUTPUT, ... 
                LAMBDA] = minimize_emissions(obj)
            %{
                :param: 
                :return SLACKX:
                :return X:
                :return FVAL:
                :return EXITFLAG:
                :return OUTPUT:
                :return LAMBDA:
            %}
                options = optimoptions('linprog','Algorithm', ...
                    'dual-simplex');
                [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog( ...
                    obj.emi_min_f, obj.A, obj.b,[],[],[],[],[],options);
                SLACK = find(~(obj.A*X == obj.b));
                BINDING = find(obj.A*X == obj.b);
                
                VAR = X'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*X;
                
                SHADOW = zeros(numel(BINDING),1);
                for binding_idx = 1:numel(BINDING)
                    new_b = obj.b;
                    new_b(BINDING(binding_idx)) = ... 
                        new_b(BINDING(binding_idx)) + 1;
                    [~, new_fval] = linprog(obj.emi_min_f, obj.A, ... 
                        new_b,[],[],[],[],[],options);
                    SHADOW(binding_idx) = new_fval - FVAL;
                end
        end
 
        function [X, FVAL, SLACK, BINDING, EXITFLAG, OUTPUT, LAMBDA] = ...
                maximize_emissions(obj)
            %{
                :param: 
                :return SLACKX:
                :return X:
                :return FVAL:
                :return EXITFLAG:
                :return OUTPUT:
                :return LAMBDA:
            %}
                [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = linprog( ...
                    obj.emi_max_f, obj.A, obj.b);
                SLACK = find(~(obj.A*X == obj.b));
                BINDING = find(obj.A*X == obj.b);
        end
        
        function [costs, emissions] = minimize_cost_and_emissions(obj)
            
            lambdas = linspace(0,1,5000);
            costs = zeros(1,numel(lambdas));
            emissions = zeros(1,numel(lambdas));
            [~,min_cost] = linprog(obj.cost_min_f, obj.A, obj.b);
            options = optimoptions('linprog','Algorithm', ...
                    'dual-simplex');
            [~,min_emi] = linprog(obj.emi_min_f, obj.A, obj.b,[],[], ...
                [],[],[],options);
            
            for idx = 1:numel(lambdas)
                lambda = lambdas(idx);
                f = (lambda/min_cost).*(obj.cost_c_o'*obj.cost_H) + ...
                    ((1-lambda)/min_emi).*(obj.g_i'*obj.emi_H);
                X = linprog(f, obj.A, obj.b,[],[],[],[],[],options);
                costs(idx) = obj.get_cost(X);
                emissions(idx) = obj.get_emissions(X);
            end
        end
        
        function cost = get_cost(obj, X)
            %{
                this function returns the cost given the decision variables
                :param X: vector of decision variables
                :return cost: total cost
            %}
                cost = obj.cost_c_o'*obj.cost_H*X;
        end
        
        function emissions = get_emissions(obj, X)
            %{
                this function returns the emissions given the decision
                variables
                :param X: vector of decision variables
                :return emissions: total emissions
            %}
            emissions = obj.g_i'*obj.emi_H*X;
        end
                
    end
        
    methods(Access = private)
                
        function construct_A_b(obj)
            %{
                this function constructs matrix A for solving the
                inequality constraints Ax<=b (A is 157x61, b is 1x157)
                Properties modified:
                    - A
                    - b
                :param:
                :return:
            %}
            
            % space to store A and b 
            obj.A = zeros(157,61);
            obj.b = NaN(157,1);
            
            % matrix to store x_it indices for each row in A 
            obj.x_ti = zeros(6,9);
            for i = 1:numel(obj.x_ti)
                obj.x_ti(i) = i; 
            end    
            % vector to store y_i indices for each row in A
            obj.y_i = NaN(1,9)';
            obj.y_i(5:8) = [55 56 57 58]';
            % vector to store z_k indices for each row in A
            obj.z_k = [59 60 61];
            
            % building matrix A and vector b
            idx = 1;
            % CONSTRAINT SET 1
            % load constraint (6 rows)
            % sum_i=1^9{x_it} + sum_k=1^3{s_kt_max} >= l_t for t={1,2,3,4,5,6}
            for t = 1:6
                obj.A(idx,obj.x_ti(t,:)) = -1;
                obj.A(idx,obj.z_k) = -obj.s_kt_max(:,t)';
                obj.b(idx) = -obj.l(t);
                idx = idx + 1;
            end
            
            % CONSTRAINT SET 2
            % instantaneous capacity constraint (54 rows)
            % x_it <= (1-o_i_u)x_i_max for i={1,2,3,4,5,6,7,8,9}, t={1,2,3,4,5,6}
            for i = 1:9
                for t = 1:6
                    obj.A(idx,obj.x_ti(t,i)) = 1;
                    obj.b(idx) = (1-obj.o_i_u(i))*obj.x_i_max(i);
                    idx = idx + 1; 
                end
            end 
            
            % CONSTRAINT SET 3
            % annual energy constraint for exisiting plants (5 rows)
            % sum_t=1^6{n_tx_it} - (1-o_i_p)8766x_i_max <= for i={1,2,3,4,9}
            for i = 1:4
                obj.A(idx,obj.x_ti(:,i)') = obj.n_t;
                obj.b(idx) = (1-obj.o_i_p(i))*obj.x_i_max(i)*8766;
                idx = idx + 1;
            end
            obj.A(idx,obj.x_ti(:,9)') = obj.n_t;
            obj.b(idx) = (1-obj.o_i_p(9))*obj.x_i_max(9)*8766;
            idx = idx + 1;
            
            % CONSTRAINT SET 4
            % annual energy constraint for new plants (4 rows)
            % sum_t=1^6{n_tx_it} - (1-o_i_p)8766y_i <= for i={5,6,7,8}
            for i = 5:8
                obj.A(idx,obj.x_ti(:,i)') = obj.n_t;
                obj.A(idx,obj.y_i(i)) = -8766*(1-obj.o_i_p(i));
                obj.b(idx) = 0;
                idx = idx + 1;
            end
            
            % CONSTRAINT SET 5
            % ainimum generation constraint (54 rows)
            % x_it >= x_i_min for i={1,2,3,4,5,6,7,8,9}, t={1,2,3,4,5,6}
            for i = 1:9
                for t = 1:6
                    obj.A(idx,obj.x_ti(t,i)) = -1;
                    obj.b(idx) = -obj.x_i_min(i);
                    idx = idx + 1;
                end
            end 
            
            % CONSTRAINT SET 6
            % lower bound on DSM program (3 rows)
            % z_k >= 0 for k={1,2,3}
            for k = 1:3
                obj.A(idx,obj.z_k(k)) = -1;
                obj.b(idx) = 0;
                idx = idx + 1;
            end
            
            % CONSTRAINT SET 7
            % upper bound on DSM program (3 rows)
            % z_k <= 1 for k={1,2,3}
            for k = 1:3
                obj.A(idx,obj.z_k(k)) = 1;
                obj.b(idx) = 1;
                idx = idx + 1;
            end
            
            % CONSTRAINT SET 8
            % new generation bounds (24 rows)
            % x_it <= (1-o_i_u)y_i for i={5,6,7,8}, t={1,2,3,4,5,6}
            for i = 5:8
                for t = 1:6
                    obj.A(idx,obj.x_ti(t,i)) = 1;
                    obj.A(idx,obj.y_i(i)) = -(1-obj.o_i_u(i));
                    obj.b(idx) = 0;
                    idx = idx + 1;
                end
            end 
            
            % CONSTRAINT SET 9
            % new generation bounds (4 rows)
            % y_i <= x_i_max for i={5,6,7,8}
            for i = 5:8
                obj.A(idx,obj.y_i(i)) = 1;
                obj.b(idx) = obj.x_i_max(i);
                idx = idx + 1;
            end
            
        end
        
        function gen_cost_obj_funs(obj)
            %{
                this function generates the cost objective functions (for
                both minimization and maximization), along with the vectors
                and matrices needed to obtain the variance
                Properties modified:
                    - cost_c_o
                    - cost_H
                    - cost_min_f
                    - cost_max_f
                    - cost_SIGMA
                :param:
                :return:
            %}
            
            % set up cost objective function matrices
            obj.cost_c_o = [obj.c_i_v; obj.c_i_c(5:8); obj.c_k_d];
            obj.cost_H = zeros(numel(obj.cost_c_o),size(obj.A,2));
            v_idx = 1;
            h_idx = 1;

            % variable cost component
            % sum_i=1^9{sum_t=1^6{c_i_v n_t x_it}}
            for i = 1:9
                obj.cost_H(v_idx,1+numel(obj.n_t)*(i-1):numel(obj.n_t)+ ... 
                    numel(obj.n_t)*(i-1)) = obj.n_t;
                v_idx = v_idx + 1;
                h_idx = h_idx + numel(obj.n_t);
            end

            % capital cost component
            % sum_i=5^8{c_i_vy_i10^3}
            for i = 5:8
                obj.cost_H(v_idx,h_idx) = 1000;
                v_idx = v_idx + 1;
                h_idx = h_idx + 1;
            end

            % demand side management cost component
            % sum_k=1^3{sum_t=1^6{c_k_d z_k s_kt}}
            for k = 1:3
                val = 0;
                for t = 1:6
                    val = val + obj.s_kt_max(k,t)*obj.n_t(t);
                end
                obj.cost_H(v_idx,h_idx) = val;
                v_idx = v_idx + 1;
                h_idx = h_idx + 1;
            end
            
            % set cost min/max vector for linprog
            obj.cost_min_f = obj.cost_c_o'*obj.cost_H;
            obj.cost_max_f = -obj.cost_c_o'*obj.cost_H;
            
            % generate variance covariance matrix for cost
            obj.cost_SIGMA = diag([obj.sigma_i_v; obj.sigma_i_c(5:8); ...
                obj.sigma_k_d]).^2;
            
        end
        
        function gen_emi_obj_funs(obj)
            %{
                this function generates the emissions objective functions 
                (for both minimization and maximization)
                Properties modified:
                    - emi_H
                    - emi_min_f
                    - emi_max_f
                :param:
                :return:
            %}
            
            % set up cost objective function matrices
            obj.emi_H = zeros(numel(obj.g_i),size(obj.A,2));
            v_idx = 1;
            
            % variable emissions component
            % sum_i=1^9{sum_t=1^6{g_i_ n_t x_it}}
            for i = 1:9
                obj.emi_H(v_idx,1+numel(obj.n_t)*(i-1):numel(obj.n_t)+ ... 
                    numel(obj.n_t)*(i-1)) = obj.n_t;
                v_idx = v_idx + 1;
            end
            
            % set cost min/max vector for linprog
            obj.emi_min_f = obj.g_i'*obj.emi_H;
            obj.emi_max_f = -obj.g_i'*obj.emi_H;
        end
        
    end
    
end