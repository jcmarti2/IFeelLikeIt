
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
        
        options
        
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
        
        var_quad_H
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
            
            % construct matrix for variance minimization
            obj.gen_var_obj_funs()
            
            obj.options = optimoptions('linprog','Algorithm', ...
                'dual-simplex','Display', 'off');
        end
           
        function [x, cost, emissions, variance, slack, binding, ...
                shadow_cost, shadow_emi, shadow_var, exitflag, output, ...
                lambda] = minimize_cost(obj)
            %{
                this function returns the results from the cost
                minimization optimization problem
                :param:
                :return x: vector of decision variables
                :return cost: optimized cost (at minimum cost [$]
                :return emissions: emissions (at minimum cost) [tons]
                :return variance: cost variance
                :return slack: vector of indices of slack variables
                :return binding: vector of indices of binding variables
                :return shadow_cost: vector of shadow prices for binding
                                variables
                :return shadow_emi: vector of shadow emissions for binding
                                    variables
                :return shadow_var: vector of shadow variance for binding
                :return exitflag: exitflag from linprog
                :return output: output from linprog
                :return lambda: lambda from linprog
            %}
                [x,cost,exitflag,output,lambda] = linprog( ...
                    obj.cost_min_f,obj.A,obj.b,[],[],[],[],[],obj.options);
                emissions = (obj.g_i'*obj.emi_H*x)*(10^-3);     
                variance = x'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*x;
                
                slack = find(obj.b - obj.A*abs(x) > 10^-10);
                binding = find(~(obj.b - obj.A*abs(x) > 10^-10));
                shadow_cost = zeros(numel(binding),1);
                shadow_emi = zeros(numel(binding),1);
                shadow_var = zeros(numel(binding),1);
                for binding_idx = 1:numel(binding)
                    new_b = obj.b;
                    if binding(binding_idx) >= 124 && ...
                            binding(binding_idx) <= 129;
                        agg_val = 0.1;
                    else
                        agg_val = 1;
                    end
                    new_b(binding(binding_idx)) = ... 
                        new_b(binding(binding_idx)) + agg_val;
                    [new_x, new_cost] = linprog(obj.cost_min_f, obj.A, ... 
                        new_b,[],[],[],[],[],obj.options);
                    shadow_cost(binding_idx) = new_cost - cost;
                    shadow_emi(binding_idx) = (obj.g_i'*obj.emi_H*new_x)*(10^-3) ...
                        - emissions;
                    shadow_var(binding_idx) = new_x'*obj.cost_H'* ... 
                        obj.cost_SIGMA*obj.cost_H*new_x;
                end
        end
        
        function [x, cost, emissions, variance, slack, binding, ...
                exitflag, output, lambda] = maximize_cost(obj)
            %{
                this function returns the results from the cost
                maximization optimization problem
                :param:
                :return x: vector of decision variables
                :return cost: optimized cost (at maximum cost) [$]
                :return emissions: emissions (at maximum cost) [tons]
                :return variance: cost variance
                :return slack: vector of indices of slack variables
                :return binding: vector of indices of binding variables
                :return exitflag: exitflag from linprog
                :return output: output from linprog
                :return lambda: lambda from linprog
            %}
                [x,cost,exitflag,output,lambda] = linprog( ...
                    obj.cost_max_f,obj.A,obj.b,[],[],[],[],[],obj.options);
                emissions = (obj.g_i'*obj.emi_H*x)*(10^-3);
                variance = x'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*x; 
                slack = find(~(obj.A*x == obj.b));
                binding = find(obj.A*x == obj.b);    
        end
        
        function [x, cost, emissions, variance, slack, binding, ...
                shadow_cost, shadow_emi, shadow_var, exitflag, output, ...
                lambda] = minimize_emissions(obj)
            %{
                this function returns the results from the emissions
                minimization optimization problem
                :param:
                :return x: vector of decision variables
                :return cost: cost (at minimum emissions) [$]
                :return emissions: optimized emissions (at minimum
                        emissions) [tons]
                :return variance: cost variance
                :return slack: vector of indices of slack variables
                :return binding: vector of indices of binding variables
                :return shadow_cost: vector of shadow prices for binding var
                :return shadow_emi: vector of shadow emissions for binding
                                    variables
                :return shadow_var: vector of shadow variance for binding
                :return exitflag: exitflag from linprog
                :return output: output from linprog
                :return lambda: lambda from linprog
            %}

                [x,emissions,exitflag,output,lambda] = linprog( ...
                    obj.emi_min_f, obj.A, obj.b,[],[],[],[],[],obj.options);
                cost = obj.cost_c_o'*obj.cost_H*x;
                emissions = emissions*(10^-3);
                variance = x'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*x;
                
                slack = find(obj.b - obj.A*abs(x) > 10^-10);
                binding = find(~(obj.b - obj.A*abs(x) > 10^-10));
                shadow_cost = zeros(numel(binding),1);
                shadow_emi = zeros(numel(binding),1);
                shadow_var = zeros(numel(binding,1));
                for binding_idx = 1:numel(binding)
                    new_b = obj.b;
                    if binding(binding_idx) >= 124 && ...
                            binding(binding_idx) <= 129;
                        agg_val = 0.1;
                    else
                        agg_val = 1;
                    end
                    new_b(binding(binding_idx)) = ... 
                        new_b(binding(binding_idx)) + agg_val;
                    [new_x, new_emi] = linprog(obj.emi_min_f, obj.A, ... 
                        new_b,[],[],[],[],[],obj.options);  
                    new_emi = new_emi*(10^-3);
                    shadow_emi(binding_idx) = new_emi - emissions;
                    shadow_cost(binding_idx) = obj.cost_c_o'* ... 
                        obj.cost_H*new_x - cost;
                    shadow_var(binding_idx) = new_x'*obj.cost_H'* ... 
                        obj.cost_SIGMA*obj.cost_H*new_x;
                end
        end
 
        function [x, cost, emissions, variance, slack, binding, ...
                exitflag, output, lambda] = maximize_emissions(obj)
            %{
                this function returns the results from the emissions
                maximization optimization problem
                :param:
                :return x: vector of decision variables
                :return cost: cost (at maximum emissions) [$]
                :return emissions: optimized emissions (at maximujm
                                    emissions) [tons]
                :return variance: cost variance
                :return slack: vector of indices of slack variables
                :return binding: vector of indices of binding variables
                :return exitflag: exitflag from linprog
                :return output: output from linprog
                :return lambda: lambda from linprog
            %}
                [x,emissions,exitflag,output,lambda] = linprog( ...
                    obj.emi_max_f, obj.A, obj.b,[],[],[],[],[],obj.options);
                emissions = emissions*(10^-3);
                cost = obj.cost_c_o'*obj.cost_H*x;
                variance = x'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*x;
                slack = find(~(obj.A*x == obj.b));
                binding = find(obj.A*x == obj.b);
        end
        
        function [costs, emissions, decisions] = minimize_cost_and_emissions(obj)
            %{
                this function computes the pareto optimality curve for
                costs and emissions
                :param:
                :return costs: vector of pareto costs [$]
                :return emissions: vector of pareto emissions [tons]
                :return decisions: matrix were each column is a decision
                vector that corresponds to each index of costs and
                emissions
            %}
            
            % intialize empty data structures
            num_steps = 1000;
            lambdas = linspace(0,1,num_steps);
            costs = zeros(1,num_steps);
            emissions = zeros(1,num_steps);
            decisions = zeros(size(obj.cost_min_f,2),num_steps);
            
            % obtain minimum emissions and minimum cost (for normalizing
            % pareto)
            [~,min_cost] = linprog(obj.cost_min_f, obj.A, obj.b,[],[], ...
                [],[],[],obj.options);
            [~,min_emi] = linprog(obj.emi_min_f, obj.A, obj.b,[],[], ...
                [],[],[],obj.options);    
            
            % compute pareto curve
            for idx = 1:num_steps
                lambda = lambdas(idx);
                % get new objective f for linprog and X as a result
                f = (lambda/min_cost).*(obj.cost_c_o'*obj.cost_H) + ...
                    ((1-lambda)/min_emi).*(obj.g_i'*obj.emi_H);
                X = linprog(f, obj.A, obj.b,[],[],[],[],[],obj.options);
                % compute new costs and emissions
                costs(idx) = obj.cost_c_o'*obj.cost_H*X;
                emissions(idx) = (obj.g_i'*obj.emi_H*X)*(10^-3);
                decisions(:,idx) = X';
            end
        end
        
        function [costs, caps] = minimize_cost_with_emissions_cap(obj)
            %{
                this function sets emissions caps and minimizes cost
                :param:
                :return costs:
                :return caps:
            %}
                    
            % minimize emissions
            [~, min_emi] = linprog(obj.emi_min_f, obj.A, obj.b,[],[], ... 
                [],[],[],obj.options);
            
            num_caps = 1000;
            
            % new A and b
            new_A = [obj.A; obj.emi_min_f];
            caps = linspace(0,1,num_caps);
            costs = zeros(numel(caps));
            for i = 1:num_caps
                new_b_row = min_emi + min_emi*caps(i);
                new_b = [obj.b; (new_b_row)];
                [~,cost] = linprog(obj.cost_min_f,new_A,new_b,[],[],[], ...
                    [],[],obj.options);
                caps(i) = new_b_row;
                costs(i) = cost;
            end
            
        end
            
        function [rates, costs, emissions] = tax_emissions(obj)
            %{
                this function computes the costs and emissions of the power
                system for various tax rates
                :param:
                :return rates: vector of rates used [%]
                :return costs: vector of costs for a given rate [$]
                :return emissions: vector of emissions for a given rate
                                   [tons]
            %}
            
            num_rates = 1000;
            max_rate = 50;
            rates = linspace(0,max_rate,num_rates);
            costs = zeros(num_rates);
            emissions = zeros(num_rates);
            
            % compute cost and emissions for each rate
            for i = 1:num_rates
                rate = rates(i);
                tax_f = obj.cost_min_f + obj.emi_min_f*(rate/100);
                [x,cost] = linprog(tax_f,obj.A,obj.b,[],[],[],[],[], ... 
                    obj.options);
                    costs(i) = cost;
                    emi = obj.g_i'*obj.emi_H*x;
                    emissions(i) = emi;   
            end
        end
        
        function [x, cost, emissions, variance, slack, binding, ...
                shadow_cost, shadow_emi, shadow_var, exitflag, output, ...
                lambda] = minimize_variance(obj)
            %{
                this function returns the results from the variance
                minimization optimization problem
                :param:
                :return x: vector of decision variables
                :return cost: cost (at minimum emissions) [$]
                :return emissions: optimized emissions (at minimum
                        emissions) [tons]
                :return variance: cost variance
                :return slack: vector of indices of slack variables
                :return binding: vector of indices of binding variables
                :return shadow_cost: vector of shadow prices for binding var
                :return shadow_emi: vector of shadow emissions for binding
                                    variables
                :return shadow_var: vector of shadow variance for binding
                :return exitflag: exitflag from linprog
                :return output: output from linprog
                :return lambda: lambda from linprog
            %}
            
            [x, variance, exitflag, output, lambda] = quadprog(...
                obj.var_quad_H,[],obj.A,obj.b);
            cost = obj.cost_c_o'*obj.cost_H*x;
            emissions = (obj.g_i'*obj.emi_H*x)*(10^-3); 
            
            slack = find(obj.b - obj.A*abs(x) > 10^-10);
            binding = find(~(obj.b - obj.A*abs(x) > 10^-10));
            shadow_cost = zeros(numel(binding),1);
            shadow_emi = zeros(numel(binding),1);
            shadow_var = zeros(numel(binding),1);
            for binding_idx = 1:numel(binding)
                new_b = obj.b;
                if binding(binding_idx) >= 124 && ...
                     binding(binding_idx) <= 129;
                        agg_val = 0.1;
                else
                    agg_val = 1;
                end
                new_b(binding(binding_idx)) = ... 
                    new_b(binding(binding_idx)) + agg_val;
                [new_x, new_var] = quadprog(obj.var_quad_H,[],obj.A,new_b);
                shadow_var(binding_idx) = new_var - variance;
                shadow_cost(binding_idx) = obj.cost_c_o'* ... 
                    obj.cost_H*new_x - cost;
                shadow_emi(binding_idx) = (obj.g_i'*obj.emi_H*new_x)*(10^-3) ...
                    - emissions;
            end
        end
         
        function [costs, variances, decisions] = minimize_cost_and_variance(obj)
            %{
                this function computes the pareto optimality curve for
                costs and variance
                :param:
                :return costs: vector of pareto costs [$]
                :return variances: vector of pareto variances [$^2]
                :return decisions: matrix were each column is a decision
                vector that corresponds to each index of costs and
                variance
            %}
            
            % intialize empty data structures
            num_steps = 1000;
            lambdas = linspace(0,1,num_steps);
            costs = zeros(1,num_steps);
            variances = zeros(1,num_steps);
            decisions = zeros(size(obj.cost_min_f,2),num_steps);
            
            % obtain minimum emissions and minimum cost (for normalizing
            % pareto)
            [~,min_cost] = linprog(obj.cost_min_f, obj.A, obj.b,[],[], ...
                [],[],[],obj.options);
            [~,min_var] = quadprog(obj.var_quad_H,[],obj.A,obj.b);    
            
            % compute pareto curve
            for idx = 1:num_steps
                lambda = lambdas(idx);
                % get new objective f for linprog and X as a result
                f = (lambda/min_cost).*(obj.cost_c_o'*obj.cost_H);
                H = ((1-lambda)/min_var).*obj.var_quad_H;
                X = quadprog(H,f,obj.A,obj.b); 
                % compute new costs and emissions
                costs(idx) = obj.cost_c_o'*obj.cost_H*X;
                variances(idx) = X'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*X;
                decisions(:,idx) = X';
            end
        end
        
        function [costs, emissions, variances, decisions] = ...
                minimize_cost_emissions_and_variance(obj)
            
            num_steps = 40;
            lambdas = linspace(0,1,num_steps);
            costs = -1*ones(1,num_steps^2);
            emissions = -1*ones(1,num_steps^2);
            variances = -1*ones(1,num_steps^2);
            decisions = -1*ones(size(obj.cost_min_f,2),num_steps^2);
            
            opt = optimoptions('quadprog','Algorithm', ... 
                'trust-region-reflective','Display', 'off');
            
            % obtain minimum emissions and minimum cost (for normalizing
            % pareto)
            [~,min_cost] = linprog(obj.cost_min_f, obj.A, obj.b,[],[], ...
                [],[],[],obj.options);
            [~,min_emi] = linprog(obj.emi_min_f, obj.A, obj.b,[],[], ...
                [],[],[],obj.options); 
            [~,min_var] = quadprog(obj.var_quad_H,[],obj.A,obj.b); 
            
            for idx_1 = 1:num_steps
                lambda = lambdas(idx_1);
                gammas = linspace(0,1-lambda,num_steps);
                for idx_2 = 1:num_steps
                    gamma = gammas(idx_2);
                    f = (gamma/min_cost).*(obj.cost_min_f) + ...
                        ((1-lambda-gamma)/min_emi).*(obj.emi_min_f);
                    H = (lambda/min_var).*obj.var_quad_H;
                    X = quadprog(H,f,obj.A,obj.b);
                    % compute new costs and emissions
                    costs((idx_1-1)*num_steps + idx_2) = obj.cost_c_o'*obj.cost_H*X;
                    emissions((idx_1-1)*num_steps + idx_2) = (obj.g_i'*obj.emi_H*X)*(10^-3);
                    variances((idx_1-1)*num_steps + idx_2) = X'*obj.cost_H'*obj.cost_SIGMA*obj.cost_H*X;
                    decisions(:,(idx_1-1)*num_steps + idx_2) = X';
                end
            end
        end
        
        function cost = get_cost(obj, x)
            %{
                this function returns the cost given the decision variables
                :param x: vector of decision variables
                :return cost: total cost [$]
            %}
                cost = obj.cost_c_o'*obj.cost_H*x;
        end
        
        function emissions = get_emissions(obj, x)
            %{
                this function returns the emissions given the decision
                variables
                :param X: vector of decision variables
                :return emissions: total emissions [tons]
            %}
            emissions = (obj.g_i'*obj.emi_H*x)*(10^-3);
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
        
        function gen_var_obj_funs(obj)
            %{
                this function generates H for variance minimization
                quadprog
                Properties modified:
                    - var_quad_H
                :param:
                :return:
            %}
            
            % set H for quadprog
            obj.var_quad_H = 2*((obj.cost_H'*obj.cost_SIGMA*obj.cost_H));
        end
        
    end
    
end