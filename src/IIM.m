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
            gamma_i = obj.compute_row_sum(obj.A,i);
        end
        
        function delta_j = get_delta_j(obj, j)
            %{
                get delta_j from matrix A at column index j
                :param j: column index for delta computation
                :return delta_j: delta from column index j
            %}
            delta_j = obj.compute_col_sum(obj.A,j);
        end
        
        function gamma_bar_i = get_gamma_bar_i(obj, i)
            %{
                get gamma_i from matrix S at row index i
                :param i: row index for gamma computation
                :return gamma_bar_i: gamma from row index i
            %}
            gamma_bar_i = obj.compute_row_sum(obj.S,i);
        end
        
        function delta_bar_j = get_delta_bar_j(obj, j)
            %{
                get delta_j from matrix S at column index j
                :param j: column index for delta computation
                :return delta_bar_j: delta from column index j
            %}
            delta_bar_j = obj.compute_col_sum(obj.S,j);
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
        
        function row_sum = compute_row_sum(~, M, i)
            %{
                compute the sum of row i of matrix M
                :param M: matrix from which row sum will be computed
                :param i: index for row sum
                :return row_sum: row sum from matrix M at row i
            %}
            row_sum = sum(M(i,:))/(length(M(i,:)) - 1);
        end
        
        function col_sum = compute_col_sum(~, M, j)
            %{
                compute the sum of column j of matrix M
                :param M: matrix from which column sum will be computed
                :param j: index for column sum
                :return col_sum: column sum from matrix M at column j
            %}
            col_sum = sum(M(:,j))/(length(M(:,j)) - 1);
        end
        
        function x_k = compute_damage_propagation(obj, f, k)   
            i = 1;
            x_k = zeros(length(obj.A),k+1);
            x_k = propagate_damage(obj, x_k, i+1, f, k);
        end
        
        function x_k = propagate_damage(obj,x_k,i,f,k)
            x_k(:,i) = obj.A*x_k(:,i-1) + f;
            if i < k
                x_k = propagate_damage(obj, x_k, i+1, f, k);
            end
        end
        
    end
    
end

