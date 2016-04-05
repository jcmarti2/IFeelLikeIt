classdef IIM
    % IIM 
    %{
        this class is an implementation of the Inoperability Input-Output
        model framework for dependency analysis
        
        authors:
            Juan Carlos Martinez
            Paul Gharzouzi
            Jimmy Chang
    %}
    %{   
        PUBLIC:
            get_A
            get_S
            get_gamma_i
            get_delta_j
            get_gamma_bar_i
            get_delta_bar_j
            get_x
    %}
    
    properties
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
           obj.A = obj.load_A(loc_A);
           obj.S = obj.compute_S();
        end
        
        function A = get_A(obj)
            %{
                get property A of IIM class instance
                :param:
                :return A: property A of IIM class instance
            %}
            A = obj.A;
        end
        
        function S = get_S(obj)
            %{
                get property S of IIM class instance
                :param:
                :return S: property S of IIM class instance
            %}
            S = obj.S;
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
        
        function x = get_x(obj, f)
            %{
                get x given f
                :param f: external disturbance vector
                :return x: degradation vector
            %}
            x = obj.compute_x(f);
        end
        
    end
        
       
    methods(Access = private)
        
        function A = load_A(obj, loc_A)
            %{
                load matrix A from loc_A .csv file
                :param loc_A: address of .csv file holding A
                :return A: matrix A
            %}
            A = loc_A;
        end
        
        function S = compute_S(obj)
            %{
                compute matrix S using property A of the IIM class
                :param:
                :return S: matrix S
            %}
        end
        
        function x = compute_x(obj, f)
            %{
                get x given f
                :param f: external disturbance vector
                :return x: degradation vector
            %}
        end
        
        function row_sum = compute_row_sum(obj, M, i)
            %{
                compute the sum of row i of matrix M
                :param M: matrix from which sum will be computed
                :param i: index for row sum
                :return row_sum: row sum from matrix M at row i
            %}
        end
        
        function col_sum = compute_col_sum(obj, M, j)
            %{
                compute the sum of column j of matrix M
                :param M: matrix from which sum will be computed
                :param j: index for column sum
                :return col_sum: column sum from matrix M at column j
            %}
        end
        
    end
    
end

