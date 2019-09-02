function [L,D,U,Aprime] = LUdecomp(A) 
    %LUdecomp performs the LU decomposition on the matrix A. It assumes A
    %is square and invertible. This function was written by Emma
    %Benjaminson and draws on information provided by: 
    %https://www.csun.edu/~panferov/math262/262_rref.pdf
    %MATLAB function rref()
    
    [m, ~] = size(A) ; % assume A is square so just need one dimension
    
    i = 1 ; % row index
    j = 1 ; % col index
    
    eps = 0.001 ; % some tolerance for small numbers that should be 0
    
    L = zeros(m,m) ; % matrix L is going to track Gaussian elimination process
    
    while i <= m && j <= m
        
        % look for max non-zero pivot element in current col
%         [M, k] = max(abs(A(i:m, j))) ; 
        [M, k] = max(A(i:m, j)) ; 
        k = k + i - 1 ; % because k will be an index with respect to A(i:m, j), not A(:,:)
        
        if M < eps
            A(i:m,j) = 0 ; 
            j = j + 1 ; 
        else 
            % swap current row with the row containing the max non-zero pivot
            % element
            A([i k], j:m) = A([k i], j:m) ; 
            
            % also perform the swap for L
            L([i k], :) = L([k i], :) ; 
            
            % save the division factor in L
            L(i,j) = 1/A(i,j) ; 
            
            % divide the i-th row by a_ij to make the pivot = 1 
            A(i,j:m) = A(i,j:m)./A(i,j) ; 
                      
            % now subtract multiples of the i-th row (which now begins with 1)
            % from all other rows. The multiple is the first non-zero value
            % of the destination row, because this will zero out that
            % element in the matrix
            for k = [1:i-1 i+1:m] % goes through every row except the i-th row
                L(k,i) = A(k,j) ; 
                L(k,k) = 1 ;
                A(k,j:m) = A(k,j:m) - A(k,j).*A(i,j:m) ; 
            end
        end
        
        i = i + 1 ; 
        j = j + 1 ; 
        
    end
    
    Aprime = A ; 
    
    % now we have A where A = DU
    D = Aprime * eye(m) ; 
    U = inv(D) * Aprime ; 
end