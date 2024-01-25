%Name: Ege KIRATLI
%ID: 2555258

%% GAUSSIAN ELIMINATION W/ PARTIAL PIVOTING

A = [5,3,2,2 ; 1,2,7,6 ; 8,5,3,7 ; 2,1,3,0]
[L,U,P] = LUPP(A)

function [L,U,P] = LUPP(A)
    n = length(A)
    L = eye(n);
    P = eye(n);
    U = A;
    
    tol = 1e-6; %define tolerance to prevent rounding errors

    for i = 1:n-1
        [~,pivotRow] = max(abs(U(i:n,i))); 
        
        pivotRow = pivotRow + i - 1; %adjust the index to the relevant submatrix
        
        U([i,pivotRow],:) = U([pivotRow,i],:); %swap the rows
        P([i,pivotRow],:) = P([pivotRow,i],:); 

        for j = i+1:n
            if abs(U(j,i)) > tol
                L(j,i) = U(j,i)/U(i,i); %factorization value
                U(j,:) = U(j,:)-L(j,i)*U(i,:); 
            end
        end
    end
end