%Name: Ege KIRATLI
%ID: 2555258

%% GAUSSIAN ELIMINATION W/ COMPLETE PIVOTING

A = [5,3,2,9 ; 1,2,7,3 ; 8,5,3,1 ; 2,1,3,0]
[L,U,P,Q] = LUCP(A)

function [L,U,P,Q] = LUCP(A)
    n = length(A)
    L = eye(n);
    P = eye(n);
    Q = eye(n);
    U = A;

    tol = 1e-6; %define tolerance to prevent rounding errors
    
    for i = 1:n-1
        [~,pivotCol] = max(abs(U(i,i:n))); 
        [~,pivotRow] = max(abs(U(i:n,i))); 
    
        pivotCol = pivotCol + i - 1;
        pivotRow = pivotRow + i - 1; %adjust the index to the relevant submatrix
        
        U([i,pivotRow],:) = U([pivotRow,i],:); %swap the rows
        P([i,pivotRow],:) = P([pivotRow,i],:); 

        U(:,[i,pivotCol]) = U(:,[pivotCol,i]); %swap the cols
        Q(:,[i,pivotCol]) = Q(:,[pivotCol,i]); 
    
        for j = i+1:n
            if abs(U(j,i)) > tol
                L(j,i) = U(j,i)/U(i,i); 
                U(j,:) = U(j,:)-L(j,i)*U(i,:); 
            end
        end
    end
end