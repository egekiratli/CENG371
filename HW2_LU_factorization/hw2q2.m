sizes = 2:300;
numSizes = length(sizes);

% Arrays to store run times and relative errors
timeLUWP = zeros(1, numSizes);
timeLUPP = zeros(1, numSizes);
timeLUCP= zeros(1, numSizes);

relErrorLUWP = zeros(1, numSizes);
relErrorLUPP = zeros(1, numSizes);
relErrorLUCP = zeros(1, numSizes);


for i = 1:numSizes
    
    An = hilb(i);

    % LU factorization without pivoting
    tic;
    [L1, U1] = LUWP(An);
    timeLUWP(i) = toc;

    relErrorLUWP(i) = norm(An - L1 * U1) / norm(An);

    % LU factorization with partial pivoting
    tic;
    [P2, L2, U2] = LUPP(An);
    timeLUPP(i) = toc;

    relErrorLUPP(i) = norm(P2 * An - L2 * U2) / norm(P2 * An);

    % LU factorization with complete pivoting
    tic;
    [P3, L3, U3, Q3] = LUCP(An);
    timeLUCP(i) = toc;

    relErrorLUCP(i) = norm(P3 * An * Q3 - L3 * U3) / norm(P3 * An * Q3);
end

% Plotting
figure;

subplot(2, 1, 1);
semilogy(sizes, timeLUWP, 'b-', 'LineWidth', 2, 'DisplayName', 'Without Pivoting');
hold on;
semilogy(sizes, timeLUPP, 'g-', 'LineWidth', 2, 'DisplayName', 'Partial Pivoting');
hold on;
semilogy(sizes, timeLUCP, 'r-', 'LineWidth', 2, 'DisplayName', 'Complete Pivoting');
xlabel('Matrix Size');
ylabel('Run Time');
legend('Location','southeast');
grid on;

subplot(2, 1, 2);
semilogy(sizes, relErrorLUWP, 'b-', 'LineWidth', 2, 'DisplayName', 'Without Pivoting');
hold on;
semilogy(sizes, relErrorLUPP, 'g-', 'LineWidth', 2, 'DisplayName', 'Partial Pivoting');
hold on;
semilogy(sizes, relErrorLUCP, 'r-', 'LineWidth', 2, 'DisplayName', 'Complete Pivoting');
xlabel('Matrix Size');
ylabel('Relative Error');
legend('Location','east');
grid on;

hold off;

condLUWP = cond(An);
condLUPP = cond(P2 * An);
condLUCP = cond(P3 * An * Q3);


%Functions Appendix

function [L,U,P,Q] = LUCP(A)
    n = length(A)
    L = eye(n);
    P = eye(n);
    Q = eye(n);
    U = A;

    tol = 1e-12; %define tolerance to prevent rounding errors
    
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

function [L,U,P] = LUPP(A)
    n = length(A)
    L = eye(n);
    P = eye(n);
    U = A;
    
    tol = 1e-12; %define tolerance to prevent rounding errors

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


function [L,U] = LUWP(A)
    n = length(A);
    L = eye(n);
    U = A;
    
    tol = 1e-12; %define tolerance to prevent rounding errors

    for i = 1:n-1
        for j = i+1:n
            if abs(U(j,i)) > tol
                L(j,i) = U(j,i)/U(i,i); %factorization value
                U(j,:) = U(j,:)-L(j,i)*U(i,:); 
            end
        end
     end
end

