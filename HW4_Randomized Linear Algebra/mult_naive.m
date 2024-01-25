function [C] =  mult_naive(A,B)
    
    [m,~] = size(A);
    [n,p] = size(B);

    AB = zeros(m,p);
    

    for i = 1:m
        for j = 1:p
            for k = 1:n
                AB(i,k) = AB(i,k) + A(i,j)*B(j,k);
            end
        end
    end
    
    C = AB
end

