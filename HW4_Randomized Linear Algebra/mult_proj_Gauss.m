function [D] = mult_proj_Gauss(A,B,c,~)
    

    [m,~] = size(A);
    [n,p] = size(B);
    

    P = randn(n,c); %random projection matrix
    for i = 1:c
        P(:,i) = P(:,i)/norm(P(:,i))
    end

    C = A*P;
    R = P'*B;

    D = C*R

end