function [D] = mult_proj_Gauss_orth(A,B,c,~)

    [m,~] = size(A);
    [n,p] = size(B);


    P = randn(n,c); %random projection matrix
    
    %Gram-Schmidt process
    for i = 1:c
        for j = 1:i-1
            P(:,i) = P(:,i) - (P(:,j)'*P(:,i)*P(:,j));
        end
        P(:,i) = P(:,i)/norm(P(:,i));
    end

    C = A*P;
    R = P'*B;


    D = C*R;

end