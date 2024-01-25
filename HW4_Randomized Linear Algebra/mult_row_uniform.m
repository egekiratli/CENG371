function [D] = mult_row_uniform(A,B,c,~)
    
    [m,~] = size(A); %this n is not used (matlab warning)
    [n,p] = size(B);

    C = zeros(m,c);
    R = zeros(c,p);

    for t = 1:c
        i = randsample(n,1); %choose random index
        pk = 1/n; %since probability is uniform
        C(:,t) = A(:,i)/sqrt(c*pk);
        R(t,:) = B(i,:)/sqrt(c*pk);
    end

    D = C*R
end