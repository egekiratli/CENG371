
function [D] = mult_row_nonuni(A,B,c,~)

    [m,~] = size(A); %this n is not used (matlab warning)
    [n,p] = size(B);

    C = zeros(m,c);
    R = zeros(c,p);
    
    normA = sum(sqrt(sum(A.^2,2))) %norm of cols of A
    normB = sum(sqrt(sum(B.^2,1))) %norm of rows of B

    for t = 1:c
        i = randsample(n,1); %choose a random index
        pi = (norm(A(:,i))*norm(B(i,:)))/(normA*normB);
        C(:,t) = A(:,i)/sqrt(c*pi);
        R(t,:) = B(i,:)/sqrt(c*pi);
    end

    D = C*R
end