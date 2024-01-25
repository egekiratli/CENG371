%Name: Ege Kıratlı
%ID: 2555258

%% Question 1

img = imread("Santa_Claus.jpg");
grayImg = double(rgb2gray(img)); %convert rbg to grayscale + uint8 to double

[U, S, V] = svd(grayImg); 

r = 1200; %it's the min(m,n) value

k1 = 1 ; S_k1 = S;
kr2 = r/2; S_kr2 = S;
kr10 = r-10; S_kr10 = S;

S_k1(k1+1:end,k1+1:end) = 0;
S_kr2(kr2+1:end,kr2+1:end) = 0;
S_kr10(kr10+1:end,kr10+1:end) = 0;


%low rank approximations for all r
S_kVal = zeros(r,1);
Ierr = zeros(r,1);

for k = 1:r
    S_k = S;
    S_k(k+1:end,k+1:end) = 0;

    I_k = U*S_k*V;
    Ierr(k) = norm(grayImg-I_k, 'fro');
    S_kSingular = diag(S_k);
    S_kVal(k) = sum(S_kSingular.^2)/sum(diag(S).^2);
end

Ierr(1)
Ierr(r/2)
Ierr(r-10)



subplot(2,1,1);
loglog(1:r, S_kVal);
title('S(k) vs k');
xlabel('k');
ylabel('S(k)');

subplot(2,1,2);
loglog(1:r, Ierr);
title('Frobenius Norm Error vs. k');
xlabel('k');
ylabel('Error');



%% Question 2

%function handles
e = @(n) randn(1,n)*.01;
t = @(n) linspace(1,n,n)/(n+1);

N = @(t) 0.3 + 2*(t) - 1.2*(t).^2 + 0.5*(t).^3;
N_obs = @(t,n) 0.3 + 2*(t+e(n)) - 1.2*(t+e(n)).^2 + 0.5*(t+e(n)).^3;

%observation sets
O_1 = [t(5);N_obs(t(5),5)];
O_2 = [t(10);N_obs(t(10),10)];
O_3 = [t(100);N_obs(t(100),100)];

%construct the A matrix where each column corresponds to a power of t
%we know N is of the form a + bt + ct^2 + dt^3
A = @(n) [ones(n,1), t(n)', t(n)'.^2, t(n)'.^3];

%observed data matrix
%we need nx1 (b-transpose) matrix for solution
b = @(n) N_obs(t(n),n)';

ck = [.3, 2, -1.2, .5];

%using backslash operator
c5 = A(5)\b(5);
c10 = A(10)\b(10);
c100 = A(100)\b(100);

%error for backslash operator
c5Err = norm(ck-c5);
c10Err = norm(ck-c10);
c100Err = norm(ck-c100);


%using least squares minimum function
c5Least = lsqminnorm(A(5),b(5));
c10Least = lsqminnorm(A(10),b(10));
c100Least = lsqminnorm(A(100),b(100));

%error for least squares function
c5LeastErr = norm(ck-c5Least);
c10LeastErr = norm(ck-c10Least);
c100LeastErr = norm(ck-c100Least);


%Note: I've read that MATLAB decides on the method based on the 
%characteristics of the matrix when using backslash operator, 
%my guess is that it uses least squares solution 
%since it's an overdetermined case. 
%However, I wanted to include the lsqminnorm function too to show my
%intent.

