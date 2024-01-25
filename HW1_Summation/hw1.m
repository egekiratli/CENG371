%Name: Ege KIRATLI
%ID: 2555258

%% Question 1
clear;
 
n = linspace(1,1000,1000);

f = @(n) n.*((n+1)./n-1)-1; %create the function

g = @(n) f(n)/eps; %create g function using machine epsilon


[roots1, count1] = findRoots(n,g)


plot(n,g(n),'.', MarkerSize=.2);
grid on
xlabel('n');
ylabel('g(n) = f/\epsilon');
text(50,-300,{'n = 10^3','(sum(g(n) == 0) = 10'}, FontSize=6)


%% Question 2

n2 = 10e6;
nums = @(n) 1 + (10e6 + 1 - n).*10e-8

% Naive Summation

%   in single precision:

startNaiveSingle = tic;

sumNaiveSingle = single(0);

for i = 1:n2
    sumNaiveSingle = sumNaiveSingle + nums(i); %add nums iteratively
end

telapsedNaiveSingle = toc(startNaiveSingle);

%   in double precision

startNaiveDouble = tic;

sumNaiveDouble = double(0);

for i = 1:n2
    sumNaiveDouble = sumNaiveDouble + nums(i); %add nums iteratively
end

telapsedNaiveDouble = toc(startNaiveDouble);


% Compensated (Kahan) Summation

%   in single precision

startCompensatedSingle = tic;

sumCompensatedSingle = single(0); %single precision
errCompensatedSingle = single(0); %var to store error

for i = 1:n2
    correctedValCompensatedSingle = single(nums(i)) - errCompensatedSingle; %find corrected value of nums
    intermediateSumCompensatedSingle = sumCompensatedSingle + correctedValCompensatedSingle; %calculate intermediate sum
    errCompensatedSingle = (intermediateSumCompensatedSingle - sumCompensatedSingle) - correctedValCompensatedSingle; %compute the error
    sumCompensatedSingle = intermediateSumCompensatedSingle; 
end

telapsedCompensatedSingle = toc(startCompensatedSingle);

%   in double precision

startCompensatedDouble = tic;

sumCompensatedDouble = double(0); % double precision
errCompensatedDouble = double(0); % var to store error

for i = 1:n2
    correctedValCompensatedDouble = nums(i) - errCompensatedDouble; % find corrected value of nums
    intermediateSumCompensatedDouble = sumCompensatedDouble + correctedValCompensatedDouble; % calculate intermediate sum
    errCompensatedDouble = (intermediateSumCompensatedDouble - sumCompensatedDouble) - correctedValCompensatedDouble; % compute the error
    sumCompensatedDouble = intermediateSumCompensatedDouble;
end
telapsedCompensatedDouble = toc(startCompensatedDouble);

% Pairwise Summation

%   in single precision

startPairwiseSingle = tic;

sumPairwiseSingle = single(0);

for i = 1:2:n2
    sumPairwiseSingle = sumPairwiseSingle + single(nums(i-1)) + single(nums(i));
end

telapsedPairwiseSingle = toc(startPairwiseSingle);

%   in double precision

startPairwiseDouble = tic;

sumPairwiseDouble = double(0);

for i = 2:2:n2
    sumPairwiseDouble = sumPairwiseDouble + double(nums(i-1)) + double(nums(i));
end

telapsedPairwiseDouble = toc(startPairwiseDouble);


%% Appendix for functions

function [roots,count] = findRoots(n, g)
    roots = [];
    count = 0;
    for i = 1:length(n)
        if g(n(i)) == 0
            count = count + 1;
            roots = [roots, n(i)];
        end
    end
end