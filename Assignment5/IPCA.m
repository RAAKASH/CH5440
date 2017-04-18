%% Iterative PCA
clc;
clear;
close all;
load('steamdata.mat');
clear std
%% Initialziation
m = 28 ;  % No of variales (ie) F1 to F5
SigmaE = eye(m); % Initial guess of Error variances (Assuming Covariances  = 0 )
iter = 0;
beta = 0.4; % For convergence (relaxation parameter)
No = 28; % No of constraints
%% IPCA code
while(iter<100)
L = SigmaE.^0.5;     % L matrix (check notes)
D = inv(L)*(Fmeas); % Transforming data
[U,S,V] = svd(D,'econ');    % Obtaining U matrix
c = (U(:,1:end)'/L); % Retransforming it for actual data set
 c = c((1):end,:);
r = c*Fmeas ;            % residual matrix
XE = r*r'/10^3 ;     % Covariance matix for errors
% X_e = [XE(1,1:5)';XE(2,2:5)';XE(3,3:5)';XE(4,4:5)';XE(5,5:5)']; % Vec of ri*ri' (only the upper triangular matrix)
% X_e = [XE(1,1:3)';XE(2,2:3)';XE(3,3:3)'];
for i=1:No
    if(i==1)
     X_e = XE(i,i:No)';
    else
     X_e = [X_e;XE(i,i:No)'];

    end
end
Xe =diag(XE);
[m,n] = size(c);
A = zeros(m*m,n*n); % For neudecker product ,Variable Initialization
for i=1:m
    for j=1:n
    A(((i-1)*m+1):i*m,((j-1)*n+1):n*j) = c(i,j)*c; % The Neudecker product of CC
    end
end
% A_dash = [A(1:5,:);A(7:10,:);A(13:15,:);A(19:20,:);A(25,:)];% Vec of ri*ri'
% A_dash = [A(1:3,:);A(5:6,:);A(6:6,:)];% Vec of ri*ri'
for i =1:No
    if(i==1)
A_dash = A((i-1)*No+i:No*i,:);% Vec of ri*ri'
    else
A_dash = [A_dash;A((i-1)*No+i:No*i,:)];% Vec of ri*ri'
    end
end
% N =[A_dash(:,1),A_dash(:,7),A_dash(:,13),A_dash(:,19),A_dash(:,25)];
% Sigma_E = (diag([(N'*N)\(N'*X_e)]));
N = [];
for i = 1:No
N = [N,A_dash(:,i*(No-1)+i)];
end
Sigma_E = (diag([(N'*N)\(N'*X_e)]));
% Sigma_E(Sigma_E<0) = 0; % Since Sigma_E values cant be lesser than 0 
Sigma_E = abs(Sigma_E);
% SigmaE = (reshape((A'*A)\(A'*X_e),4,4))
SigmaE = (beta*SigmaE+(1-beta)*Sigma_E);
iter =iter +1

% pause
end
%% Using No of contraints = 5 , Errors were calculated (Didnt converge for no of constraints = 4,3) 
c = U(:,(n-2):n)'/L;
-(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
plot(log(diag(S).^2))