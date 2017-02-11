%% Assignment 2 
%% Question 1
clc
clear
close
%% Data initialization /scaling
data_1 = xlsread('windedata.xlsx','White Wine','A2:L4899');
data_2 = xlsread('windedata.xlsx','Red Wine','A2:L1600');
dat_1 = (data_1 - mean(data_1))./(zeros(size(data_1))+std(data_1));
dat_2 = (data_2 - mean(data_2))./(zeros(size(data_2))+std(data_2));
D_1 = dat_1;
D_2 = dat_2;
[m1,n1] = size(data_1);
[m2,n2] = size(data_2);
%%  Scaling of data in the y column
dat_1(:,12) = dat_1(:,12)*std(data_1(:,12)) + mean(data_1(:,12));
dat_2(:,12) = dat_2(:,12)*std(data_2(:,12)) + mean(data_2(:,12));
%% OLS / TLS
%% OLS 
A1 = [dat_1(:,1:11) ,ones(m1,1)];
A2 = [dat_2(:,1:11) ,ones(m2,1)];
b1 = dat_1(:,12);
b2 = dat_2(:,12);
theta1 = ((A1(1:3430,:))'*(A1(1:3430,:)))\(A1(1:3430,:)' *b1(1:3430));
theta2 = ((A2(1:1200,:))'*(A2(1:1200,:)))\(A2(1:1200,:)' *b2(1:1200));

RMSE1 = rms(A1(3431:end,:)*theta1 -b1(3431:end)); % Scale of 10
RMSE2 = rms(A2(1201:end,:)*theta2 -b2(1201:end)); % Scale of 10

%% TLS - Ist attempt
% n=11;
% Sz_1 = (dat_1(1:3430,:)'*dat_1(1:3430,:));
% [V_1,D_1] =eig(Sz_1);
% reg_1 = V_1./V_1(end,:);
% reg_1 = -reg_1(1:11,12); 
% beta_1 =  mean(dat_1(1:3430,12));
% RMSE_TLS_1 = (sum((dat_1(3431:end,1:11)*reg_1-dat_1(3431:end,12)+beta_1).^2/m1))^0.5;
% D_1= diag(D_1);
% 
% Sz_2 = (dat_2(1:1200,:)'*dat_2(1:1200,:));
% [V_2,D_2] =eig(Sz_2);
% reg_2 = V_2./V_2(end,:);
% reg_2 = -reg_2(1:11,12); 
% beta_2 =  mean(dat_2(1:1200,12));
% RMSE_TLS_2 = (sum((dat_2(1201:end,1:11)*reg_2-dat_2(1201:end,12)+beta_2).^2/m1))^0.5;
% D_2= diag(D_2);
%% TLS - 2nd attempt

std1 = std(data_1(:,12));
mean1 = mean(data_1(:,12));
[U_1,S_1,V_1] =svd(D_1(1:3430,:));
S_1 = diag(S_1); %
reg_1 = V_1(:,end);
reg_1 = (-reg_1/reg_1(end));
reg_1 = reg_1(1:(end-1));
d_1 = D_1 + (-D_1*V_1(:,end))*(V_1(:,end)') ; % Transforming Data to new basis
X_1 = d_1(:,1:(end-1));

RMSE_TLS_1 = rms((X_1(3431:end,:)*reg_1-D_1(3431:end,12))*std1); % Scale of 10
%RMSE_TLS_1 = rms(X_1(3431:end,:)*reg_1-D_1(3431:end,12)); 


% Test to check if TLS algo has been implemented properly
% RMSE_TLS_1 = (sum((X_1(1:3430,:)*reg_1-d_1(1:3430,end)).^2/m1))^0.5;

std2 = std(data_2(:,12));
[U_2,S_2,V_2] =svd(D_2(1:1200,:));
S_2 = diag(S_2); %
reg_2 = V_2(:,end);
reg_2 = (-reg_2/reg_2(end));
d_2 = D_2+(-D_2*V_2(:,end))*(V_2(:,end)') ;
X_2 = d_2(:,1:(end-1));
reg_2 = reg_2(1:(end-1));
RMSE_TLS_2 = rms((X_2(1201:end,:)*reg_2-D_2(1201:end,12))*std2); % Scale of 10
%RMSE_TLS_2 = rms(X_2(1201:end,:)*reg_2-D_2(1201:end,12));

% Test to check if TLS algo has been implemented properly
% RMSE_TLS_2 = (sum((X_2(1:1200,:)*reg_2-d_2(1:1200,end)).^2/m1))^0.5;
pause;
%% Question 2

clc
clear
close
%% Data initialization /scaling
data_1 = xlsread('temperature_global.xlsx','temperature_global','C8:G38');
[m1,n1] = size(data_1);
dat_1 = (data_1 - mean(data_1))./(zeros(size(data_1))+std(data_1));
% Correlation 
sum(dat_1(:,1).*dat_1(:,end))
sum(dat_1(:,2).*dat_1(:,end))
sum(dat_1(:,3).*dat_1(:,end))
sum(dat_1(:,4).*dat_1(:,end))
%%  Scaling of data in the y column
[reg_1,RMSE_1,e1] = TLS(data_1,0.7);
[reg_1_OLS,RMSE_1_OLS] = OLS(data_1,0.7);

%% Question 3

clc
clear
close
%% 
S = [7,21,34;21,64,102;34,102,186];
[V,D] = eig(S);
D = diag(D);
%% 2nd question 
(D(3))/sum(D) > 0.95 
%% Constraint equations
[U,Z,V] = svd(S);
U(:,2:3)'
%% Scores
Score= [10.1 , 73, 135.5]*U(:,1);
%%  Solving for mass , length using two equations
A = U(:,2:3)';
Ans = -[A(:,1),A(:,3)]\((A(:,2)*73))
%%
B = U(:,3);
mass = -[73,135.5]*B(2:3)/(B(1));

