function [ reg_1,RMSE_TLS_1,e1 ] = TLS( data_1,per)
%% Hopefully TLS

std1 = std(data_1(:,end));
[m,n] = size(data_1);
m1 = round(per*m); 
D_1 = (data_1 - mean(data_1))./(zeros(size(data_1))+std(data_1));
[U1,S1,V_1] =svd(D_1(1:m1,:));
reg_1 = V_1(:,end);
reg_1 = (-reg_1/reg_1(end));
reg_1 = reg_1(1:(end-1));
%% 
d_1 = D_1 + (-D_1*V_1(:,end))*(V_1(:,end)') ; % Transforming Data to new basis
X_1 = d_1(:,1:(n-1));
e1 = (X_1((m1+1):end,:)*reg_1-D_1((m1+1):end,n))*std1;
RMSE_TLS_1 = rms((X_1((m1+1):end,:)*reg_1-D_1((m1+1):end,n))*std1); % Scale of 10

end

