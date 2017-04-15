function [ theta1,RMSE1 ] = OLS( data_1,per)


std1 = std(data_1(:,1:(end-1)));
[m,n] = size(data_1);
dat_1 = (data_1(:,1:end-1) - mean(data_1(:,1:end-1)))./(std1);

x = round(per*m); 

% A1 = [dat_1(:,1:(n-1)) ,ones(m,1)];
% A2 = [data_1(:,1:(n-1)) ,ones(m,1)];
A1 = [dat_1(:,1:(n-1))];
b1 = data_1(:,n);
theta1 = (((A1(1:x,:))'*(A1(1:x,:)))\(A1(1:x,:)' *b1(1:x)))'./std1;

RMSE1 = rms(((data_1((x+1):end,1:(n-1))-mean(data_1(:,1:end-1)))*theta1' -b1((x+1):end))); 
% RMSE1 = 0;

end

