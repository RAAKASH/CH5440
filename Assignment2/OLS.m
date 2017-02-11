function [ theta1,RMSE1 ] = OLS( data_1,per)


std1 = std(data_1(:,end));
[m,n] = size(data_1);
dat_1 = (data_1 - mean(data_1))./(zeros(size(data_1))+std(data_1));

x = round(per*m); 

A1 = [dat_1(:,1:(n-1)) ,ones(m,1)];
b1 = dat_1(:,n);
theta1 = ((A1(1:x,:))'*(A1(1:x,:)))\(A1(1:x,:)' *b1(1:x));

RMSE1 = rms((A1((x+1):end,:)*theta1 -b1((x+1):end))*std1); % Scale of 10


end

