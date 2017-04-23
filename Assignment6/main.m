clc ;
clear ;
load('vpdata.mat');
%% 
m = mean(temp);
s = std(temp) ; 
data = (temp(1:70) - m)/s;
w = 10^.5;  % SCale in gaussian kernel
%% 
n = length(data);
for i=1:70
    for j=1:70
        phi(i,j) = gauss( data(i),data(j),w );
    end
end
%% Reorienting data

[U,lamda] = eig(phi);
U = fliplr(U);
lamda = (diag(fliplr(diag(lamda)')));
%% Acquiring data for PRESS
 for i = 71:100
    for j = 1:70
     phi1(i-70,j) = gauss( (temp(i)-m)/s,data(j),w );
    end
 end
for k = 1:n
T  = (lamda(1:k,1:k)^0.5)\(U(:,(1:k))'*phi');
beta = (psat(1:70)'*T')/(T*T');
T1 = (lamda(1:k,(1:k))^0.5)\U(:,(1:k))'*phi1';
err(k,:) = (psat(71:100)'- beta*T1);
% pause
end
Err =((sum(err.^2,2)/30).^0.5);
no = find(Err(1:21) == min(Err(1:21)));

%% PLotting Press
close all
plot(Err);
xlabel('No of Principal Components');
ylabel('RMSE');
title('RMSE vs PCs');
pause;
close all
plot(err(no,:));
xlabel('Last 30 samples');
legend('For the no of PC corresponding to the least RMSE');
ylabel('RMSE');
title('RMSE vs last 30 samples');
pause;
min(Err)
max(abs(err(no,:)))
%% 55 ,100
clear err1
clear phi2
t = [40:55,56:80]';
% t = [55,100]';
n = length(t);
for i = 1:n
    for j = 1:70
     phi2(i,j) = gauss( (t(i)-m)/s,data(j),w );
    end
end
k = no;
T  = (lamda(1:k,1:k)^0.5)\(U(:,(1:k))'*phi');
beta = (psat(1:70)'*T')/(T*T');
T1 = (lamda(1:no,1:no)^0.5)\U(:,1:no)'*phi2';
err1(1,:) = (p(t)'- beta*T1)
close all;
plot(t,beta*T1,40:80,p(40:80))
xlabel('Temperature')
ylabel('P_{sat}');
title('Modeling comparison');
legend('Modelled','True');
%% 
 