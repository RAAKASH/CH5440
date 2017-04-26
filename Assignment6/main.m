clc ;
clear ;
load('vpdata.mat');
pause
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
%% Question 2
clc
clear all
load('arx.mat')
%% Constructing data Matrix
y = ymeas;
u = umeas;
D = [y(:,2:end);y(:,1:end-1);u(:,1:end-1)];
%% OLS Estimation
D = D';
b = D(:,1);
A = D(:,2:3);
reg1 = (A'*A)\(A'*b);
%% TLS Estimation
[U,S,V] = svd(D','econ');
reg2 = -U(end,1:end-1)/U(end,end);
%% Part D - Constructing stacked data
Dat  = [];
Dat1 = [];
[m,n] = size(y);
k  = 10; % Stack Size
for i=1:(n-k)
Dat  = [Dat,flipud(y(i:i+k)')];
Dat1 = [Dat1 ,flipud(u(i:i+k-1)')];
end

%% Assuming order = 1
% here k = m (from class notes)
ord  = 1; 
ind = k-ord+1; % No fo constarint equations
FinD = [Dat;Dat1];
[U,S,V] = svd(FinD,'econ'); 
Constraint = U(:,(end-ind+1):end)';
A = Constraint(:,1:k+1);
B = -Constraint(:,k+2:end);
[m,n] =size(A);
% Structure of A
A_class = zeros(m,n);
for i=1:n-1
A_class(i,i:i+1) = [1,1];
end
%Obtaing rotation matrix for A
for i=1:n-1
T = A(:,[1:(i-1),(i+2):n]);
b = -T(i,:)';
a = T([1:(i-1),(i+1):(m)],:)';
reg =  ((a'*a)\(a'*b))';
R = [reg(1:(i-1)),1,reg(i:m-1)];
M(i,:) = R;
end
A_new = M*A;
B_new = M*B;
a0 = mean(diag(A_new));
S=0;
for i =1:m
   S = S + ((A_new(i,i+1)));
end
a1 = S/(m-1);
b1 = mean(diag(B_new));
a1=-a1/a0;
b1 = -b1/a0;
a0 = -a0/a0;
% Compiling results
res = [[-1,reg1'];[-1,reg2];[a0,a1,b1]];

