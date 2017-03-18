%% Notes
% Note : When using data representation of convention used in class (d) ,
% covariance matrix is found by cov(d').

%-------------------------------------------------------------------------------------------------%
%% Variables loading
clear ;

close all;
load('Inorfull.mat');
d = [DATA,CONC]' ; % Class convention if transpose is used
[U,S,V] = svd(d);
test = zeros(size(d));
for i =1:length(diag(S))
test = test +U(:,i)*S(i,i)*V(:,i)';
end
rms(rms(test-d))
% pause;
%% LOOCV 1-a
X = d(1:(end-3),:);
[~,S,~] = svd(X);

s = diag(S);
n= length(s);
[r,m] =size(d);     % No of rows,columns
err = zeros(m,n); %LOOCV for different number of constraints.
ERROR = zeros(1,n); % Final set of Errors using LOOCV for different number of errors.
%% n constaints to 3 constraints 1-a
clc;
fprintf('Computing CV of original data \n');
for j = 1:m  % j represents the CV 
  X =  [d(1:(end-3),1:(j-1)),d(1:(end-3),(j+1):end)]'; % X data contra class convention for regression
  
  [U,S,V] = svd(X);
  [mu,nu] = size(U);

for i = 1:(n) % i represents the no of scores/PCs taken into consideration
  A =  U;
  B =  V(:,1:i)';
  T = A*S;
  T = T(:,1:i);
  b =  [d((end-2):end,1:(j-1)),d((end-2):end,(j+1):end)]'; %Y data for training
  constraints = (T'*T)\(T'*b); % Obtaining regression model
  err(j,i) = rms(d((end-2):end,j)'-d(1:(end-3),j)'*B'*constraints); %Predicting
  fprintf('No of constraints : %d \n',i);
end
end
for k = 1:(n)    % matrix becomes singular at no of scores =30
ERROR(k) = rms(err(:,k)); 
end
fprintf('Minimum error for unaltered data \n'); 
min(ERROR(1:end))
find(ERROR==min(ERROR(1:end)))
xlswrite('Error file',err,'Alldata');
xlswrite('Mean Error file',ERROR,'Alldata');
%% Plotting  1-a
fprintf('Printing Graphs');
close all;
plot(1:(n-1),ERROR(1:(n-1)));
xlabel('No of PCs');
ylabel('Errors');
title('LOOCV vs No of PCs (1-a)');
print('LOOCV vs No of PCs (1-a)','-dpng');
pause
%%   1-b ,Finding new set of data by averaging
clc
fprintf('Averaged Data set computation \n');
summation=zeros(r,1);
d_new = zeros(r,m/5);

for i=1:m
    summation = summation + d(:,i);
    if(mod(i,5)==0)
    d_new(:,i/5) = summation/5;
    summation=zeros(r,1);
    end
end
X =  d_new(1:(end-3),:);
[~,S,~] = svd(X);
n = length(diag(S));
err = zeros(m/5,n);
%% Computing LOOCV 1-b
fprintf('LOOCV on Averaged Data set  \n');
for j = 1:m/5
     
  X =  [d_new(1:(end-3),1:(j-1)),d_new(1:(end-3),(j+1):end)]';
  [U,S,V] = svd(X);
  [mu,nu] = size(U);

for i = 1:(n)
 A =  U;
  B =  V(:,1:i)';
  T = A*S;
  T = T(:,1:i);
  b =  [d_new((end-2):end,1:(j-1)),d_new((end-2):end,(j+1):end)]';
  constraints = (T'*T)\(T'*b);
  err(j,i) = rms(d_new((end-2):end,j)'-d_new(1:(end-3),j)'*B'*constraints);
  fprintf('No of constraints : %d \n',i);
end
end
for k = 1:(n)
ERROR1(k) = rms(err(:,k)); 
end
fprintf('Min Error \n');

min(ERROR1)
find(ERROR1==min(ERROR1(1:end)))
xlswrite('Error file',err,'Avgdata');
xlswrite('Mean Error file',ERROR1,'Avgdata');
%% Plotting 1-b
fprintf('Plotting graph of LOOCV \n');
close all;
plot(1:n,ERROR1);
xlabel('No of PCs');
ylabel('Errors');
title('LOOCV vs No of PCs (1-b)');
print('LOOCV vs No of PCs (1-b)','-dpng');
pause
%%   1-c ,Finding standard deviation
clc
fprintf('Computing standard deviation \n');
dev = zeros(r-3,m/5);
for i=1:(m/5)
    dev(:,i) = std(d(1:(end-3),((i-1)*5+1):(i*5))')';
end

%% Plotting standard deviations for various wavelengths , mixtures 1-c
close all;
clc ;
xaxis =298+2*(1:(r-3)); 
plot(xaxis,dev(:,1),'r',xaxis,dev(:,2),'b',xaxis,dev(:,3),'y',xaxis,dev(:,4),'g');
xlabel('Wavelengths in nm');
ylabel('Standard deviations of absorbance matrix');
legend('Mixture -1','Mixture-2','Mixture-3','Mixture-4');
title('Std vs Wavelength');
fprintf('It can be seen that the standard deviations are of the same \n trend for all mixtures but vary with wavelengh');
pause;
close all;
fprintf('\n \n Hence averaging deviation over all mixtures ,printing new graph');
stddev = mean(dev')';
plot(xaxis,stddev,'r');
xlabel('Wavelengths in nm');
ylabel('Standard deviations of absorbance matrix');
legend('Averaged');
title('Std vs Wavelength (Generic trend)');

stddev(r-2:r,1) = 1;
pause;
%% Scaling abosrbance matrix 1-c
d_std = ((d_new')./(stddev'))';
X =  d_std(1:(end-3),:);
[U,S,V] = svd(X);
n = length(diag(S));
%% 1-c
for j = 1:m/5
     
  X =  [d_std(1:(end-3),1:(j-1)),d_std(1:(end-3),(j+1):end)]';
  [U,S,V] = svd(X);
  [mu,nu] = size(U);

for i = 1:(n)
 A =  U;
  B =  V(:,1:i)';
  T = A*S;
  T = T(:,1:i);
  b =  [d_std((end-2):end,1:(j-1)),d_std((end-2):end,(j+1):end)]';
  constraints = (T'*T)\(T'*b);
  err(j,i) = rms(d_std((end-2):end,j)'-d_std(1:(end-3),j)'*B'*constraints);
fprintf('No of constraints : %d \n',i);
end
end
for k = 1:(n)
ERROR2(k) = rms(err(:,k)); 
end
fprintf('Min error \n');
min(ERROR2)
find(ERROR2 == min(ERROR2))
xlswrite('Error file',err,'stddata');
xlswrite('Mean Error file',ERROR2,'stddata');
%% Plotting  1-c
close all;
plot(1:(n),ERROR2);
xlabel('No of PCs');
ylabel('Errors');
title('LOOCV vs No of PCs (1-c)');
print('LOOCV vs No of PCs (1-c)','-dpng');
pause
%%  1-a using the First datapoints of the 130 samples
d_1a = zeros(size(d_std));
for i =1:m/5
d_1a(:,i) = d(:,(i-1)*5+1); 
end
X = d_1a(1:(end-3),:);
[~,S,~] = svd(X);

s = diag(S);
n= length(s);

err = zeros(m/5,n); %LOOCV for different number of constraints.

%% n constaints to 3 constraints 1-a
clc;
fprintf('Computing CV of original data using the first data point of the 5samples alone \n');
for j = 1:m/5  % j represents the CV 
  X =  [d_1a(1:(end-3),1:(j-1)),d_1a(1:(end-3),(j+1):end)]'; % X data contra class convention for regression
  
  [U,S,V] = svd(X);
  [mu,nu] = size(U);

for i = 1:(n) % i represents the no of scores/PCs taken into consideration
  A =  U;
  B =  V(:,1:i)';
  T = A*S;
  T = T(:,1:i);
  b =  [d_1a((end-2):end,1:(j-1)),d_1a((end-2):end,(j+1):end)]'; %Y data for training
  constraints = (T'*T)\(T'*b); % Obtaining regression model
  err(j,i) = rms(d_1a((end-2):end,j)'-d_1a(1:(end-3),j)'*B'*constraints); %Predicting
  fprintf('No of constraints : %d \n',i);
end
end
for k = 1:(n)    % matrix becomes singular at no of scores =30
ERROR1a(k) = rms(err(:,k)); 
end
fprintf('Minimum error for unaltered data \n'); 
min(ERROR1a(1:end))
find(ERROR1a==min(ERROR1a(1:end)))
xlswrite('Error file',err,'Firstdata');
xlswrite('Mean Error file',ERROR1a,'Firstdata');
%% Plotting  1-a
fprintf('Printing Graphs');
close all;
plot(1:(n-1),ERROR1a(1:(n-1)));
xlabel('No of PCs');
ylabel('Errors');
title('LOOCV vs No of PCs (1-a)');
print('LOOCV vs No of PCs (1-a) final','-dpng');
pause

%--------------------------------------------------------------------------------------------------%
%% Question 2 - 1 data 1
clc;
clear;
close all;
load('flowdata1.mat');
clear std
%% Measured data
fprintf('Flow data 1 \n');
data = Fmeas;
[U,S,~] = svd((data)');
s = diag(S);
n = length(s);
c = U(:,(n-2):n)'; % 2 independent variables F1,F2
constraints_PCA_meas1 = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
c = Atrue;
constraints_real = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
plot(1:length(s),s)
xlabel('ith Singular value');
ylabel('Singular Values');
title('Singular values plot');
print('Singular value plot 2-a ,1','-dpng');
pause
%% Real data
data = Ftrue;
[U,S,V] = svd((data)'); % Try mean centering the data , The results seem better when it isnt cenetered
s = diag(S);
n = length(s);
c = U(:,(n-2):n)'; % 2 independent variables F1,F2
constraints_PCA_true = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
c = Atrue;
constraints_real = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
plot(1:length(s),s);
clc
%% Auto scale of measured data
fprintf('Auto scaled measurements of flow 1 \n');
data = (Fmeas-mean(Fmeas))./std(Fmeas); 
[U,S,~] = svd(data');
s = diag(S);
n = length(s);
c = U(:,(n-2):n)'; % 2 independent variables F1,F2
constraints_PCA_meas_auto1 = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
c = Atrue;
constraints_real = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
plot(1:length(s),s)
xlabel('ith Singular value');
ylabel('Singular Values');
title('Singular values plot');
print('Singular value plot 2-b,1','-dpng');
A=constraints_PCA_meas_auto1*Ftrue(:,1:2)';
error = sum(sum(Ftrue(:,3:end) - A'));% As can be seen the results are terrible
pause  
%% Finding the independent variables using condition number
data = Fmeas;
[U,S,~] = svd(data');
c = U(:,(n-2):n)';
[~,e,~] = svd(c(:,3:end));
e = diag(e);
condition_number1 = max(e)/min(e);
[~,e,~] = svd(c(:,1:3));
e = diag(e);
condition_number2 = max(e)/min(e);
[~,e,~] = svd(c(:,2:4));
e = diag(e);
condition_number3 = max(e)/min(e);

% Run for all combinations
v = nchoosek(1:5,3);
for i =1:10
[~,e,~] = svd(c(:,v(i,:)));
e = diag(e);
condition_number(i) = max(e)/min(e);
end
Dependent_var = v(condition_number ==min(condition_number),:);

%% flow data 2
clc;
clear;
close all;
load('flowdata2.mat');
clear std
%% Measured data
% data = Fmeas./std(Fmeas); 
data = Fmeas;
[U,S,~] = svd(data');
s = diag(S);
n = length(s);
c = U(:,(n-2):n)'; % 2 independent variables F1,F2
constraints_PCA_meas2 = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
c = Atrue;
constraints_real = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
plot(1:length(s),s)
xlabel('ith Singular value');
ylabel('Singular Values');
title('Singular values plot');
print('Singular value plot 2-a ,2','-dpng');
pause
%% scaled data
 data = (Fmeas-mean(Fmeas))./std(Fmeas); 
[U,S,~] = svd(data');
s = diag(S);
n = length(s);
c = U(:,(n-2):n)'; % 2 independent variables F1,F2
constraints_PCA_meas_auto2 = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
c = Atrue;
constraints_real = -(c(:,3:end)'*c(:,3:end))\(c(:,3:end)'*c(:,1:2))
plot(1:length(s),s)
xlabel('ith Singular value');
ylabel('Singular Values');
title('Singular values plot');
print('Singular value plot 2-b ,2','-dpng');
pause
%% Finding the independent variables using condition number
data = Fmeas;
[U,S,~] = svd(data');
n =length(diag(S));
c = U(:,(n-2):n)';
[~,e,~] = svd(c(:,3:end));
e = diag(e);
condition_number1 = max(e)/min(e);
[~,e,~] = svd(c(:,1:3));
e = diag(e);
condition_number2 = max(e)/min(e);
[~,e,~] = svd(c(:,2:4));
e = diag(e);
condition_number3 = max(e)/min(e);
% Run for all combinations
v = nchoosek(1:5,3);
for i =1:10
  ind =  (1:5~=v(i,1)).*(1:5~=v(i,2)).*(1:5~=v(i,3));
  [~,e,~] = svd(-(c(:,v(i,:))'*c(:,v(i,:)))\(c(:,v(i,:))'*c(:,ind==1)));
  e = diag(e);
condition_number(i) = max(e)/min(e);
end
Dependent_var = v(condition_number ==min(condition_number),:);

