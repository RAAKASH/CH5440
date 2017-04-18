%% Loading data
clear all
load('Inorfull.mat')
pause
%% First samples of the 5
d = DATA(1:5:end,:)';% class convention
[U,S,V] = svd(d,'econ');
%% Scree plot to determine rank (approx)
close all;
plot(diag(S).^2);
pause
close all;
% Found to be close to 10
d_new = U(:,1:20)*S(1:20,1:20)*V(:,1:20)';%98.89 VARIANCE CAPTURED for n = 15
d_new(d_new<0) = 0; %Setting it to 0
[U,S,V] = svd(d_new,'econ');    
%% NMF
% Note CONC : Cr 0.0764 ,Ni 0.1965, Co .1720
k = 3;  
[W,H]=nmf(d_new,k,abs(U(:,1:k)*S(1:k,1:k)),abs(V(:,1:k)'),'verbose',0,'type','plain'); % Using sparse with alpha,beta = 0 is same as plain
% [W,H]= nnmf(d_new,3); % Matlab script to check accuracy
c1 = sum(sum(abs(d_new-W*H)));
c2 = max(max(abs(d_new-W*H)));
[~,n] = size(d_new);
real =  [PureNi',PureCr',PureCo'];
r = real - mean(real);
w = W -mean(W);
for i=1:k
    for j=1:k
    ex(i,j) = sum(w(:,i).*r(:,j))/(sum(w(:,i).^2)*sum(r(:,j).^2))^0.5;
    end
end

%% Avg data
d=0;
for i=1:25
d(i,1:176) = mean(DATA(((i-1)*5 + 1):((i-1)*5 + 5),:));
end
d=d';
[U,S,V]=svd(d);
d_new = U(:,1:15)*S(1:15,1:15)*V(:,1:15)';%98.89 VARIANCE CAPTURED
d_new(d_new<0) = 0; %Setting it to 0

%% NMF
% Note CONC : Cr 0.0764 ,Ni 0.1965, Co .1720
k = 3;
[W,H]=nmf(d_new,k,abs(U(:,1:k)*S(1:k,1:k)),abs(V(:,1:k)'),'verbose',0,'type','plain'); % Using sparse with alpha,beta = 0 is same as plain
% [W,H]=nnmf(d_new,3);
sum(sum(abs(d_new-W*H)));
max(max(abs(d_new-W*H)));

ext = (W\d_new);
[~,n] = size(d_new);
real =  [PureNi',PureCr',PureCo'];
r = real - mean(real);
w = W -mean(W);
for i=1:k
    for j=1:k
    ex1(i,j) = sum(w(:,i).*r(:,j))/(sum(w(:,i).^2)*sum(r(:,j).^2))^0.5;
    end
end
%% Q -2 
clear all;
load('ncadata.mat')
%%
k=3;
[U,S,V] = svd(measabs,'econ');
% Denoise
meas_new = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
[U,S,V] = svd(meas_new,'econ');
B = U(:,1:k);
% struct (U(:,1:k)*M) = struct (B)  
%  Structure= [1,1,0,1,1,1,0 
%              1,0,1,0,1,0,1
%              0,1,1,1,0,1,1 ]'
M = eye(3);
A = [B(3,2:3);B(7,2:3)];
b = -[B(3,1);B(7,1)];
M(2:3,1) = (A'*A)\(A'*b);
A = [B(2,[1,3]);B(4,[1,3]);B(6,[1,3])];
b = -[B(2,2);B(4,2);B(6,2)];
M([1,3],2) = (A'*A)\(A'*b);
A = [B(1,1:2);B(5,1:2)];
b = -[B(1,3);B(5,3)];
M(1:2,3) = (A'*A)\(A'*b); % Rotation matrix generated
alpha = -1;
Est_Pure_Spec = alpha*(M\(S(1:k,1:k)*V(:,1:k)')); % alpha =scale
e = Est_Pure_Spec' - mean(Est_Pure_Spec');
p = pureabs' -mean(pureabs');  
for i = 1:k
    for j = 1:k
     ex(i,j) = sum((e(:,i).*p(:,j)))/(sum(e(:,i).^2)*sum(p(:,j).^2))^0.5;
    end
end
k = 3;
[U,S,V] = svd(measabs,'econ');
meas_new1 = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';

[ B,P ] = NCAfun(meas_new1 );
Est_Pure_Spec1 = P ; % alpha =scale
e1 = Est_Pure_Spec1' - mean(Est_Pure_Spec1');
for i = 1:k
    for j = 1:k
     ex1(i,j) = sum((e1(:,i).*p(:,j)))/(sum(e1(:,i).^2)*sum(p(:,j).^2))^0.5;
    end
end