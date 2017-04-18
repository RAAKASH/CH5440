clear all;
load('ncadata.mat')
clc
Z = measabs;
%%
m = 7;
n=3;
B0 =[1,1,0,1,1,1,0,1,0,1,0,1,0,1,0,1,1,1,0,1,1 ];
B = B0;
%%
clear Aeq Error
err = 10;
i=0;
while(err>10^-100)
B_dash = reshape(B,m,n);
P = (B_dash'*B_dash)\(B_dash'*Z);

B0 = B;
P0 = P;
 fun = @(B)Obj(B,P,Z);
A = [];
b = [];
A_eq = [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0; ...
        0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0 ;...
        0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;...
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0];

b_eq =  zeros(7,1);
B = fmincon(fun,B,A,b,A_eq,b_eq);

err = rms(rms((B-B0)./B0));
err1 = rms(rms((P-P0)./P0));
i=i+1;
err
err1
pause;
end