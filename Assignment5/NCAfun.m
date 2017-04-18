function [ B_dash,P ] = NCAfun( Z )
%NCA using Alternating Least Squares
% Interfaced with fminunc
[U,S,V] = svd(Z);
k=3;
B = U(:,1:k);
B0 = B;
P = S(1:k,1:k)*V(:,1:k)';
%%
m = 7;
n= 3;
% B0 =[5,6,0,9,10,4,0,3,0,2,0,1,0,8,0,15,2,1,0,3,9 ];
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
B_dash = reshape(B,m,n);
P = (B_dash'*B_dash)\(B_dash'*Z);

err = rms(rms((B-B0)));
err1 = rms(rms((P-P0)));
i=i+1;
err
err1
% pause;
end

end

