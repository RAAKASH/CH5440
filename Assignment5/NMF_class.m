function [ B,P ,err ] = NMF_class( A )
%Executes NMF
[U,S,V] = svd(A,'econ');
 B = abs(U*S);% First Iteration
 P = abs(V'); % Second Iteration
[m,n] = size(B);
[m1,n1] = size(P);
B1 = B; % Temp matrices
P1 = P; % Temp matrices
err = 10;
while(err>10^-15)
    Z1 = B1*P1*P1';
    Z  = B1'*B1*P1; 
    for p=1:m1
        for j=1:n1
            P1(p,j) = P(p,j)*(sum(B(:,p).*A(:,j))/Z(p,j));
        end
    end
    for i=1:m
        for j=1:n
            B1(i,j) = B(i,j)*(sum(A(i,:).*P1(j,:))/Z1(i,j));
        end
    end
     err = rms(rms(P-P1)) + rms(rms(B-B1))
    
    B = B1;
    P = P1;
    
     pause
end

end

