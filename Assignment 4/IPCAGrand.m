function [ Smaster,Sdmaster ] = IPCAGrand( C)
%Finds the error variances for all posssible values for the number of
%contraints
%
[m,n]=size(C);
SigmaE = eye(m); % Initial guess of Error variances (Assuming Covariances  = 0 )
iter = 0;
beta = 0.4; % For convergence (relaxation parameter)
No = m; % No of constraints
Smaster = []; % Contains master error variances of all transformed data 

for j = 28:-1:1
No = j;
    for i = 1:100
     L = SigmaE.^0.5;     % L matrix (check notes)
     D = L\(C);
     [U,S,V] = svd(D,'econ');
    
     c = (U(:,1:end)'/(L));
     c = c((m-No +1):end,:);

    [sd,f] = stdest(c,C);
    if(f==0)
        break;
    end
    SigmaE = diag(sd.^2);
    i
   
    end
    if(f==0)
        break;
    end
   Smaster(j,:) = diag(S).^2'/n; 
   Sdmaster(j,:) = sd';
end

end

