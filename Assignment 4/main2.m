%% Question 2 with IPCA
%% 2 -a Data matrix construction
 clear
 C = [];
 load('steamdatamiss.mat')
 clear std;
 [m,n] = size(Fmeas); 
 for j = 1:n
     f=0;
     for i=1:m
        if(isnan(Fmeas(i,j)))
         f=1;
         continue
        end
     end
     if(f==0)
     C = [C,Fmeas(:,j)];
     end
     
 end
 %%  ipca
 [m,n] = size(C); 
SigmaE = eye(m); % Initial guess of Error variances (Assuming Covariances  = 0 )
iter = 0;
beta = 0.4; % For convergence (relaxation parameter)
No = m; % No of constraints
Smaster = []; % Contains master error variances of all transformed data 

 
 [ Smaster,Sdmaster ] = IPCAGrand( C); % IPCA
%% Results
No = 11;
Sigma_E =diag(Sdmaster(No,:).^2);
L = diag(Sdmaster(No,:));
[U,S,V] = svd(L\C);
c = (U(:,1:end)'/(L));
c = c((m-No +1):end,:);
Reg =-(c(:,18:end)'*c(:,18:end))\(c(:,18:end)'*c(:,1:17));
Reg1 =-(Atrue(:,18:end)'*Atrue(:,18:end))\(Atrue(:,18:end)'*Atrue(:,1:17));
Res = abs(Reg1-Reg);
 sum(sum(Res))
 max(max((Res)))
 No_Const(Smaster)
 %% 2 -b Data matrix construction (Mean Imputed)
 [m,n] = size(Fmeas);
 t = isnan(Fmeas);
 for j = 1:n
     f=0;
     for i=1:m
        if(isnan(Fmeas(i,j)))
        C(i,j) =mean(Fmeas(t(i,:))); 
         continue
        else
        C(i,j) =Fmeas(i,j);
        end
     end
     
 end
 
 [ Smaster,Sdmaster ] = IPCAGrand( C); % IPCA
 %% Results
 No =11;
L = diag(Sdmaster(No,:));
[U] = svd(L\C);
c = (U(:,1:end)'/(L));
c = c((m-No +1):end,:);
 Reg =-(c(:,18:end)'*c(:,18:end))\(c(:,18:end)'*c(:,1:17));
Reg1 =-(Atrue(:,18:end)'*Atrue(:,18:end))\(Atrue(:,18:end)'*Atrue(:,1:17));
Res = abs(Reg1-Reg);
sum(sum(Res))
max(max(Res))
 No_Const(Smaster)
%% 2 -c PCA imputed
err = 10;
while(err>10^-15)
C0 =C;
[U,S,V]= svd(C,'econ');
C_new = U(:,1:17)*S(1:17,1:17)*(V(:,1:17))';
C(t) = C_new(t);
err = rms(rms((C - C0)./C0));
end
[ Smaster,Sdmaster ] = IPCAGrand( C); % IPCA
%%
No = 11;
 L = diag(Sdmaster(No,:));
[U,S,V] = svd(L\C);
c = (U(:,1:end)'/(L));
c = c((m-No +1):end,:);
Reg =-(c(:,18:end)'*c(:,18:end))\(c(:,18:end)'*c(:,1:17));
Reg1 =-(Atrue(:,18:end)'*Atrue(:,18:end))\(Atrue(:,18:end)'*Atrue(:,1:17));
Res = abs(Reg1-Reg);
sum(sum(Res))
max(max(Res))
 No_Const(Smaster)