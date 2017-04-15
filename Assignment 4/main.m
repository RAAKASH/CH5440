clear all;
load('steamdata.mat');

%% Initialziation
[m,n] = size(Fmeas); 
SigmaE = eye(m); % Initial guess of Error variances (Assuming Covariances  = 0 )
iter = 0;
beta = 0.4; % For convergence (relaxation parameter)
No = m; % No of constraints
Smaster = []; % Contains master error variances of all transformed data 
[U,S,V] = svd(diag(std)\Fmeas,'econ');
%% 1-a
for j = 28:-1:1
No = j;
    for i = 1:100
     L = SigmaE.^0.5;     % L matrix (check notes)
     D = L\(Fmeas);
     [U,S,V] = svd(D,'econ');
    
     c = (U(:,1:end)'/(L));
     c = c((m-No +1):end,:);

    [sd,f] = stdest(c,Fmeas);
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
%% 1-b
No = 11;
Sigma_E =diag(Sdmaster(No,:).^2);
L = diag(Sdmaster(No,:));
[U,S,V] = svd(L\Fmeas);
c = (U(:,1:end)'/(L));
c = c((m-No +1):end,:);
Reg =-(c(:,18:end)'*c(:,18:end))\(c(:,18:end)'*c(:,1:17));
Reg1 =-(Atrue(:,18:end)'*Atrue(:,18:end))\(Atrue(:,18:end)'*Atrue(:,1:17));
Res = abs(Reg1-Reg);
 sum(sum(Res))
 max(max((Res)))
 %% 2 -a
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
 %% Auto scale
 Data =  (C'./std(C'))';
 [U,S,V] = svd(Data,'econ');
 const = U'/diag(std(C')) ;
 c = const(18:end,:);
 Reg =-(c(:,18:end)'*c(:,18:end))\(c(:,18:end)'*c(:,1:17));
Reg1 =-(Atrue(:,18:end)'*Atrue(:,18:end))\(Atrue(:,18:end)'*Atrue(:,1:17));
Res = abs(Reg1-Reg);
sum(sum(Res))
max(max(Res))
 %% 2 -b
 [m,n] =size(Fmeas);
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
  Data =  (C'./std(C'))';
 [U,S,V] = svd(Data,'econ');
 const = U'/diag(std(C')) ;
  c = const(18:end,:);
 Reg =-(c(:,18:end)'*c(:,18:end))\(c(:,18:end)'*c(:,1:17));
Reg1 =-(Atrue(:,18:end)'*Atrue(:,18:end))\(Atrue(:,18:end)'*Atrue(:,1:17));
Res = abs(Reg1-Reg);
sum(sum(Res))
max(max(Res))
%% 2 -c
C_dash =C;
C_new = C;
err = 10;
while(err>10^-10)
C0 = C;
[U,S,V]= svd(C,'econ');
C_new = U(:,1:17)*S(1:17,1:17)*(V(:,1:17))';
C(t) = C_new(t);
err = rms(rms((C - C0)./C0));
end

 Data =  (C'./std(C'))';
 [U,S,V] = svd(Data,'econ');
 const = U'/diag(std(C')) ;
 c = const(18:end,:);
 Reg =-(c(:,18:end)'*c(:,18:end))\(c(:,18:end)'*c(:,1:17));
Reg1 =-(Atrue(:,18:end)'*Atrue(:,18:end))\(Atrue(:,18:end)'*Atrue(:,1:17));
Res = abs(Reg1-Reg);
sum(sum(Res))
max(max(Res))