%% Data set initialisation
d1 = ...
[10	8.04
8	6.95
13	7.58
9	8.81
11	8.33
14	9.96
6	7.24
4	4.26
12	10.84
7	4.82
5	5.68];

d2= ...
[10	9.14
8	8.14
13	8.74
9	8.77
11	9.26
14	8.1
6	6.13
4	3.1
12	9.13
7	7.26
5	4.74];

d3 = ...
[10	7.46
8	6.77
13	12.74
9	7.11
11	7.81
14	8.84
6	6.08
4	5.39
12	8.15
7	6.42
5	5.73];

d4 = ...
[8	6.58
8	5.76
8	7.71
8	8.84
8	8.47
8	7.04
8	5.25
19	12.5
8	5.56
8	7.91
8	6.89];
%% Separating of dependent , independent variables
x1 = d1(:,1);
y1 = d1(:,2);
x2 = d2(:,1);
y2 = d2(:,2);
x3 = d3(:,1);
y3 = d3(:,2);
x4 = d4(:,1);
y4 = d4(:,2);
m = length(x1);

x1 = [x1,ones(m,1)];
x2 = [x2,ones(m,1)];
x3 = [x3,ones(m,1)];
x4 = [x4,ones(m,1)];
X2 = [d2(:,1).^2,x2];
%% finding fit
theta1 = (x1'*x1)\(x1'*y1);
theta2 = (x2'*x2)\(x2'*y2);
Theta2 = (X2'*X2)\(X2'*y2);
theta3 = (x3'*x3)\(x3'*y3);
theta4 = (x4'*x4)\(x4'*y4);
%% Errors
norm(y1-x1*theta1) % ERROR DATASET 1
norm(y2-x2*theta2) % ERROR DATASET 2
norm(y2-X2*Theta2) % ERROR DATASET 2 ,Quad fit
norm(y3-x3*theta3) % ERROR DATASET 3
norm(y4-x4*theta4) % ERROR DATASET 4
