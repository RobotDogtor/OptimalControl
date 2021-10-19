
clear
clc

%% 
n = 8;
m = 1;

%n size
x = (rand(n,1)-0.5)*10;
xr = (rand(n,1)-0.5)*10;
e = x-xr;
A = (rand(n,n)-0.5)*10;
Ar = diag(-rand(n,1)*10);
R1 = rand(n,n)*10;
P = lyap(Ar,R1);

%m size
r = (rand(m,1)-0.5)*5;
R2 = diag(rand(m,1)*10);

%both
B = (rand(n,m)-0.5)*5;
Br = (rand(n,m)-0.5)*10;

coeff1 = e'*P*B*inv(R2);
coeff2 = 2*e'*P*(A-Ar)*x - 2*e'*P*Br*r;
L2 = (coeff1\coeff2)';
test0 = coeff1*L2' - coeff2

coeff3 = 0.5*B*inv(R2)*L2';
coeff4 = (A-Ar)*x - Br*r;
test3 = e'*P*coeff3 - e'*P*coeff4
test4 = e'*P*(coeff3 - coeff4)
