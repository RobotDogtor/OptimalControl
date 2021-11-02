clear
clc

%% 
n = 2;
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

S = B*inv(R2)*B';

P2 = icare(A,S,R1)