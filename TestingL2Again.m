
clear
clc
close all

for i = 1:10000
    %%
    n = 3;
    m = 4;

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
    Br = (rand(n,m))*5;

%     eig(Br*Br'*P^2)


    %% 
    xd = (rand(n,1)-0.5)*10;
    dxd = (rand(n,1)-0.5)*10;

    xddiff = (Ar*xd-dxd);
    Bi = 2*Br'*P'*(e)*xddiff';
    Bi = (e'*P*Br)'*xddiff'; %guarantees this term is pos
    Bi = (e'*P*Br)'*xddiff'*abs(2*e'*P*xddiff) *(1/abs(e'*P*Br*(e'*P*Br)')) *(1/(xddiff'*xddiff));
    Bi = (e'*P*Br)'*xddiff'*abs(2*e'*P*xddiff);
%     Bi = xddiff'*2*(abs(e'*P*xddiff)/(e'*P*Br*xddiff'*xddiff));
    % Bi = Br'*dot(e'*P,xddiff)
    % Bi = (1/dot(xddiff,xddiff))*xddiff'
    % Bi*xddiff
    result(i) = 2*e'*P*xddiff - e'*P*Br*Bi*xddiff;
    checking(i) = result(i)>0.0001;
    if(checking(i)>0)
        x=2;
    end

end

sum(checking)
% dot(e'*P,xddiff)
% 2*e'*P*xddiff
plot(1:length(result),result,'.')
ylim([-100 100])
