clear 
clc
close all
format compact
format short g



%% system
n = 2;
m = 1;

vec = @(Mat) reshape(Mat,size(Mat,1)*size(Mat,2),1);
invvec = @(vec) reshape(vec,n,length(vec)/n);

%Real system: dx = Astar*x + Bstar*u
% Astar = (rand(n,n)-0.5)*10;
% Bstar = (rand(n,m)-0.5)*5;
Astar = [0 1;
         -2 -12];
bstar = 7;
Bstar = [0; bstar];
thAhat_star = vec(Astar);
thBhat_star = vec(bstar);

%reference system: dxr = ar*xr + br*r
wn = 4;
zeta = 0.5;
Ar = [0 1;
      -wn^2 -2*zeta*wn];
Br = [0 wn^2]';

%% Optimal Adaptive
%tuning
Qa = diag(ones(1,n*n)*5);
Qb = diag(ones(1,m)*5);
R1 = diag(ones(1,n)*1); %Lyapunov derivative for error dynamics
R2 = diag(ones(1,m)*100); %Weight Matrix for L = ... u'R2u ...
R2inv = inv(R2);

%Initial Values
x0 = [1 0]';
xr0 = [0 0]';
e0 = x0-xr0;
Ahat0 = Astar;
Bhat0 = bstar;
thAhat0 = vec(Ahat0);
thBhat0 = vec(Bhat0);
% solve are for X:  A'*X + X*A - X*B*X + C = 0   X = are(A, B, C)
P0 = lyap(Ar,R1);
V20 =(x0-xr0)'*P0*(x0-xr0) + (thAhat0-thAhat_star)'*inv(Qa)*(thAhat0-thAhat_star) +  (thBhat0-thBhat_star)'*inv(Qb)*(thBhat0-thBhat_star);
V30 = V20;
xs0 = [x0; xr0; thAhat0; thBhat0; V20; V30];
%time
t_end = 1;
tspan = [0 t_end];
%Simulate
opts = odeset('RelTol',1e-6);
tic
[t_ode,xs] = ode45(@(t,y) odefcn_optimal(t,y,Astar,Bstar,Ar,Br,Qa,Qb,R1,R2inv,n,m),tspan,xs0,opts);
toc
% Results
t_ode = t_ode';
x = xs(:,1:n)';
xr = xs(:,n+1:2*n)';
thAhat = xs(:,2*n+1:2*n+n*n)';
thBhat = xs(:,2*n+n*n+1:2*n+n*n+m)';
V2 = xs(:,end-1)';
V3 = xs(:,end)';
for i = 1:length(t_ode)
    P = lyap(Ar,R1);
    V1(i) = (x(:,i)-xr(:,i))'*P*(x(:,i)-xr(:,i)) + ...
            (thAhat(:,i)-thAhat_star)'*inv(Qa)*(thAhat(:,i)-thAhat_star) + ...
            (thBhat(:,i)-thBhat_star)'*inv(Qb)*(thBhat(:,i)-thBhat_star);
end
checkIfV123Failed(V1,V2,V3,0.0001)

[xd, dxd, ddxd] = xd_fcn(t_ode);    dt = 0.00000001;
%r= (dxd - Ar.*xd)./br;
% u = thx.*x + thr.*r;

% Plot
tplotrange = [0,t_end];
plotEveryNPoints = ceil(length(x)/100000);
plotGeneralAdaptiveSystemPerformance(t_ode,xr,xd,x,V1,V2,V3,thAhat,thBhat,vec(Astar),vec(bstar),tplotrange,plotEveryNPoints)

%% Functions 
function [xd,dxd,ddxd] = xd_fcn(t)
    %Desired Trajectory
    desx = @(tcurr) 1*sin(2*tcurr);% -2*cos(1.3*tcurr) + 5*sin(0.5*tcurr) + 0.5*cos(5*tcurr); 
    
    dt = 0.00001;
    xd = desx(t);
    xdnext = desx(t+dt);
    xdprev = desx(t-dt);
    dxd = (xdnext-xd)/dt;
    dxdprev = (xd-xdprev)/dt;
    ddxd = (dxd-dxdprev)/dt;
end

function dxs = odefcn_optimal(t,y,Astar,Bstar,Ar,Br,Qa,Qb,R1,R2inv,n,m)
    if t-round(t*2)/2<0.0001
        disp(['t = ' num2str(t)])
    end
    vec = @(Mat) reshape(Mat,size(Mat,1)*size(Mat,2),1);
    invvec = @(vec) reshape(vec,n,length(vec)/n);
    %identify simulation states
    x = y(1:n);  xr = y(n+1:2*n);  thAhat = y(2*n+1:2*n+n*n);  
    thBhat = y(2*n+n*n+1:2*n+n*n+m);  V2 = y(end-1);  V3 = y(end);
    Ahat = invvec(thAhat);  Bhat = thBhat;
    %find r
    [xd,dxd,ddxd] = xd_fcn(t);
    r = xd + dxd*Ar(2,2)/Ar(2,1) - ddxd/Ar(2,1);
    %find u
    if Bhat==0
        Bhat =0.1;
        thBhat = 0.1;
    end
    e = x-xr;
    S = [0 0; 0 2*Bhat*R2inv*Bhat];
    P = lyap(Ar,R1);
%     c1 = 2*e'*P*((Ahat-Ar)*x-Br*r-0.5*S*P*e);
    c1 = 2*e'*P*((Ahat-Ar)*x-Br*r);
    c2 = e'*P(:,2)*Bhat*R2inv;
    if abs(c2)<=0.000001
        L2T = c1/100000;
    else
        L2T = c1/c2;
    end
    u = -R2inv*(Bhat'*P(2,:)*e+0.5*L2T);
    %calculate derivatives
    dx = Astar*x + Bstar*u;
    dxr = Ar*xr+ Br*r;
    dthAhat = Qa*kron(x,eye(length(x)))*P*e;
    dthBhat = Qb*u*P(2,:)*e;
    
    thAhat_tilde = thAhat - vec(Astar);
    thBhat_tilde = thBhat - Bstar(2,1);
    %V1 is base lyap function after all is done, V3 from simple derivative,
    %V2 after substitution of adaptive laws into derivative
    dV2 = 2*e'*P*(Astar*x + Bstar*u - Ar*xr - Br*r) + 2*thAhat_tilde'*inv(Qa)*dthAhat + 2*thBhat_tilde'*inv(Qb)*dthBhat;
    dV3 = 2*e'*P*(Ahat*x + [0;Bhat]*u -Ar*xr-Br*r);
    dxs = [dx; dxr; dthAhat; dthBhat; dV2; dV3];
end