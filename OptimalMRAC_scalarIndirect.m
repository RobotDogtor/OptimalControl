clear 
clc
close all
format compact
format short g

%scalar system: dx = astar*x + bstar*u
J = 0.0026;
B = 0.00057;
L = 0.0045;
R = 0.5;
K = 0.56;
astar = -(B*R+K^2)/(R*J);
bstar = K/(R*J);

astar = astar;
bstar = bstar;

%reference system: dxr = ar*xr + br*r
ar = -30;
br = 10;

%From compatability equations:
thx_star = (ar-astar)/bstar;
thr_star = br/bstar;

%% nonOptimal adaptive control
gam_x = 11;
gam_r = 12;

x0 = 3;
xr0 = 3;
thx0 = 0;
thr0 = 0;
V20 =0.5*(x0-xr0)^2 + (abs(bstar)/(2*gam_x))*(thx0-thx_star)^2 +  (abs(bstar)/(2*gam_r))*(thr0-thr_star)^2;
V30 = V20;
xs0 = [x0 xr0 thx0 thr0 V20 V30];
%time
t_end = 10;
tspan = [0 t_end];
%Simulate
opts = odeset('RelTol',1e-7);
tic
[t_ode,xs] = ode45(@(t,y) odefcn_nonoptimal(t,y,astar,bstar,ar,br,gam_x,gam_r,thx_star,thr_star),tspan,xs0,opts);
toc
% Results
t_ode = t_ode';
x = xs(:,1)';
xr = xs(:,2)';
thx = xs(:,3)';
thr = xs(:,4)';
V2 = xs(:,5)';
V3 = xs(:,6)';
V1 = 0.5*(x-xr).^2 + (abs(bstar)/(2*gam_x))*(thx-thx_star).^2 +  (abs(bstar)/(2*gam_r))*(thr-thr_star).^2;
checkIfV123Failed(V1,V2,V3,0.0001)

xd = xd_fcn(t_ode);    dt = 0.00000001;
dxd = (xd_fcn(t_ode+dt)-xd)./dt;
r= (dxd - ar.*xd)./br;
u = thx.*x + thr.*r;

% Plot
tplotrange = [0,t_end];
plotAdaptiveSystemPerformance(t_ode,xr,xd,x,V1,V2,V3,thx,thr,thx_star,thr_star,tplotrange,1)


%% Optimal Adaptive
gam_x = 1000;
gam_r = 100;

x0 = 3;
xr0 = 3;
ahat0 = -100;
bhat0 = 100; %if bhat goes 0, fails need projection
Q = 300;
P = lyap(ar,Q);
R2 = 10;
V20 =(x0-xr0)^2*P + (1/(gam_x))*(ahat0-astar)^2 +  (1/(gam_r))*(bhat0-bstar)^2;
V30 = V20;
xs0 = [x0 xr0 ahat0 bhat0 V20 V30];
%time
t_end = 10;
tspan = [0 t_end];
%Simulate
opts = odeset('RelTol',1e-6);
tic
[t_ode,xs] = ode45(@(t,y) odefcn_optimal(t,y,astar,bstar,ar,br,gam_x,gam_r,inv(R2),P),tspan,xs0,opts);
toc
% Results
t_ode = t_ode';
x = xs(:,1)';
xr = xs(:,2)';
ahat = xs(:,3)';
bhat = xs(:,4)';
V2 = xs(:,5)';
V3 = xs(:,6)';
for i = 1:length(t_ode)
    V1_(i) = (x(:,i)-xr(:,i))'*P*(x(:,i)-xr(:,i)) + (ahat(:,i)-astar)'*(1/gam_x)*(ahat(:,i)-astar) +  (bhat(:,i)-bstar)'*(1/gam_r)*(bhat(:,i)-bstar);
end
checkIfV123Failed(V1_,V2,V3,0.0001)

xd = xd_fcn(t_ode);    dt = 0.00000001;
dxd = (xd_fcn(t_ode+dt)-xd)./dt;
r= (dxd - ar.*xd)./br;
% u = thx.*x + thr.*r;

% Plot
tplotrange = [0,t_end];
plotEveryNPoints = ceil(length(x)/100000);
plotAdaptiveSystemPerformance(t_ode,xr,xd,x,V1_,V2,V3,ahat,bhat,astar,bstar,tplotrange,plotEveryNPoints)

%% Functions 
function xd = xd_fcn(t)
    %Desired Trajectory
    xd = 1*sin(2*t) -2*cos(1.3*t) + 5*sin(0.5*t) + 0.5*cos(5*t); 
end

function dxs = odefcn_nonoptimal(t,y,astar,bstar,ar,br,gam_x,gam_r,thx_star,thr_star)
    %find r
    dt = 0.000001;
    xdcurr = xd_fcn(t);
    xdnext = xd_fcn(t+dt);
    dxd = (xdnext-xdcurr)/dt;
    r= (dxd - ar*xdcurr)/br;
    
    %identify simulation states
    x = y(1);  xr = y(2);  thx = y(3);  thr = y(4);  V2 = y(5);  V3 = y(6);
    
    %compute derivatives
    u = thx*x + thr*r;
    e = x-xr;
    dx = astar*x + bstar*u;
    dxr = ar*xr+ br*r;
    dthx = -sign(bstar)*gam_x*x*e;
    dthr = -sign(bstar)*gam_r*r*e;
    thxti = thx - thx_star;
    thrti = thr - thr_star;
    dV2 = ar*e^2 + thxti*abs(bstar)*(sign(bstar)*e*x+dthx/gam_x) + ...
        thrti*abs(bstar)*(sign(bstar)*e*r+dthr/gam_r);
    dV3 = e*(dx-dxr) + abs(bstar)*thxti*dthx/gam_x + abs(bstar)*thrti*dthr/gam_r;
    dxs = [dx; dxr; dthx; dthr; dV2; dV3];
end

function dxs = odefcn_optimal(t,y,astar,bstar,ar,br,gam_x,gam_r,R2inv,P)
    %find r
    dt = 0.0001;
    xdcurr = xd_fcn(t);
    xdnext = xd_fcn(t+dt);
    dxd = (xdnext-xdcurr)/dt;
    r= (dxd - ar*xdcurr)/br;
    
    %identify simulation states
    x = y(1);  xr = y(2);  ahat = y(3);  bhat = y(4);  V2 = y(5);  V3 = y(6);
    
    %compute derivatives
    e = x-xr;
    diff = (ahat-ar)*x-br*r;
    u = -diff/bhat;
    dx = astar*x + bstar*u;
    dxr = ar*xr+ br*r;
    dahat = gam_x*x*P*e;
    dbhat = gam_r*u*P*e;
    
    ahat_tilde = ahat - astar;
    bhat_tilde = bhat - bstar;
    %V1 is base lyap function after all is done, V3 from simple derivative,
    %V2 after substitution of adaptive laws into derivative
    dV2 = 2*e*P*(astar*x+bstar*u-ar*xr-br*r) + 2*ahat_tilde'*(1/gam_x)*dahat + 2*bhat_tilde'*(1/gam_r)*dbhat;
    dV3 = 2*e*P*(ahat*x + bhat*u -ar*xr-br*r);
    dxs = [dx; dxr; dahat; dbhat; dV2; dV3];
end