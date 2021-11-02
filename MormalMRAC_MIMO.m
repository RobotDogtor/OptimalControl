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
Bstar = [-1; -7];

%reference system: dxr = ar*xr + br*r
Ar = [0 1;
      -2 -8];
Br = [0 -1]';

%% Optimal Adaptive
%tuning
Qa = diag(ones(1,n*n)*0.5);
Qb = diag(ones(1,n*m)*0.5);
Rr1 = diag(ones(1,n)*10); %Lyapunov derivative for reference model control
Rr2 = diag(ones(1,m)*1);
Rr2inv = inv(Rr2);
R1 = diag(ones(1,n)*10); %Lyapunov derivative for error dynamics
R2 = diag(ones(1,m)*1); %Weight Matrix for L = ... u'R2u ...
R2inv = inv(R2);

%Initial Values
x0 = [1 0]';
xr0 = [2 0]';
e0 = x0-xr0;
Ahat0 = [0 1; -6 -7];
Bhat0 = [3 5]';
thAhat0 = vec(Ahat0);
thBhat0 = vec(Bhat0);
% solve are for X:  A'*X + X*A - X*B*X + C = 0   X = are(A, B, C)
P0 = are(Ar,2*Bhat0*inv(R2)*Bhat0',R1);
V20 =(x0-xr0)'*P0*(x0-xr0) + (thAhat0-thAhat_star)'*inv(Qa)*(thAhat0-thAhat_star) +  (thBhat0-thBhat_star)'*inv(Qb)*(thBhat0-thBhat_star);
V30 = V20;
xs0 = [x0; xr0; thAhat0; thBhat0; V20; V30];
%time
t_end = 10;
tspan = [0 t_end];
%Simulate
opts = odeset('RelTol',1e-6);
tic
[t_ode,xs] = ode45(@(t,y) odefcn_starController(t,y,Astar,Bstar,Ar,Br,Qa,Qb,Rr1,Rr2inv,R1,R2inv,n,m),tspan,xs0,opts);
toc
% Results
t_ode = t_ode';
x = xs(:,1:n)';
xr = xs(:,n+1:2*n)';
thAhat = xs(:,2*n+1:2*n+n*n)';
thBhat = xs(:,2*n+n*n+1:2*n+n*n+n*m)';
V2 = xs(:,end-1)';
V3 = xs(:,end)';
for i = 1:length(t_ode)
    P = are(Ar,2*invvec(thBhat(:,i))*inv(R2)*invvec(thBhat(:,i))',R1);
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
plotGeneralAdaptiveSystemPerformance(t_ode,xr,xd,x,V1,V2,V3,thAhat,thBhat,vec(Astar),vec(Bstar),tplotrange,plotEveryNPoints)

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

function dxs = odefcn_optimal(t,y,Astar,Bstar,Ar,Br,Qa,Qb,Rr1,Rr2inv,R1,R2inv,n,m)
    if t-round(t*2)/2<0.0001
        disp(['t = ' num2str(t)])
    end
    vec = @(Mat) reshape(Mat,size(Mat,1)*size(Mat,2),1);
    invvec = @(vec) reshape(vec,n,length(vec)/n);
    %identify simulation states
    x = y(1:n);  xr = y(n+1:2*n);  thAhat = y(2*n+1:2*n+n*n);  
    thBhat = y(2*n+n*n+1:2*n+n*n+n*m);  V2 = y(end-1);  V3 = y(end);
    Ahat = invvec(thAhat);  Bhat = invvec(thBhat);
    %find r
    [xd_,dxd_,ddxd_] = xd_fcn(t);
    xd = [xd_; dxd_];
    dxd = [dxd_; ddxd_];
%     Pr = are(Ar,2*Br*Rr2inv*Br',Rr1);
%     er = xr-xd;
%     cr1 = (er'*Pr*Br*Rr2inv);   cr2 = (2*er'*Pr*(Ar*xd-dxd));
    S = 2*Br*Rr2inv*Br';
    Pr = lyap(Ar,Rr1);
    er = xr-xd;
    cr1 = (er'*Pr*Br*Rr2inv);   cr2 = (2*er'*Pr*(Ar*xd-dxd - 0.5*S*Pr*er));
    if norm(cr1)<=0.001
        L2rT = zeros(1,size(Bhat,2));
    else
        L2rT = cr1\cr2;
    end
    r = -Rr2inv*(Br'*Pr*er + 0.5*L2rT);
    %find u
    e = x-xr;
    S = 2*Bhat*R2inv*Bhat';
    P = lyap(Ar,S,R1);
    c1 = Bhat*R2inv;
    c2 = 2*((Ahat-Ar)*x - Br*r);
    if norm(c1)<=0.0001
        L2T = zeros(1,size(Bhat,2))+100;
    else
        L2T = ([1 0]*c1\[1 0]*c2);
    end
    u = -0.5*R2inv*(2*Bhat'*P*e+L2T);
    %calculate derivatives
    dx = Astar*x + Bstar*u;
    dxr = Ar*xr+ Br*r;
    dthAhat = Qa*kron(x,eye(length(x)))*P*e;
    dthBhat = Qb*kron(u,eye(length(u)))*P*e;
    
    thAhat_tilde = thAhat - vec(Astar);
    thBhat_tilde = thBhat - vec(Bstar);
    %V1 is base lyap function after all is done, V3 from simple derivative,
    %V2 after substitution of adaptive laws into derivative
    dV2 = 2*e'*P*(Astar*x + Bstar*u - Ar*xr - Br*r) + 2*thAhat_tilde'*inv(Qa)*dthAhat + 2*thBhat_tilde'*inv(Qb)*dthBhat;
    dV3 = 2*e'*P*(Ahat*x + Bhat*u -Ar*xr-Br*r);
    dxs = [dx; dxr; dthAhat; dthBhat; dV2; dV3];
end

%%
function dxs = odefcn_starController(t,y,Astar,Bstar,Ar,Br,Qa,Qb,Rr1,Rr2inv,R1,R2inv,n,m)
    if t-round(t*2)/2<0.0001
        disp(['t = ' num2str(t)])
    end
    vec = @(Mat) reshape(Mat,size(Mat,1)*size(Mat,2),1);
    invvec = @(vec) reshape(vec,n,length(vec)/n);
    %identify simulation states
    x = y(1:n);  xr = y(n+1:2*n);  thAhat = y(2*n+1:2*n+n*n);  
    thBhat = y(2*n+n*n+1:2*n+n*n+n*m);  V2 = y(end-1);  V3 = y(end);
    Ahat = invvec(thAhat);  Bhat = invvec(thBhat);
    %find r
    [xd_,dxd_,ddxd_] = xd_fcn(t);
    xd = [xd_; dxd_];
    dxd = [dxd_; ddxd_];
    S = 2*Br*Rr2inv*Br';
    Pr = are(Ar,2*Br*Rr2inv*Br',Rr1);
    er = xr-xd;
%     cr1 = (er'*Pr*Br*Rr2inv);   cr2 = (2*er'*Pr*(Ar*xd-dxd));
%     
% %     Pr = lyap(Ar,Rr1);
% %     er = xr-xd;
% %     cr1 = (er'*Pr*Br*Rr2inv);   cr2 = (2*er'*Pr*(Ar*xd-dxd - 0.5*S*Pr*er));
%     if norm(cr1)<=0.001
%         L2rT = zeros(1,size(Bhat,2));
%     else
%         L2rT = cr1\cr2;
%     end
    xddiff = (Ar*xd-dxd);
    Bi = (er'*Pr*Br)'*xddiff'*abs(4*er'*Pr*xddiff)*(1/(er'*Pr*Br*(er'*Pr*Br)')) *(1/(xddiff'*xddiff));
    L2rT = inv(Rr2inv)*Bi*(Ar*xd-dxd);
    r = -Rr2inv*(Br'*Pr*er + 0.5*L2rT);

    %find u
    e = x-xr;
    S = 2*Bstar*R2inv*Bstar';
    P = are(Ar,S,R1);
    c1 = Bstar*R2inv;
    c2 = 2*((Astar-Ar)*x - Br*r);
    if norm(c1)<=0.0001
        L2T = zeros(1,size(Bstar,2))+1000;
    else
        L2T = ([1 0]*c1\[1 0]*c2);
    end
    u = -0.5*R2inv*(2*Bstar'*P*e+L2T);
    %calculate derivatives
    dx = Astar*x + Bstar*u;
    dxr = Ar*xr+ Br*r;
    dthAhat = 0*Qa*kron(x,eye(length(x)))*P*e;
    dthBhat = 0*Qb*kron(u,eye(length(u)))*P*e;
    
    thAhat_tilde = thAhat - vec(Astar);
    thBhat_tilde = thBhat - vec(Bstar);
    %V1 is base lyap function after all is done, V3 from simple derivative,
    %V2 after substitution of adaptive laws into derivative
    dV2 = 2*e'*P*(Astar*x + Bstar*u - Ar*xr - Br*r) + 2*thAhat_tilde'*inv(Qa)*dthAhat + 2*thBhat_tilde'*inv(Qb)*dthBhat;
    dV3 = 2*e'*P*(Ahat*x + Bhat*u -Ar*xr-Br*r);
    dxs = [dx*0; dxr; dthAhat*0; dthBhat*0; dV2; dV3];
end