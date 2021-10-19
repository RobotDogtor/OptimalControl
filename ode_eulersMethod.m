function [t_ode,xs] = ode_eulersMethod(derivativeFunction,tspan,xs0,dt)
t_ode = tspan(1):dt:tspan(end);
nt = length(t_ode);
n = length(xs0);
xs = zeros(nt,n);
xs(1,:) = xs0;
for i = 1:nt-1
    xs(i+1,:) = xs(i,:) + derivativeFunction(t_ode(i),xs(i,:))'*dt;
end

end