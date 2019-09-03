function [ thdot ] = InitialConditionsV2(th, r, rdot,Param )
% Take the input rdot, r and ME to output thetadot
% At the point th = 0

L0 = Param.L0;
fr = Param.fr1;
k1 = Param.kA;
ME = Param.ME_0;
m = Param.m;
g = 9.81;

[ ~,thC,dCoP,~,dxc,dyc,~,~ ] = findCoP( th,Param,1 );

L = L0 + r;
dL = 1;

% Local position of CoP
xc =  fr*sin(thC);
yc = -fr*cos(thC) + fr;
% Mass 1
Y = -(-xc)*sin(th) + (L-yc)*cos(th);
% Single derivatives
dXdth = dCoP + sin(th)*xc - dxc*cos(th) - dyc*sin(th) + cos(th)*(L - yc);
dYdth = cos(th)*xc - dyc*cos(th) + dxc*sin(th) - sin(th)*(L - yc);
dXdr = dL*sin(th);
dYdr = dL*cos(th);

% Use a quadratic equation from ME = KE + PE
% a thdot^2 + b thdot + c = 0

a = dXdth^2 + dYdth^2;

b = 2*rdot*(dXdth*dXdr + dYdth*dYdr);

PE = 1/2*k1*r^2 + Y*m*g;

c = (dXdr^2 + dYdr^2)*rdot^2 + 2*(PE - ME)/m;

x = roots([a b c]);    

% if min(x) > 0
%     warning('Check the initial conditions')
%     r
%     rdot
%     ME
% end
thdot = max(x);

end

