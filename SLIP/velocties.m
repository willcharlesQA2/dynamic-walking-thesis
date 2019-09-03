function [ Xdot,Ydot ] = velocties( th,r,thdot,rdot,Param )
%Calculate Xdot and Ydot for values of th,r,thdot,rdot etc.

L0 = Param.L0;
fr = Param.fr;


[ ~,thC,dCoP,~,dxc,dyc,~,~ ] = findCoP( th,Param );

L = L0 + r;
dL = 1;

% Local position of CoP
xc =  fr*sin(thC);
yc = -fr*cos(thC) + fr;
% Mass 1
% Y = -(-xc)*sin(th) + (L-yc)*cos(th);
% Single derivatives
dXdth = dCoP + sin(th)*xc - dxc*cos(th) - dyc*sin(th) + cos(th)*(L - yc);
dYdth = cos(th)*xc - dyc*cos(th) + dxc*sin(th) - sin(th)*(L - yc);
dXdr = dL*sin(th);
dYdr = dL*cos(th);

Xdot = dXdth*thdot + dXdr*rdot;
Ydot = dYdth*thdot + dYdr*rdot;
end

