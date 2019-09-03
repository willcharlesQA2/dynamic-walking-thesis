function [ phi2,r2,phi2dot,r2dot ] = newConditions( u2,Param )
%Finds Initial Conditions for the next step
th      = u2(end,1);
r       = u2(end,2);
thdot   = u2(end,3);
rdot    = u2(end,4);

% For rear foot
[ CoP,thC,dCoP,ddCoP,dxc,dyc,ddxc,ddyc ] = findCoP( th,Param,1 );


xc =  Param.fr1*sin(thC);
yc = -Param.fr1*cos(thC) + Param.fr1;
L = Param.L0 + u2(end,2);
X =  (-xc)*cos(th) + (L-yc)*sin(th) + CoP;
Y = -(-xc)*sin(th) + (L-yc)*cos(th);

if isnan(Y)
    warning('Y is NaN')
end
if X >= Param.IC
%     warning('X is past IC point')
    x0 = [0,pi/2]; % used to be [-pi/2 , pi/6]
else
    x0 = [-pi/2,0]; % used to be [-pi/2 , pi/6]
end
% r2 = L2 - L0;

% figure(10)
% plot([X,Param.IC,Param.IC-X])
% % u2(2)
% drawnow
phi2 = fzero(@(th) fsw(th,X,Y,Param),x0);

% For front foot
[ CoP2,thC2,dCoP2,ddCoP2,dxc2,dyc2,ddxc2,ddyc2 ] = findCoP( phi2,Param,2 );
xc2 =  Param.fr2*sin(thC2);
yc2 = -Param.fr2*cos(thC2) + Param.fr2;

r2 = (Y - xc2*sin(phi2))/cos(phi2) + yc2 - Param.L0;


% r2 = (Y^2 + (Param.IC - X)^2)^0.5 - Param.L0;
% phi2 = -atan((Param.IC-X)/Y);




L = Param.L0 + r;
dL = 1;

L2 = Param.L0 + r2;
dL2 = 1;

% Single derivatives
dX1dth = dCoP + sin(th)*xc - dxc*cos(th) - dyc*sin(th) + cos(th)*(L - yc);
dY1dth = cos(th)*xc - dyc*cos(th) + dxc*sin(th) - sin(th)*(L - yc);
dX1dr = dL*sin(th);
dY1dr = dL*cos(th);

% Single derivatives
dX1dphi2 = dCoP2 + sin(phi2)*xc2 - dxc2*cos(phi2) - dyc2*sin(phi2) + cos(phi2)*(L2 - yc2);
dY1dphi2 = cos(phi2)*xc2 - dyc2*cos(phi2) + dxc2*sin(phi2) - sin(phi2)*(L2 - yc2);
dX1dr2 = dL2*sin(phi2);
dY1dr2 = dL2*cos(phi2);

% So using Xdot- = Xdot+ and Ydot- + Ydot+
% Xdot = L*qdot
Lb = [dX1dth , dX1dr ; dY1dth,  dY1dr ];
La = [dX1dphi2, dX1dr2; dY1dphi2, dY1dr2];
qdot = [thdot;rdot];

[qplusdot] = La\(Lb*qdot);

phi2dot = qplusdot(1);
r2dot  = qplusdot(2);
end


function out = fsw(phi2,X,Y,Param)
[ CoP2,thC2,dCoP2,ddCoP2,dxc2,dyc2,ddxc2,ddyc2 ] = findCoP( phi2,Param,2 );
xc2 =  Param.fr2*sin(thC2);
yc2 = -Param.fr2*cos(thC2) + Param.fr2;

out = CoP2 + Param.IC - xc2*cos(phi2) + (sin(phi2)*(Y - xc2*sin(phi2)))/cos(phi2) - X;
if isnan(out)
    warning('WTF?')
end
end
