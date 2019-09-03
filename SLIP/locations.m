function [X,Y,Xdot,Ydot] = locations(t,u,Param,step)
% Plot walker as a diagram

if nargin == 0
    % Original leg length
    L0 = 1;
    % Foot length
    fl = 0.2;
    fr = 0.3; % foot radius
    theta = -0.5;
    r = 0;
else
    L0 = Param.L0;  % Natural leg length
    fr = Param.fr1;  % foot radius
    Aa = Param.alpha;
end

% fl = L0*thetaEnd;
% Maximum Theta, when leg becomes locked
Th = Param.Th;
Tt = Param.Tt;


theta   = u(:,1);
r       = u(:,2);
thdot   = u(:,3);
rdot    = u(:,4);

L = L0 + r;
dL = 1;

CoP = zeros(size(theta));
dCoP = zeros(size(theta));
dxc = zeros(size(theta));
dyc = zeros(size(theta));

% Locked rocker
for i = 1:size(theta,1);
    [ CoP(i,:),~,dCoP(i,:),~,dxc(i,:),dyc(i,:),~,~ ] = findCoP( theta(i,:),Param,1 );
end

% CoP in local coordinates
xc =  fr*sin(theta);
yc = -fr*cos(theta) + fr;

X =  (-xc).*cos(theta) + (L-yc).*sin(theta) + CoP;
Y = -(-xc).*sin(theta) + (L-yc).*cos(theta);

% Mass 1
% Single derivatives
dXdth = dCoP + sin(theta).*xc - dxc.*cos(theta) - dyc.*sin(theta) + cos(theta).*(L - yc);
dYdth = cos(theta).*xc - dyc.*cos(theta) + dxc.*sin(theta) - sin(theta).*(L - yc);
dXdr = dL*sin(theta);
dYdr = dL*cos(theta);

Xdot = dXdth.*thdot + dXdr.*rdot;
Ydot = dYdth.*thdot + dYdr.*rdot;
end

