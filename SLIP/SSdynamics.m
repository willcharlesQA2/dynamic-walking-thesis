function [ du ] = SSdynamics( t,u,Param )
% Dynamics for walking with the spring attached to the CoP.
th      = u(1);
r       = u(2);
thdot   = u(3);
rdot    = u(4);

qdot = [thdot;rdot];

L0 = Param.L0;  % Natural leg length
k  = Param.k1;   % Leg Stiffness (N/m)
m1 = Param.m;   % mass (kg)
fr = Param.fr1;  % foot radius
g = 9.81;


[ ~,thC,dCoP,ddCoP,dxc,dyc,ddxc,ddyc ] = findCoP( th,Param,1 );

L = L0 + r;
dL = 1;
ddL = 0;

% Local position of CoP
xc =  fr*sin(thC);
yc = -fr*cos(thC) + fr;


% Mass 1
% Single derivatives
dX1dth = dCoP + sin(th)*xc - dxc*cos(th) - dyc*sin(th) + cos(th)*(L - yc);
dY1dth = cos(th)*xc - dyc*cos(th) + dxc*sin(th) - sin(th)*(L - yc);
dX1dr = dL*sin(th);
dY1dr = dL*cos(th);

% Double Derivatives
ddX1dthdth = ddCoP + cos(th)*xc - 2*dyc*cos(th) + 2*dxc*sin(th) - ddxc*cos(th) - sin(th)*(L - yc) - ddyc*sin(th);
ddY1dthdth = 2*dxc*cos(th) - sin(th)*xc + 2*dyc*sin(th) - cos(th)*(L - yc) - ddyc*cos(th) + ddxc*sin(th);
ddX1dthdr = dL*cos(th);
ddY1dthdr = -dL*sin(th);
ddX1drdr = ddL*sin(th);
ddY1drdr = ddL*cos(th);


%% Derivative repeats
ddX1drdth = ddX1dthdr;
ddY1drdth = ddY1dthdr;

% Theta for mass matrix
Thetathth = + m1*(dX1dth*dX1dth + dY1dth*dY1dth);
Thetathr = + m1*(dX1dth*dX1dr + dY1dth*dY1dr);
Thetarr = + m1*(dX1dr*dX1dr + dY1dr*dY1dr);

dThetaththdth = + m1*(ddX1dthdth*dX1dth + dX1dth*ddX1dthdth + ddY1dthdth*dY1dth + dY1dth*ddY1dthdth);
dThetathrdth = + m1*(ddX1dthdth*dX1dr + dX1dth*ddX1drdth + ddY1dthdth*dY1dr + dY1dth*ddY1drdth);

dThetarrdth = + m1*(ddX1drdth*dX1dr + dX1dr*ddX1drdth + ddY1drdth*dY1dr + dY1dr*ddY1drdth);
dThetaththdr = + m1*(ddX1dthdr*dX1dth + dX1dth*ddX1dthdr + ddY1dthdr*dY1dth + dY1dth*ddY1dthdr);
dThetathrdr = + m1*(ddX1dthdr*dX1dr + dX1dth*ddX1drdr + ddY1dthdr*dY1dr + dY1dth*ddY1drdr);
dThetarrdr = + m1*(ddX1drdr*dX1dr + dX1dr*ddX1drdr + ddY1drdr*dY1dr + dY1dr*ddY1drdr);

Thetarth = Thetathr;

dThetarthdth = dThetathrdth;
dThetarthdr = dThetathrdr;

M = [	Thetathth,	Thetathr,;
		Thetarth,	Thetarr ];

N = [	 + dThetaththdth*thdot - 1/2*dThetaththdth*thdot	 + dThetathrdth*rdot - 1/2*dThetathrdth*rdot,	 + dThetaththdr*thdot - 1/2*dThetarthdth*thdot	 + dThetathrdr*rdot - 1/2*dThetarrdth*rdot,;
		 + dThetarthdth*thdot  - 1/2*dThetaththdr*thdot      + dThetarrdth*rdot  - 1/2*dThetathrdr*rdot,	 + dThetarthdr*thdot  - 1/2*dThetarthdr*thdot	 + dThetarrdr*rdot  - 1/2*dThetarrdr*rdot   ];

G = [ m1*g*dY1dth; m1*g*dY1dr + k*r ];
     
% PE = m1*g*Y1 + 1/2*k*r^2;

U = M\(-N*qdot-G);


du = [ thdot; rdot;
        U(1);U(2)];
    
end

