function [ du ] = DSdynamics( t,u,Param )
% Dynamics for walking with the spring attached to the CoP.
th      = u(1);
r       = u(2);
thdot   = u(3);
rdot    = u(4);

qdot = [thdot;rdot];

L0 = Param.L0;  % Natural leg length
k1  = Param.k1;   % Leg Stiffnes (N/m)
k2  = Param.k2;
m1 = Param.m;   % mass (kg)
fr1 = Param.fr1;  % foot radius of rear leg
fr2 = Param.fr2;  % foot radius of front leg
IC = Param.IC;
g = 9.81;


[ CoP,thC,dCoP,ddCoP,dxc,dyc,ddxc,ddyc ] = findCoP( th,Param,1 );

L = L0 + r;
dL = 1;
ddL = 0;

% Local position of CoP
xc =  fr1*sin(thC);
yc = -fr1*cos(thC) + fr1;

X1 =  (-xc)*cos(th) + (L-yc)*sin(th) + CoP;
Y1 = -(-xc)*sin(th) + (L-yc)*cos(th);

% if X1 > Param.IC
%     du = NaN(4,1);
% %     error('WTF')
%  

%% Find th2 and r2 from th1 and r1
[ th2,r2,th2dot,r2dot ] = newConditions( u',Param );
[ CoP2,thC2,dCoP2,ddCoP2,dxc2,dyc2,ddxc2,ddyc2 ] = findCoP( th2,Param,2 );
xc2 =  Param.fr2*sin(thC2);
yc2 = -Param.fr2*cos(thC2) + Param.fr2;

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

     
% dSPR1 = (k*dY1dth*(yc2 - L0 + (Y1 - xc2*sin(th2))/cos(th2)))/cos(th2);
% dSPR2 = (k*dY1dr* (yc2 - L0 + (Y1 - xc2*sin(th2))/cos(th2)))/cos(th2);  

% when leg is locked 
% 0 = (Y - xc2*sin(th2) - (sin(th2)^2*(xc2*sin(th2) - Y))/cos(th2)^2)*dth2dth1 + (sin(th2)*dYdth1)/cos(th2) - dXdth1

% 0 = (Y1 - xc2*sin(th2) - (sin(th2)^2*(xc2*sin(th2) - Y1))/cos(th2)^2)*dth2dr1 + (sin(th2)*dY1dr)/cos(th2) - dX1dr

dth2dth1 = - ((sin(th2)*dY1dth)/cos(th2) - dX1dth) / (Y1 - xc2*sin(th2) - (sin(th2)^2*(xc2*sin(th2) - Y1))/cos(th2)^2);
dth2dr1 = -( (sin(th2)*dY1dr)/cos(th2) - dX1dr) / (Y1 - xc2*sin(th2) - (sin(th2)^2*(xc2*sin(th2) - Y1))/cos(th2)^2);

%% I think this is where to change to k2: (Also fix ME function)
dSPR1 = -k2*((dY1dth - xc2*cos(th2)*dth2dth1)/cos(th2) - (sin(th2)*(xc2*sin(th2) - Y1)*dth2dth1)/cos(th2)^2)*(L0 - yc2 + (xc2*sin(th2) - Y1)/cos(th2));
dSPR2 = -k2*((dY1dr - xc2*cos(th2)*dth2dr1)/cos(th2) - (sin(th2)*(xc2*sin(th2) - Y1)*dth2dr1)/cos(th2)^2)*(L0 - yc2 + (xc2*sin(th2) - Y1)/cos(th2));

% dSPR1 = k*r2*r2dot/thdot;
% dSPR2 = k*r2*r2dot/rdot;

G = [ m1*g*dY1dth       +  dSPR1 ;
      m1*g*dY1dr + k1*r  +  dSPR2 ];
     
  
% PE = m1*g*Y1 + 1/2*k*r^2;

U = M\(-N*qdot-G);


du = [ thdot; rdot;
        U(1);U(2)];
    
end

