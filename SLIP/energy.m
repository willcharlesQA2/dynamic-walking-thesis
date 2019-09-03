function [ME,KE,PE] = energy(t,u,Param,Cond)

    L0 = Param.L0;
    fr1 = Param.fr1;
    fr2 = Param.fr2;
    m1 = Param.m;
    k1 = Param.k1;
    k2 = Param.k2;
    Aa = Param.Aa;
    g = 9.81;
    
    
    % Preallocate to increase speed
        KE = zeros(size(t));
        PE = zeros(size(t));
for i = 1:size(t)
    % DoFs
    th = u(i,1);
    r = u(i,2);
    thdot = u(i,3);
    rdot = u(i,4);

    
qdot = [thdot;rdot];


[ CoP,thC,dCoP,~,dxc,dyc,~,~ ] = findCoP( th,Param,1 );

% If cond is locked then spring has no stiffness
switch Cond
    case {'SS'}
        % Front Spring Energy
        k2 = 0;
        IC = 0;
    case {'DS'}
        k2 = k2;
        IC = Param.IC;
    otherwise
        error('Case not recognised')
end


L = L0 + r;
dL = 1;

% Local position of CoP
xc =  fr1*sin(thC);
yc = -fr1*cos(thC) + fr1;

X =  (-xc)*cos(th) + (L-yc)*sin(th) + CoP;
Y = -(-xc)*sin(th) + (L-yc)*cos(th);

% Mass 1
% Single derivatives
dX1dth = dCoP + sin(th)*xc - dxc*cos(th) - dyc*sin(th) + cos(th)*(L - yc);
dY1dth = cos(th)*xc - dyc*cos(th) + dxc*sin(th) - sin(th)*(L - yc);
dX1dr = dL*sin(th);
dY1dr = dL*cos(th);


% Theta for mass matrix
Thetathth = + m1*(dX1dth*dX1dth + dY1dth*dY1dth);
Thetathr = + m1*(dX1dth*dX1dr + dY1dth*dY1dr);
Thetarr = + m1*(dX1dr*dX1dr + dY1dr*dY1dr);
Thetarth = Thetathr;


M = [	Thetathth,	Thetathr,;
		Thetarth,	Thetarr ];


KE(i) = 1/2*(qdot'*M*qdot);

%     KE2(i) = 1/2*Theta11*q1dot^2 + 1/2*Theta22*q2dot^2  ...
%          + Theta12*q1dot*q2dot;
    
% Potential energies (no slope) & Spring(ankle) = 1/2*kA*(phi - theta)^2
 [ th2,r2t,th2dot,r2dot ] = newConditions( [th,r,thdot,rdot],Param );
 

[ CoP2,thC2,dCoP2,ddCoP2,dxc2,dyc2,ddxc2,ddyc2 ] = findCoP( th2,Param,2 );
xc2 =  Param.fr2*sin(thC2);
yc2 = -Param.fr2*cos(thC2) + Param.fr2;


 
L2 = (Y - xc2*sin(th2))/cos(th2) + yc2;

r2 = L2 - L0;

PE(i) = m1*g*Y + 1/2*k1*r^2 + 1/2*k2*r2^2;



end

% Mechanical energy
ME = KE + PE;

end