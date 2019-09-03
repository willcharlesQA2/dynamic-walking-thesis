function [ Xddot,Yddot ] = SSaccelerations( t,u,Param )
%Recalculates Xddot and Yddot from t and u.

% Define parameters
L0 = Param.L0;
fr1 = Param.fr1;

%Prelocate to increase speed
du = zeros(size(t,1),4);
Xddot = zeros(size(t));
Yddot = zeros(size(t));

for i = 1:size(t)
    [ du(i,:) ] = SSdynamics( t(i),u(i,:),Param );
    
    th      = u(i,1);
    r       = u(i,2);
    thdot   = u(i,3);
    rdot    = u(i,4);
%     thdot   = du(i,1);
%     rdot    = du(i,2);
    thddot  = du(i,3);
    rddot   = du(i,4);
    
    [ ~,thC,dCoP,ddCoP,dxc,dyc,ddxc,ddyc ] = findCoP( th,Param,1);
    
    L = L0 + r;
    dL = 1;
    ddL = 0;
    
    % Local position of CoP
    xc =  fr1*sin(thC);
    yc = -fr1*cos(thC) + fr1;
    
    
    % Mass 1
    % Single derivatives
    dXdth = dCoP + sin(th)*xc - dxc*cos(th) - dyc*sin(th) + cos(th)*(L - yc);
    dYdth = cos(th)*xc - dyc*cos(th) + dxc*sin(th) - sin(th)*(L - yc);
    dXdr = dL*sin(th);
    dYdr = dL*cos(th);
    
    % Double Derivatives
    ddXdthdth = ddCoP + cos(th)*xc - 2*dyc*cos(th) + 2*dxc*sin(th) - ddxc*cos(th) - sin(th)*(L - yc) - ddyc*sin(th);
    ddYdthdth = 2*dxc*cos(th) - sin(th)*xc + 2*dyc*sin(th) - cos(th)*(L - yc) - ddyc*cos(th) + ddxc*sin(th);
    ddXdthdr = dL*cos(th);
    ddYdthdr = -dL*sin(th);
    ddXdrdr = ddL*sin(th);
    ddYdrdr = ddL*cos(th);
    
    Xddot(i) = dXdth*thddot + dXdr*rddot + ddXdthdth*thdot^2 + 2*ddXdthdr*thdot*rdot + ddXdrdr*rdot^2;
    Yddot(i) = dYdth*thddot + dYdr*rddot + ddYdthdth*thdot^2 + 2*ddYdthdr*thdot*rdot + ddYdrdr*rdot^2;
end


end

