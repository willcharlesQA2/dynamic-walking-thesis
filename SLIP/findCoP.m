function [ CoP,thC,dCoP,ddCoP,dxc,dyc,ddxc,ddyc ] = findCoP( th,Param,leg )
%Calculate CoP (in [X,Y] global coordinates) and derivatives
Th = Param.Th;
Tt = Param.Tt;
if leg == 1
    fr = Param.fr1;
elseif leg == 2
    fr = Param.fr2;
else
    error('Please select 1 or 2')
end

% if theta is less than theta of heel contact
if th <= Th
    CoP = fr*Th;
    dCoP = 0;
    ddCoP = 0;
    thC = Th;
    
    dxc = 0;
    dyc = 0;

    ddxc = 0;
    ddyc = 0;
% if theta is more than theta of toe contact
elseif th >= Tt
    CoP = fr*Tt;
    dCoP = 0;
    ddCoP = 0;
    thC = Tt;
    
    dxc = 0;
    dyc = 0;

    ddxc = 0;
    ddyc =  0;
else
    CoP = fr*th;
    dCoP = fr;
    ddCoP = 0;
    thC = th;
    
    dxc = fr*cos(thC);
    dyc = fr*sin(thC);

    ddxc = -fr*sin(thC);
    ddyc =  fr*cos(thC);
%     error('wtf? no heel or toe')
end

end

