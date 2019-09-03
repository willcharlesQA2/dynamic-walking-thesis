function [ FX,FY ] = GRFLOA2( t1,u1,t2,u2,Param )
%Using the 'line of action' calculate the GRF in terms of X and Y
% Single support and double support

% th      = u(1);
% r       = u(2);
% thdot   = u(3);
% rdot    = u(4);



L0 = Param.L0;  % Natural leg length
k1  = Param.k1;   % Leg Stiffnes (N/m)
k2 = Param.k2;
m1 = Param.m;   % mass (kg)
fr1 = Param.fr1;  % foot radius
fr2 = Param.fr2;
IC = Param.IC;
g = 9.81;

for i = 1:size(t1)
    [ CoPSS(i,1),th1SS(i,1),~,~,~,~,~,~ ] = findCoP( u1(i,1),Param,1 );
end
for i = 1:size(t2)
    [ CoPDS1(i,1),th1DS(i,1),~,~,~,~,~,~ ] = findCoP( u2(i,1),Param,1 );
    [ th2(i,1),r2(i,1),~,~ ] = newConditions( u2(i,:),Param );
    [ CoPDS2(i,1),th2DS(i,1),~,~,~,~,~,~ ] = findCoP( th2(i,1),Param,2 );
end

%loop through answers
rLoop = {u1(:,2),u2(:,2),r2};
thLoop = {u1(:,1),u2(:,1),th2};
tLoop = {t1,t2,t2};
uLoop = {u1,u2,u2};
CoPLoop = {CoPSS,CoPDS1,CoPDS2};
thCLoop = {th1SS,th1DS,th2DS};
ICLoop = {0,0,IC};
kLoop = {k1,k1,k2};
frLoop = {fr1,fr1,fr2};

for loop = 1:3
    
    th = thLoop{loop};
    r = rLoop{loop};
    CoP = CoPLoop{loop};
    thC = thCLoop{loop};
    k = kLoop{loop};
    fr = frLoop{loop};
%     FP = ICLoop{loop};
    
    L = L0 + r;


    % Local position of CoP
    xc =  fr*sin(thC);
    yc = -fr*cos(thC) + fr;

    X =  (-xc).*cos(th) + (L-yc).*sin(th) + CoP;
    Y = -(-xc).*sin(th) + (L-yc).*cos(th);


    % If there are problems, consider if X > IC
    % No need to add Param.IC, as X2 is calculated from th2
    thla1 = atan((X-CoP)./Y);

    Fy = k*r;

    Fx = -Fy.*tan(th-thla1);

    FR = sqrt(Fx.^2+Fy.^2);

    FX{loop} = FR.*sin(thla1);
    FY{loop} = FR.*cos(thla1);
end

% GRFStance = [FX{1}
    

end

