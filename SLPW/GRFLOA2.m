function [ GRF ] = GRFLOA2( t1,u1,t2,u2,Param,step )
%Uses the 'Line of Action' theory to find the ground reaction forces. This
%function is intended for the 3 masses scenario.
% The Centre of Mass is found firce and then the Line of Action is assumed
% to go through the CoM.
% th1 = u1(:,1);
% r1 = u1(:,3);
% r1dot = u1(:,6);

[ Param ] = switchLeg( Param,step,1 );

m1 = Param.m1;
m2 = Param.m2;
m3 = Param.m3;

% IC = Param.IC;

% alpha = Param.alpha;
k1 = Param.k1;
k2 = Param.k2;
c1 = 0;
c2 = 0;
% 
% th = u(:,1);
% r = u(:,2);
% 
[SX1,SY1,SX2,SY2,SX3,SY3,SXc,SYc,ICnextStep,sthS] = locations(t1,u1,Param,step,'SS');
[DX1,DY1,DX2,DY2,DX3,DY3,DXc,DYc,IC2,sthD] = locations(t2,u2,Param,step,'DS');

for i = 1:size(t2,1)
    [xth2,yth2,dxth2,dyth2,ddxth2,ddyth2]=xth_yth(u2(i,1)+u2(i,2),2,0,Param);
    [ sthD2(i,1), ~,~ ] = arcLength(u2(i,1)+u2(i,2),1,2,Param);
end
% P1 = [SX1,SY1];
% P2 = [SX2,SY2];
% P3 = [SX3,SY3];
m = [m1,m2,m3];
% r = [P1,P2,P3];

mrSS = [m1*SX1+m2*SX2+m3*SX3,m1*SY1+m2*SY2+m3*SY3];
mrDS = [m2*DX1+m1*DX2+m3*DX3,m2*DY1+m1*DY2+m3*DY3];
M = sum(m);

% coordinates of CoM
R_SS = 1/M*mrSS;
R_DS = 1/M*mrDS;

%loop through answers
rLoop = {u1(:,3),u2(:,3),u2(:,4)};
rdotLoop = {u1(:,6),u2(:,7),u2(:,7)};
thLoop = {u1(:,1),u2(:,1),u2(:,1)+u2(:,2)};
tLoop = {t1,t2,t2};
uLoop = {u1,u2,u2};
CoPLoop = {sthS,sthD,sthD2+IC2};
% thCLoop = {th1SS,th1DS,th2DS};
% ICLoop = {0,0,IC};
kLoop = {k1,k2,k1};
cLoop = {c1,c1,c2};
RLoop = {R_SS,R_DS,R_DS};

% F1 = Stance leg   SS
% F2 = Trailing leg DS
% F3 = Lead leg     DS

for loop = 1:3

    r = rLoop{loop};
    r1dot = rdotLoop{loop};
    th = thLoop{loop};   % Might change when leg is locked?
    % tLoop
    % uLoop
    CoP = CoPLoop{loop};
    k = kLoop{loop};
    c = cLoop{loop};
    R = RLoop{loop};


%       % Roll-over shape
% [xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(q1,1,0,Param);
% 
% % sth = 0;
% 
% [ sth, dsth, ddsth ] = arcLength(q1,dxth,dyth,ddxth,ddyth,S1th,S4th,Param.h,1,Param);

 
    % Line of action between CoP and CoM
    thla = atan((R(:,1)-CoP)./R(:,2));

    Fl = k.*r + c.*r1dot;

    Fth = -Fl.*tan(th-thla);

    FX{loop} = (Fl.^2+Fth.^2).^0.5 .*sin(thla);

    FY{loop} = (Fl.^2+Fth.^2).^0.5 .*cos(thla);
% 
end

GRF = struct('FX1', FX{1},'FY1', FY{1}...
            ,'FX2T',FX{2},'FY2T',FY{2}...
            ,'FX2L',FX{3},'FY2L',FY{3}...
            ,'t1'  ,t1   ,'t2'  ,t2   ...
            ,'step',step    );

end

