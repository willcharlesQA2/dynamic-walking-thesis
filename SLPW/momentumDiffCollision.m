function [ TotalMomChange,hmmg ] = momentumDiffCollision( u1End,u2Start,Param )
%Calculates momentum before and after collision
step = 1;
[ Param ] = switchLeg( Param,step,1 );



m1 = Param.m1;
m2 = Param.m2;
m3 = Param.m3;

% bif = 1;

% should be bif,number
% bif
qB      = u1End(:,1:3);
qdotB   = u1End(:,4:6);
% bif + 1
qA      = u2Start(:,1:4);
qdotA   = u2Start(:,5:8);

% Calculate first derivatives for single-support Pre-imapct and
% double-support post impact

q1 = qB(1);
q2 = qB(2);
q3 = qB(3);

% 
[~,~,dxth,dyth,ddxth,ddyth]=xth_yth(q1,1,0,Param);

[ sth, dsth, ddsth ] = arcLength(q1,1,1,Param);
% % 
% %% Roll-over shape for front leg
[xth2,yth2,dxth2,dyth2,ddxth2,ddyth2]=xth_yth(q1+q2,2,0,Param);
% 
% % sth = 0;
[ s2, dsth2, ddsth2 ] = arcLength((q1+q2),1,2,Param);
% 
% [dX1dq1Pre, dY1dq1Pre, dX1dq2Pre, dY1dq2Pre, dX1dq3Pre, dY1dq3Pre, dX1dq4Pre, dY1dq4Pre, ...
%     dX2dq1Pre, dY2dq1Pre, dX2dq2Pre, dY2dq2Pre, dX2dq3Pre, dY2dq3Pre, dX2dq4Pre, dY2dq4Pre, ...
%     dX3dq1Pre, dY3dq1Pre, dX3dq2Pre, dY3dq2Pre, dX3dq3Pre, dY3dq3Pre, dX3dq4Pre, dY3dq4Pre ] = ...
%     velocitiesSS( q1,q2,q3,q4,1);
% 
% [dX1dq1Post, dY1dq1Post, dX1dq2Post, dY1dq2Post, dX1dq3Post, dY1dq3Post, dX1dq4Post, dY1dq4Post, ...
%     dX2dq1Post, dY2dq1Post, dX2dq2Post, dY2dq2Post, dX2dq3Post, dY2dq3Post, dX2dq4Post, dY2dq4Post, ...
%     dX3dq1Post, dY3dq1Post, dX3dq2Post, dY3dq2Post, dX3dq3Post, dY3dq3Post, dX3dq4Post, dY3dq4Post] = ...
%     velocitiesSS( q1,q2,q3,q4,1);

[  X1, Y1, X2, Y2, X3, Y3,...
    dX1dq1, dY1dq1, dX1dq2, dY1dq2, dX1dq3, dY1dq3, dX1dq4, dY1dq4,...
    dX2dq1, dY2dq1, dX2dq2, dY2dq2, dX2dq3, dY2dq3, dX2dq4, dY2dq4,...
    dX3dq1, dY3dq1, dX3dq2, dY3dq2, dX3dq3, dY3dq3, dX3dq4, dY3dq4]...
    = derivativesSS([q1,q2,q3],Param,'V');

% Swing leg foot contact not needed
%  [  Xc, Yc,...
%     dXcdq1, dYcdq1, dXcdq2, dYcdq2, dXcdq3, dYcdq3, dXcdq4, dYcdq4 ]...
%     = derivativesPc([q1,q2,q3,q4,q1dot,q2dot,q3dot,q4dot],Param);

%% Velocities
q1dot = qdotB(1);
q2dot = qdotB(2);
q3dot = qdotB(3);

X1dotB = dX1dq1*q1dot + dX1dq2*q2dot + dX1dq3*q3dot;
Y1dotB = dY1dq1*q1dot + dY1dq2*q2dot + dY1dq3*q3dot;

X2dotB = dX2dq1*q1dot + dX2dq2*q2dot + dX2dq3*q3dot;
Y2dotB = dY2dq1*q1dot + dY2dq2*q2dot + dY2dq3*q3dot;

X3dotB = dX3dq1*q1dot + dX3dq2*q2dot + dX3dq3*q3dot;
Y3dotB = dY3dq1*q1dot + dY3dq2*q2dot + dY3dq3*q3dot;

%% B
q1dot = qdotA(1);
q2dot = qdotA(2);
q3dot = qdotA(3);

X1dotA = dX1dq1*q1dot + dX1dq2*q2dot + dX1dq3*q3dot;
Y1dotA = dY1dq1*q1dot + dY1dq2*q2dot + dY1dq3*q3dot;

X2dotA = dX2dq1*q1dot + dX2dq2*q2dot + dX2dq3*q3dot;
Y2dotA = dY2dq1*q1dot + dY2dq2*q2dot + dY2dq3*q3dot;

X3dotA = dX3dq1*q1dot + dX3dq2*q2dot + dX3dq3*q3dot;
Y3dotA = dY3dq1*q1dot + dY3dq2*q2dot + dY3dq3*q3dot;

%% Find Centre of Mass
P1 = [X1,Y1];
P2 = [X2,Y2];
P3 = [X3,Y3];

P = [X1,Y1;
	 X2,Y2;
	 X3,Y3 ];

m = [m1,m2,m3];
r = [P1,P2,P3];

mr = [m1*X1+m2*X2+m3*X3,m1*Y1+m2*Y2+m3*Y3];
M = sum(m);

% coordinates of CoM
R = 1/M*mr;
Rx = R(1);
Ry = R(2);
P1 = 0;

% Momentum of CoM 
% VB(Mass,X/Y)
VB =    [X1dotB,Y1dotB;
         X2dotB,Y2dotB;
         X3dotB,Y3dotB];

VA =    [X1dotA,Y1dotA;
         X2dotA,Y2dotA;
         X3dotA,Y3dotA];

% MomXB = [   m(1)*VB(1,1);
%             m(2)*VB(2,1);
%             m(3)*VB(3,1)    ];
%         
% MomYB = [   m(1)*VB(1,2);
%             m(2)*VB(2,2);
%             m(3)*VB(3,2)    ];

%%[MomXB,MomYB]
MomB = [m;m]'.*VB;
MomA = [m;m]'.*VA;

%% Difference After - Before
% XDiff = MomXA - MomXB;
% YDiff = MomYA - MomYB;

% [XDiff,YDiff]
MomDiff = MomA - MomB;
SumDiff = sum(MomDiff,1);

%% Angle of direction of momentum difference
%% Angle of Centre of mass
% Rx -sth and Ry - 0
comg = atan(Ry/(Rx-sth));
hmmg = atan(SumDiff(2)/SumDiff(1));

%%Angle of individual CoMs
ng = atan(P(:,2)./(P(:,1)-sth));

% Momentum difference in direction of CoM
% Ch = XDiff*cos(ng) + YDiff*sin(ng);
% Momentum difference tangential to direction of CoM
% In = XDiff*sin(ng) - YDiff*cos(ng);

Ch = MomDiff(:,1).*cos(ng) + MomDiff(:,2).*sin(ng);
In = MomDiff(:,1).*sin(ng) - MomDiff(:,2).*cos(ng);

%% Or total momentum
Ch2 = SumDiff(1).*cos(comg) + SumDiff(2).*sin(comg);
In2 = SumDiff(1).*sin(comg) - SumDiff(2).*cos(comg);

Ch3 = SumDiff(1).*cos(hmmg) + SumDiff(2).*sin(hmmg);
In3 = SumDiff(1).*sin(hmmg) - SumDiff(2).*cos(hmmg);

TotalMomChange = sqrt(SumDiff(1)^2+SumDiff(2)^2);
% 
% SumMomB = sum(MomB,1);
% SumMomA = sum(MomA,1);
% 
% InB = SumMomB(1).*sin(comg) - SumMomB(2).*cos(comg);
% InA = SumMomA(1).*sin(comg) - SumMomA(2).*cos(comg);
%% Check Angular Momentum is conserved
mass = 1;
for mass = 1:3
    angMomB(mass) = m(mass)*(VB(mass,2)*(P(mass,1)-sth) - VB(mass,1)*P(mass,2));
    angMomA(mass) = m(mass)*(VA(mass,2)*(P(mass,1)-sth) - VA(mass,1)*P(mass,2));
end
end
