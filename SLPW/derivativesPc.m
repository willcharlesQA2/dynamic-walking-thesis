function [  Xc, Yc,...
    dXcdq1, dYcdq1, dXcdq2, dYcdq2, dXcdq3, dYcdq3, dXcdq4, dYcdq4,...
    ddXcdq1dq1, ddYcdq1dq1, ddXcdq1dq2, ddYcdq1dq2, ddXcdq1dq3, ...
    ddYcdq1dq3, ddXcdq1dq4, ddYcdq1dq4, ddXcdq2dq2, ddYcdq2dq2, ...
    ddXcdq2dq3, ddYcdq2dq3, ddXcdq2dq4, ddYcdq2dq4, ddXcdq3dq3, ...
    ddYcdq3dq3, ddXcdq3dq4, ddYcdq3dq4, ddXcdq4dq4, ddYcdq4dq4 ]...
    = derivativesPc(input,Param)
% Takes inputs for postions and velocities and calculates the derivatives
% and double derivatives for each mass.
% Inputs:   input   - q1,q2,q3,q4,q1dot,q2dot,q3dot,q4dot
%           Param   - model parameters


if nargin >= 2
    q = input(1:4);

    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    
end
% Initial conditions
if nargin == 0
    q1 = -1.2;
    q2 = 1.3;
    q3 = 0.01;
    q4 = 0.04;
    q1dot = -1.2;
    q2dot = 1.3;
    q3dot = 0.01;
    q4dot = 0.04;
    
    qdot = [q1dot;q2dot;q3dot;q4dot];
end
% % Parameters
% global Param

LH1 = Param.LH1;
xd2 = 0;
yd2 = Param.Lr2-Param.LH1;


S1th = Param.S1th;
S4th = Param.S4th;

S1th2 = Param.S1th2;
S4th2 = Param.S4th2;

h = Param.h; 

% FP2 = 2;
% xq2=0.4;
% yq2 = 0.003;
% sq2=0.02;

% Roll-over shape
[xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(q1,1,0,Param);

% sth = 0;

[ sth, dsth, ddsth ] = arcLength(q1,1,1,Param);

%% Roll-over shape for front leg
[xth2,yth2,dxth2,dyth2,ddxth2,ddyth2]=xth_yth(q1+q2,2,0,Param);

[ ~, dsth2, ddsth2 ] = arcLength((q1+q2),1,2,Param);

dxth2dth1 = dxth2;
dxth2dth2 = dxth2;
ddxth2ddth1 = ddxth2;
ddxth2dth1dth2 = ddxth2;
ddxth2ddth2 = ddxth2;

dyth2dth1 = dyth2;
dyth2dth2 = dyth2;
ddyth2ddth1 = ddyth2;
ddyth2dth1dth2 = ddyth2;
ddyth2ddth2 = ddyth2;

% springs
yH1 = LH1 + q3;
dyH1 = 1;
ddyH1 = 0;

%% Swing leg spring
yH2 = LH1 + q4;
dyH2 = 1;
ddyH2 = 0;

%% Swing leg x position
xH2 = 0;

%% Xc, Yc point of collision. Lagrangian multiplier Xc,Yc = 0...
%% xc
dxH2dr2 = 0;
dyH2dr2 = dyH2;
ddxH2ddr2 = 0;
ddyH2ddr2 = ddyH2;


%%
xc =  (xth2-xH2)*cos(q2) + (yth2-yH2)*sin(q2);

dxcdth1 = dxth2dth1*cos(q2) + dyth2dth1*sin(q2);

ddxcddth1 = ddxth2ddth1*cos(q2) + ddyth2ddth1*sin(q2);

dxcdth2 = dxth2dth2*cos(q2) - (xth2-xH2)*sin(q2) + dyth2dth2*sin(q2) + (yth2-yH2)*cos(q2);

ddxcddth2 = ddxth2ddth2*cos(q2) - 2*dxth2dth2*sin(q2) - (xth2-xH2)*cos(q2) + ddyth2ddth2*sin(q2) + 2*dyth2dth2*cos(q2) - (yth2-yH2)*sin(q2);

dxcdr2 = -dxH2dr2*cos(q2) - dyH2dr2*sin(q2);

ddxcddr2 = -ddxH2ddr2*cos(q2) - ddyH2ddr2*sin(q2);

ddxcdth1dth2 = ddxth2dth1dth2*cos(q2) - dxth2dth1*sin(q2) + ddyth2dth1dth2*sin(q2) + dyth2dth1*cos(q2);

ddxcdth2dr2 = +dxH2dr2*sin(q2) - dyH2dr2*cos(q2);       

ddxcdth1dr2 = 0;

dxcdr1 = 0;

ddxcddr1 = 0;

ddxcdth1dr1 = 0;

ddxcdth2dr1 = 0;

ddxcdr1dr2 = 0;

%% yc
yc = -(xth2-xH2)*sin(q2) + (yth2-yH2)*cos(q2) + yH1;

dycdth1 = -dxth2dth1*sin(q2) + dyth2dth1*cos(q2);

ddycddth1 = -ddxth2ddth1*sin(q2) + ddyth2ddth1*cos(q2);

dycdth2 = -dxth2dth2*sin(q2) - (xth2-xH2)*cos(q2) + dyth2dth2*cos(q2) - (yth2-yH2)*sin(q2);

ddycddth2 = -ddxth2ddth2*sin(q2) - 2*dxth2dth2*cos(q2) + (xth2-xH2)*sin(q2) + ddyth2ddth2*cos(q2) - 2*dyth2dth2*sin(q2) - (yth2-yH2)*cos(q2);

dycdr2 = +dxH2dr2*sin(q2) - dyH2dr2*cos(q2);

ddycddr2 = +ddxH2ddr2*sin(q2) - ddyH2ddr2*cos(q2);

ddycdth1dth2 = -ddxth2dth1dth2*sin(q2) - dxth2dth1*cos(q2) + ddyth2dth1dth2*cos(q2) - dyth2dth1*sin(q2);

ddycdth2dr2 = +dxH2dr2*cos(q2) + dyH2dr2*sin(q2);

dycdr1 = dyH1;

ddycddr1 = ddyH1;

ddycdth1dr1 = 0;

ddycdth1dr2 = 0;

ddycdth2dr1 = 0;

ddycdr1dr2 = 0;

%% Xc - You add sth as this is working backwards. sth2 would need to be 
% taken away to find OC (opposite foot 0 position).
Xc = (xc-xth)*cos(q1) + (yc-yth)*sin(q1) + sth;

dXcdq1 = (dxcdth1 - dxth)*cos(q1) - (xc-xth)*sin(q1) + (dycdth1-dyth)*sin(q1) + (yc-yth)*cos(q1) + dsth;

ddXcdq1dq1 = (ddxcddth1 - ddxth)*cos(q1) - 2*(dxcdth1-dxth)*sin(q1) - (xc-xth)*cos(q1) + (ddycddth1 - ddyth)*sin(q1) + 2*(dycdth1-dyth)*cos(q1) - (yc-yth)*sin(q1) + ddsth;

dXcdq2 = dxcdth2*cos(q1) + dycdth2*sin(q1);

ddXcdq2dq2 = ddxcddth2*cos(q1) + ddycddth2*sin(q1);

dXcdq3 = dxcdr1*cos(q1) + dycdr1*sin(q1);

ddXcdq3dq3 = ddxcddr1*cos(q1) + ddycddr1*sin(q1);

dXcdq4 = dxcdr2*cos(q1) + dycdr2*sin(q1);

ddXcdq4dq4 = ddxcddr2*cos(q1) + ddycddr2*sin(q1);

ddXcdq1dq2 = ddxcdth1dth2*cos(q1) - dxcdth2*sin(q1) + ddycdth1dth2*sin(q1) + dycdth2*cos(q1);

ddXcdq1dq3 = ddxcdth1dr1*cos(q1) - dxcdr1*sin(q1) + ddycdth1dr1*sin(q1) + dycdr1*cos(q1);

ddXcdq1dq4 = ddxcdth1dr2*cos(q1) - dxcdr2*sin(q1) + ddycdth1dr2*sin(q1) + dycdr2*cos(q1);

ddXcdq2dq3 = ddxcdth2dr1*cos(q1) + ddycdth2dr1*sin(q1);

ddXcdq2dq4 = ddxcdth2dr2*cos(q1) + ddycdth2dr2*sin(q1);

ddXcdq3dq4 = ddxcdr1dr2*cos(q1) + ddycdr1dr2*sin(q1);

%% Yc
Yc = -(xc-xth)*sin(q1) + (yc-yth)*cos(q1);

dYcdq1 = -(dxcdth1-dxth)*sin(q1) - (xc-xth)*cos(q1) + (dycdth1-dyth)*cos(q1) - (yc-yth)*sin(q1);

ddYcdq1dq1 = -(ddxcddth1 - ddxth)*sin(q1) - 2*(dxcdth1-dxth)*cos(q1) + (xc-xth)*sin(q1) + (ddycddth1 - ddyth)*cos(q1) - 2*(dycdth1-dyth)*sin(q1) - (yc-yth)*cos(q1);

dYcdq2 = -dxcdth2*sin(q1) + dycdth2*cos(q1);

ddYcdq2dq2 = -ddxcddth2*sin(q1) + ddycddth2*cos(q1);

ddYcdq1dq2 = -ddxcdth1dth2*sin(q1) - dxcdth2*cos(q1) + ddycdth1dth2*cos(q1) - dycdth2*sin(q1);

dYcdq3 = -dxcdr1*sin(q1) + dycdr1*cos(q1);

ddYcdq3dq3 = -ddxcddr1*sin(q1) + ddycddr1*cos(q1);

ddYcdq1dq3 = -ddxcdth1dr1*sin(q1) - dxcdr1*cos(q1) + ddycdth1dr1*cos(q1) - dycdr1*sin(q1);

dYcdq4 = -dxcdr2*sin(q1) + dycdr2*cos(q1);

ddYcdq4dq4 = -ddxcddr2*sin(q1) + ddycddr2*cos(q1);

ddYcdq1dq4 = -ddxcdth1dr2*sin(q1) - dxcdr2*cos(q1) + ddycdth1dr2*cos(q1) - dycdr2*sin(q1);

ddYcdq2dq3 = -ddxcdth2dr1*sin(q1) + ddycdth2dr1*cos(q1);

ddYcdq2dq4 = -ddxcdth2dr2*sin(q1) + ddycdth2dr2*cos(q1);

ddYcdq3dq4 = -ddxcdr1dr2*sin(q1) + ddycdr1dr2*cos(q1);

end


