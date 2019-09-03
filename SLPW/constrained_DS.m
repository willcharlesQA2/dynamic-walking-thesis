function [q2,q4,q2dot,q4dot] = constrained_DS(q1,q3,q1dot,q3dot)
% For a constrained double support system, finds q2 & q4 from q1 & q3

warning('Is this still being used? constrained_DS');

if nargin == 0
    q1 = 0.2293;
    q3 = -0.2430;
    q1dot = 2.5575;
    q3dot = 0.1922;
    Xst = 0.84776;
end

global Param

Lr1 = Param.LrA;
Lr2 = Param.LrB;
LH1 = Param.LH1;


S1th = Param.S1th;
S4th = Param.S4th;
S1th2 = Param.S1th2;
S4th2 = Param.S4th2;
h = Param.h;

%% Roll-over shape for rear leg
[xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(q1,1,0,Param);

% sth = 0;

[ sth1, dsth, ddsth ] = arcLength(q1,dxth,dyth,ddxth,ddyth,S1th,S4th,h,0,Param);

% springs
yr1 = Lr1 + q3;
yH1 = LH1 + q3;
dyr1 = 1;
dyH1 = 1;
dyr2 = 1;
dyH2 = 1;
dxH2dr2 = 0;
dyH2dr2 = dyH2;


% X3 = -xth*cos(q1) + (yH1 - yth)*sin(q1) + sth;
Y3 = xth*sin(q1) + (yH1 - yth)*cos(q1);

%------------------------------------------------------------------------------REMEMBER-----%
% REMEMBER!!!! th2 and s2 need to be corrected for golden_DS
% q2 = golden_DS(q1,q3);

%%
q2 = golden_DSYc(q1,q3);

%% Roll-over shape for front leg
[xth2,yth2,dxth2,dyth2,ddxth2,ddyth2]=xth_yth(q1+q2,1,0,Param);

% sth = 0;

% %%CALCULATE ARC LENGTH
% if (q1+q2)<=S1th2          %start point theta(S1)
%     th_arc=S1th2;
% elseif (q1+q2)>=S4th2
%     th_arc=S4th2;      %end point theta(S4)
% else
%     th_arc=(q1+q2);
% end
% 
% % sth = Arc length i.e. s_{\theta}
% if th_arc == S1th2 || th_arc==S4th2
%     dsth2=0; ddsth2=0;
% else
%     dsth2=(dxth2^2+dyth2^2)^0.5;
%     ddsth2=(dxth2*ddxth2+dyth2*ddyth2)/(dxth2^2+dyth2^2)^0.5;
% end
dxth2dth1 = dxth2;
dxth2dth2 = dxth2;
% ddxth2ddth1 = ddxth2;
% ddxth2ddth2 = ddxth2;
% ddxth2dth1dth2 = ddxth2;
dyth2dth1 = dyth2;
dyth2dth2 = dyth2;
% ddyth2ddth1 = ddyth2;
% ddyth2ddth2 = ddyth2;
% ddyth2dth1dth2 = ddyth2;


%% q4
q4 = 0;
% phi2 = q1+q2;
% q4 = (Y3/cos(phi2) - xth2*tan(phi2)) + yth2 - LH1;
% 
% if abs(q4) > 1e-10
%     error('q4 ~= 0')
% end
% 
yH2 = LH1 + q4;
% 
xH2 = 0;

%% qdots

yc = -(xth2-xH2)*sin(q2) + (yth2-yH2)*cos(q2) + yH1;
dycdth1 = -dxth2dth1*sin(q2) + dyth2dth1*cos(q2);
dycdth2 = -dxth2dth2*sin(q2) - (xth2-xH2)*cos(q2) + dyth2dth2*cos(q2) - (yth2-yH2)*sin(q2);
dycdr2 = +dxH2dr2*sin(q2) - dyH2dr2*cos(q2);
dycdr1 = dyH1;

xc =  (xth2-xH2)*cos(q2) + (yth2-yH2)*sin(q2);
dxcdth1 = dxth2dth1*cos(q2) + dyth2dth1*sin(q2);
dxcdth2 = dxth2dth2*cos(q2) - (xth2-xH2)*sin(q2) + dyth2dth2*sin(q2) + (yth2-yH2)*cos(q2);
dxcdr2 = -dxH2dr2*cos(q2) - dyH2dr2*sin(q2);
dxcdr1 = 0;

dYcdq1 = -(dxcdth1-dxth)*sin(q1) - (xc-xth)*cos(q1) + (dycdth1-dyth)*cos(q1) - (yc-yth)*sin(q1);
dYcdq2 = -dxcdth2*sin(q1) + dycdth2*cos(q1);
dYcdq3 = -dxcdr1*sin(q1) + dycdr1*cos(q1);
dYcdq4 = -dxcdr2*sin(q1) + dycdr2*cos(q1);

dXcdq1 = (dxcdth1 - dxth)*cos(q1) - (xc-xth)*sin(q1) + (dycdth1-dyth)*sin(q1) + (yc-yth)*cos(q1) + dsth;
dXcdq2 = dxcdth2*cos(q1) + dycdth2*sin(q1);
dXcdq3 = dxcdr1*cos(q1) + dycdr1*sin(q1);
dXcdq4 = dxcdr2*cos(q1) + dycdr2*sin(q1);

 q2dot = ((dYcdq1*q1dot+dYcdq3*q3dot)/dYcdq4 - (dXcdq1*q1dot+dXcdq3*q3dot)/dXcdq4) / ...
     (dXcdq2/dXcdq4 - dYcdq2/dYcdq4);
 
 q4dot = -(dYcdq1*q1dot + dYcdq2*q2dot + dYcdq3*q3dot)/dYcdq4;