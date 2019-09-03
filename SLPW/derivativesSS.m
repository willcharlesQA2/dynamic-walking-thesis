function [  X1, Y1, X2, Y2, X3, Y3,...
    dX1dq1, dY1dq1, dX1dq2, dY1dq2, dX1dq3, dY1dq3, dX1dq4, dY1dq4,...
    dX2dq1, dY2dq1, dX2dq2, dY2dq2, dX2dq3, dY2dq3, dX2dq4, dY2dq4,...
    dX3dq1, dY3dq1, dX3dq2, dY3dq2, dX3dq3, dY3dq3, dX3dq4, dY3dq4,...
    ddX1dq1dq1, ddY1dq1dq1, ddX1dq1dq2, ddY1dq1dq2, ddX1dq1dq3, ddY1dq1dq3, ddX1dq1dq4, ddY1dq1dq4, ddX1dq2dq2, ddY1dq2dq2, ddX1dq2dq3, ddY1dq2dq3, ddX1dq2dq4, ddY1dq2dq4, ddX1dq3dq3, ddY1dq3dq3, ddX1dq3dq4, ddY1dq3dq4, ddX1dq4dq4, ddY1dq4dq4, ...
    ddX2dq1dq1, ddY2dq1dq1, ddX2dq1dq2, ddY2dq1dq2, ddX2dq1dq3, ddY2dq1dq3, ddX2dq1dq4, ddY2dq1dq4, ddX2dq2dq2, ddY2dq2dq2, ddX2dq2dq3, ddY2dq2dq3, ddX2dq2dq4, ddY2dq2dq4, ddX2dq3dq3, ddY2dq3dq3, ddX2dq3dq4, ddY2dq3dq4, ddX2dq4dq4, ddY2dq4dq4, ...
    ddX3dq1dq1, ddY3dq1dq1, ddX3dq1dq2, ddY3dq1dq2, ddX3dq1dq3, ddY3dq1dq3, ddX3dq1dq4, ddY3dq1dq4, ddX3dq2dq2, ddY3dq2dq2, ddX3dq2dq3, ddY3dq2dq3, ddX3dq2dq4, ddY3dq2dq4, ddX3dq3dq3, ddY3dq3dq3, ddX3dq3dq4, ddY3dq3dq4, ddX3dq4dq4, ddY3dq4dq4 ]...
    = derivativesSS(input,Param,output)
% Takes inputs for postions and velocities and calculates the derivatives
% and double derivatives for each mass.
% Inputs:   input   - q1,q2,q3,q4,q1dot,q2dot,q3dot,q4dot
%           Param   - model parameters
%           output  - 'V' = single derivatives (Velocities)
%                   - 'A' = double derivatives (Accelerations)
%                   - 'B' = both derivatives (Velocities and Accelerations)

if nargin >= 2

% QDOTS NOT NEEDED FOR THIS. IF THEY ARE, INCLUDE THIS:
%     q1 = input(1);
%     q2 = input(2);
%     q3 = input(3);
% 
% if size(input,1) == 8
%     q4 = input(4);
%     q1dot = input(5);
%     q2dot = input(6);
%     q3dot = input(7);
%     q4dot = input(8);
% else
%     q4 = 0;
%     q1dot = input(4);
%     q2dot = input(5);
%     q3dot = input(6);
%     q4dot = 0;
% end    
    q = input(1:3);
%     qdot = input(4:6);
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = 0;
    
%     q1dot = qdot(1);
%     q2dot = qdot(2);
%     q3dot = qdot(3);
%     q4dot = 0;
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

Lr1 = Param.Lr1;
LH1 = Param.LH1;
xd2 = 0;
yd2 = Param.Lr2-Param.LH1;


S1th = Param.S1th;
S4th = Param.S4th;

% FP2 = 2;
% xq2=0.4;
% yq2 = 0.003;
% sq2=0.02;

% Roll-over shape
[xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(q1,1,0,Param);

% sth = 0;

[ sth, dsth, ddsth ] = arcLength(q1,1,1,Param);


% springs
yr1 = Lr1 + q3;
yH1 = LH1 + q3;
dyr1 = 1;
dyH1 = 1;
ddyH1 = 0;
ddyr1 = 0;



if strcmp(output,'V') == 1 || strcmp(output,'B') == 1
    X3 = -xth*cos(q1) + (yH1 - yth)*sin(q1) + sth;
    Y3 = xth*sin(q1) + (yH1 - yth)*cos(q1);
    
    X1 = -xth*cos(q1) + (yr1 - yth)*sin(q1) + sth;
    Y1 = xth*sin(q1) + (yr1 - yth)*cos(q1);
    
    xp2 =    (xd2)*cos(q2)+(yd2)*sin(q2);
    yp2 =    -(xd2)*sin(q2)+(yd2)*cos(q2)+yH1;
    
    X2 = (xp2-xth)*cos(q1) + (yp2 - yth)*sin(q1) + sth;
    Y2=-(xp2-xth)*sin(q1)+(yp2-yth)*cos(q1);
    
    %% Single Derivatives
    % dq1
    dX3dq1 = -dxth*cos(q1)+xth*sin(q1)-dyth*sin(q1)+(yH1-yth)*cos(q1)+dsth;
    dY3dq1 = dxth*sin(q1)+xth*cos(q1)-dyth*cos(q1)-(yH1-yth)*sin(q1);
    
    dX1dq1 = -dxth*cos(q1)+xth*sin(q1)-dyth*sin(q1)+(yr1-yth)*cos(q1)+dsth;
    dY1dq1 = dxth*sin(q1)+xth*cos(q1)-dyth*cos(q1)-(yr1-yth)*sin(q1);
    
    % dq2
    dX3dq2 = 0;
    dY3dq2 = 0;
    
    dX1dq2 = 0;
    dY1dq2 = 0;
    
    % dq3
    dX3dq3 = dyH1*sin(q1);
    dY3dq3 = dyH1*cos(q1);
    
    dX1dq3 = dyr1*sin(q1);
    dY1dq3 = dyr1*cos(q1);
    
    % dq4
    dX3dq4 = 0;
    dY3dq4 = 0;
    
    dX1dq4 = 0;
    dY1dq4 = 0;
    
    %% X2
    %% Single Derivatives
    dxp2 = -xd2*sin(q2) + yd2*cos(q2);
    dyp2dth2 = -xd2*cos(q2) - yd2*sin(q2);
    dyp2dr1 = dyH1;
    
    dX2dq1 = -dxth*cos(q1) - (xp2-xth)*sin(q1) - dyth*sin(q1) + (yp2-yth)*cos(q1) + dsth;
    dY2dq1 = +dxth*sin(q1) - (xp2-xth)*cos(q1) - dyth*cos(q1) - (yp2-yth)*sin(q1);
    
    dX2dq2 = dxp2*cos(q1) + dyp2dth2*sin(q1);
    dY2dq2 = -dxp2*sin(q1) + dyp2dth2*cos(q1);
    
    dX2dq3 = dyp2dr1*sin(q1);
    dY2dq3 = dyp2dr1*cos(q1);
    
    dX2dq4 = 0;
    dY2dq4 = 0;
    
end    
    %% Double Derivatives
if strcmp(output,'A') == 1 || strcmp(output,'B') == 1
    % dq1dq1
    ddX3dq1dq1 = -ddxth*cos(q1) + 2*dxth*sin(q1) + xth*cos(q1) ...
        - ddyth*sin(q1) - 2*dyth*cos(q1) - (yH1-yth)*sin(q1) + ddsth;
    ddY3dq1dq1 = ddxth*sin(q1) + 2*dxth*cos(q1) - xth*sin(q1) ...
        -ddyth*cos(q1) + 2*dyth*sin(q1) -(yH1-yth)*cos(q1);
    
    ddX1dq1dq1 = -ddxth*cos(q1) + 2*dxth*sin(q1) + xth*cos(q1) ...
        - ddyth*sin(q1) - 2*dyth*cos(q1) - (yr1-yth)*sin(q1) + ddsth;
    ddY1dq1dq1 = ddxth*sin(q1) + 2*dxth*cos(q1) - xth*sin(q1) ...
        -ddyth*cos(q1) + 2*dyth*sin(q1) -(yr1-yth)*cos(q1);
    
    % dq1dq2
    ddX3dq1dq2 = 0;
    ddY3dq1dq2 = 0;
    
    ddX1dq1dq2 = 0;
    ddY1dq1dq2 = 0;
    
    % dq1dq3
    ddX3dq1dq3 = dyH1*cos(q1);         % CHECK THESE
    ddY3dq1dq3 = -dyH1*sin(q1);
    
    ddX1dq1dq3 = dyr1*cos(q1);         % CHECK THESE
    ddY1dq1dq3 = -dyr1*sin(q1);
    
    % dq1dq4
    ddX3dq1dq4 = 0;
    ddY3dq1dq4 = 0;
    
    ddX1dq1dq4 = 0;
    ddY1dq1dq4 = 0;
    
    % dq2dq2
    ddX3dq2dq2 = 0;
    ddY3dq2dq2 = 0;
    
    ddX1dq2dq2 = 0;
    ddY1dq2dq2 = 0;
    
    % dq2dq3
    ddX3dq2dq3 = 0;
    ddY3dq2dq3 = 0;
    
    ddX1dq2dq3 = 0;
    ddY1dq2dq3 = 0;
    
    % dq2dq4
    ddX3dq2dq4 = 0;
    ddY3dq2dq4 = 0;
    
    ddX1dq2dq4 = 0;
    ddY1dq2dq4 = 0;
    
    % dq3dq3
    ddX3dq3dq3 = ddyH1*sin(q1);
    ddY3dq3dq3 = ddyH1*cos(q1);
    
    ddX1dq3dq3 = ddyr1*sin(q1);
    ddY1dq3dq3 = ddyr1*cos(q1);
    
    % dq3dq4
    ddX3dq3dq4 = 0;
    ddY3dq3dq4 = 0;
    
    ddX1dq3dq4 = 0;
    ddY1dq3dq4 = 0;
    
    % dq4dq4
    ddX3dq4dq4 = 0;
    ddY3dq4dq4 = 0;
    
    ddX1dq4dq4 = 0;
    ddY1dq4dq4 = 0;
    
    
    
    %% X2
    ddxp2 = -xd2*cos(q2) - yd2*sin(q2);
    ddyp2ddth2 = xd2*sin(q2) - yd2*cos(q2);
    ddyp2ddr1 = ddyH1;
    ddyp2dth2dr1 = 0;
    
    ddX2dq1dq1 = -ddxth*cos(q1) + 2*dxth*sin(q1) - (xp2-xth)*cos(q1) - ddyth*sin(q1) - 2*dyth*cos(q1) - (yp2-yth)*sin(q1) + ddsth;
    ddY2dq1dq1 = +ddxth*sin(q1) + 2*dxth*cos(q1) + (xp2-xth)*sin(q1) - ddyth*cos(q1) + 2*dyth*sin(q1) - (yp2-yth)*cos(q1);
    
    ddX2dq2dq2 = ddxp2*cos(q1) + ddyp2ddth2*sin(q1);
    ddY2dq2dq2 = -ddxp2*sin(q1) + ddyp2ddth2*cos(q1);
    
    ddX2dq3dq3 = ddyp2ddr1*sin(q1);
    ddY2dq3dq3 = ddyp2ddr1*cos(q1);
    
    ddX2dq1dq2 = -dxp2*sin(q1) + dyp2dth2*cos(q1);
    ddY2dq1dq2 = -dxp2*cos(q1) - dyp2dth2*sin(q1);
    
    ddX2dq1dq3 = dyp2dr1*cos(q1);
    ddY2dq1dq3 = -dyp2dr1*sin(q1);
    
    ddX2dq2dq3 = ddyp2dth2dr1*sin(q1);
    ddY2dq2dq3 = ddyp2dth2dr1*cos(q1);
    
    ddX2dq1dq4 = 0;
    ddY2dq1dq4 = 0;
    
    ddX2dq2dq4 = 0;
    ddY2dq2dq4 = 0;
    
    ddX2dq3dq4 = 0;
    ddY2dq3dq4 = 0;
    
    ddX2dq4dq4 = 0;
    ddY2dq4dq4 = 0;
end
% else
%     warning('Ouput string not recognised in derivative function. Please use ''A'',''B'' or ''V''')

%% If you need just accelerations, copy the following code
% [  X1, Y1, X2, Y2, X3, Y3,...
%     ~,  ~,  ~,  ~,  ~,  ~,  ~,  ~,  ...
%     ~,  ~,  ~,  ~,  ~,  ~,  ~,  ~,  ...
%     ~,  ~,  ~,  ~,  ~,  ~,  ~,  ~,  ...
%     ddX1dq1dq1, ddY1dq1dq1, ddX1dq1dq2, ddY1dq1dq2, ddX1dq1dq3, ddY1dq1dq3, ddX1dq1dq4, ddY1dq1dq4, ddX1dq2dq2, ddY1dq2dq2, ddX1dq2dq3, ddY1dq2dq3, ddX1dq2dq4, ddY1dq2dq4, ddX1dq3dq3, ddY1dq3dq3, ddX1dq3dq4, ddY1dq3dq4, ddX1dq4dq4, ddY1dq4dq4, ...
%     ddX2dq1dq1, ddY2dq1dq1, ddX2dq1dq2, ddY2dq1dq2, ddX2dq1dq3, ddY2dq1dq3, ddX2dq1dq4, ddY2dq1dq4, ddX2dq2dq2, ddY2dq2dq2, ddX2dq2dq3, ddY2dq2dq3, ddX2dq2dq4, ddY2dq2dq4, ddX2dq3dq3, ddY2dq3dq3, ddX2dq3dq4, ddY2dq3dq4, ddX2dq4dq4, ddY2dq4dq4, ...
%     ddX3dq1dq1, ddY3dq1dq1, ddX3dq1dq2, ddY3dq1dq2, ddX3dq1dq3, ddY3dq1dq3, ddX3dq1dq4, ddY3dq1dq4, ddX3dq2dq2, ddY3dq2dq2, ddX3dq2dq3, ddY3dq2dq3, ddX3dq2dq4, ddY3dq2dq4, ddX3dq3dq3, ddY3dq3dq3, ddX3dq3dq4, ddY3dq3dq4, ddX3dq4dq4, ddY3dq4dq4 ]...
%     = derivativesSS(input,Param,'A');
end


