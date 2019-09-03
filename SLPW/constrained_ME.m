function [ q2,q1dot,q2dot,q4dot ] = constrained_ME( q1,q3,q3dot,ME_0,Param ) 
% Finds initial conditions from a given th1,r1,r1dot,ME
% find theta_1 first? roots of two values remember
% output - th2,th1dot,th2dot,r2dot

step = 1;

[ Param ] = switchLeg( Param,step,0 );

Param.xd2 = 0;
Param.yd2 = Param.Lr2-Param.LH1;
Param.xd1 = 0;
Param.yd1 = Param.Lr1-Param.LH1;

m1 = Param.m1;
m2 = Param.m2;
m3 = Param.m3;
k1 = Param.k1;
k2 = Param.k2;
alpha = Param.alpha;
% ME_0 = Param.ME_0;

q4 = 0;
q2 = golden_DSYc(q1,q3,Param);

[  X1, Y1, X2, Y2, X3, Y3,...
    dX1dq1, dY1dq1, dX1dq2, dY1dq2, dX1dq3, dY1dq3, dX1dq4, dY1dq4,...
    dX2dq1, dY2dq1, dX2dq2, dY2dq2, dX2dq3, dY2dq3, dX2dq4, dY2dq4,...
    dX3dq1, dY3dq1, dX3dq2, dY3dq2, dX3dq3, dY3dq3, dX3dq4, dY3dq4]...
    = derivativesSS([q1,q2,q3,q4],Param,'V');

[  Xc, Yc,...
    dXcdq1, dYcdq1, dXcdq2, dYcdq2, dXcdq3, dYcdq3, dXcdq4, dYcdq4 ]...
    = derivativesPc([q1,q2,q3,q4],Param);

[xth2,yth2,dxth2,dyth2,ddxth2,ddyth2]=xth_yth(q1+q2,2,0,Param);

[ s2, dsth2, ddsth2 ] = arcLength((q1+q2),1,2,Param);

 XdotM = [ dXcdq1-dsth2, dXcdq2-dsth2, dXcdq3, dXcdq4 ];
 YdotM = [ dYcdq1, dYcdq2, dYcdq3, dYcdq4 ];

 % XdotM*Qpost = 0
 % YdotM*Qpost = 0
 % roots??
 
 g = 9.81;
 
PE = g*(m1*(X1*sin(alpha)+Y1*cos(alpha)) + m2*(X2*sin(alpha)+Y2*cos(alpha)) ...
+ m3*(X3*sin(alpha)+Y3*cos(alpha))) + 1/2*k1*q3^2 + 1/2*k2*q4^2;

Theta11 = m1*(dX1dq1^2 + dY1dq1^2) + m2*(dX2dq1^2 + dY2dq1^2) + m3*(dX3dq1^2 + dY3dq1^2);
Theta22 = m1*(dX1dq2^2 + dY1dq2^2) + m2*(dX2dq2^2 + dY2dq2^2) + m3*(dX3dq2^2 + dY3dq2^2);
Theta33 = m1*(dX1dq3^2 + dY1dq3^2) + m2*(dX2dq3^2 + dY2dq3^2) + m3*(dX3dq3^2 + dY3dq3^2);
Theta44 = m1*(dX1dq4^2 + dY1dq4^2) + m2*(dX2dq4^2 + dY2dq4^2) + m3*(dX3dq4^2 + dY3dq4^2);

Theta12 = m1*(dX1dq1*dX1dq2 + dY1dq1*dY1dq2) + m2*(dX2dq1*dX2dq2 + dY2dq1*dY2dq2) + m3*(dX3dq1*dX3dq2 + dY3dq1*dY3dq2);
Theta13 = m1*(dX1dq1*dX1dq3 + dY1dq1*dY1dq3) + m2*(dX2dq1*dX2dq3 + dY2dq1*dY2dq3) + m3*(dX3dq1*dX3dq3 + dY3dq1*dY3dq3);
Theta14 = m1*(dX1dq1*dX1dq4 + dY1dq1*dY1dq4) + m2*(dX2dq1*dX2dq4 + dY2dq1*dY2dq4) + m3*(dX3dq1*dX3dq4 + dY3dq1*dY3dq4);

Theta23 = m1*(dX1dq2*dX1dq3 + dY1dq2*dY1dq3) + m2*(dX2dq2*dX2dq3 + dY2dq2*dY2dq3) + m3*(dX3dq2*dX3dq3 + dY3dq2*dY3dq3);
Theta24 = m1*(dX1dq2*dX1dq4 + dY1dq2*dY1dq4) + m2*(dX2dq2*dX2dq4 + dY2dq2*dY2dq4) + m3*(dX3dq2*dX3dq4 + dY3dq2*dY3dq4);

Theta34 = m1*(dX1dq3*dX1dq4 + dY1dq3*dY1dq4) + m2*(dX2dq3*dX2dq4 + dY2dq3*dY2dq4) + m3*(dX3dq3*dX3dq4 + dY3dq3*dY3dq4);

MassMatrix = [  Theta11, Theta12, Theta13, Theta14;
                Theta12, Theta22, Theta23, Theta24;
                Theta13, Theta23, Theta33, Theta34;
                Theta14, Theta24, Theta34, Theta44   ];
                
% a = 1/2*Theta11;
% b = Theta12*q2dot + Theta13*q3dot;
% c = 1/2*Theta22*q2dot^2 + 1/2*Theta33*q3dot^2 + Theta23*q2dot*q3dot + PE;

% KE_A = 1/2*Qpost'*MassMatrix*Qpost;

% th1dot,th2dot,r2dot,
x0 = [1.5 -0.4 -0.8];
options = optimset('Display','off');
[x,fval] = fsolve(@(x)limitcycle(x,XdotM,YdotM,MassMatrix,PE,ME_0,q3dot),x0,options);

q1dot = x(1);
q2dot = x(2);
q4dot = x(3);
end

function [Fout] = limitcycle(x,XdotM,YdotM,MassMatrix,PE,ME_0,q3dot)
    
    qdot = [x(1:2),q3dot,x(3)];
    F1 = XdotM*qdot';
    F2 = YdotM*qdot';
    ME_A = 1/2*qdot*MassMatrix*qdot' + PE;
    % ME_After - ME_Before
    F3 = ME_A - ME_0;
    
    Fout = [F1, F2, F3];
end


