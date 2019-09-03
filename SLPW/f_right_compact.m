function du = compact_f_right(t,input,Param)

if nargin >= 2
    q = input(1:3);
    qdot = input(4:6);
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = 0;

    q1dot = qdot(1);
    q2dot = qdot(2);
    q3dot = qdot(3);
    q4dot = 0;
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


m1 = Param.m1;
m2 = Param.m2;
m3 = Param.m3;

k1 = Param.k1;
k2 = Param.k2;
c1 = Param.c1;
c2 = Param.c2;

alpha = Param.alpha;


[  X1, Y1, X2, Y2, X3, Y3,...
    dX1dq1, dY1dq1, dX1dq2, dY1dq2, dX1dq3, dY1dq3, dX1dq4, dY1dq4,...
    dX2dq1, dY2dq1, dX2dq2, dY2dq2, dX2dq3, dY2dq3, dX2dq4, dY2dq4,...
    dX3dq1, dY3dq1, dX3dq2, dY3dq2, dX3dq3, dY3dq3, dX3dq4, dY3dq4,...
    ddX1dq1dq1, ddY1dq1dq1, ddX1dq1dq2, ddY1dq1dq2, ddX1dq1dq3, ddY1dq1dq3, ddX1dq1dq4, ddY1dq1dq4, ddX1dq2dq2, ddY1dq2dq2, ddX1dq2dq3, ddY1dq2dq3, ddX1dq2dq4, ddY1dq2dq4, ddX1dq3dq3, ddY1dq3dq3, ddX1dq3dq4, ddY1dq3dq4, ddX1dq4dq4, ddY1dq4dq4, ...
    ddX2dq1dq1, ddY2dq1dq1, ddX2dq1dq2, ddY2dq1dq2, ddX2dq1dq3, ddY2dq1dq3, ddX2dq1dq4, ddY2dq1dq4, ddX2dq2dq2, ddY2dq2dq2, ddX2dq2dq3, ddY2dq2dq3, ddX2dq2dq4, ddY2dq2dq4, ddX2dq3dq3, ddY2dq3dq3, ddX2dq3dq4, ddY2dq3dq4, ddX2dq4dq4, ddY2dq4dq4, ...
    ddX3dq1dq1, ddY3dq1dq1, ddX3dq1dq2, ddY3dq1dq2, ddX3dq1dq3, ddY3dq1dq3, ddX3dq1dq4, ddY3dq1dq4, ddX3dq2dq2, ddY3dq2dq2, ddX3dq2dq3, ddY3dq2dq3, ddX3dq2dq4, ddY3dq2dq4, ddX3dq3dq3, ddY3dq3dq3, ddX3dq3dq4, ddY3dq3dq4, ddX3dq4dq4, ddY3dq4dq4 ]...
    = derivativesSS(input,Param,'B');





%% Repeats
ddX1dq2dq1 = ddX1dq1dq2;
ddX1dq3dq1 = ddX1dq1dq3;
ddX1dq4dq1 = ddX1dq1dq4;
ddX1dq3dq2 = ddX1dq2dq3;
ddX1dq4dq2 = ddX1dq2dq4;
ddX1dq4dq3 = ddX1dq3dq4;

ddX2dq2dq1 = ddX2dq1dq2;
ddX2dq3dq1 = ddX2dq1dq3;
ddX2dq4dq1 = ddX2dq1dq4;
ddX2dq3dq2 = ddX2dq2dq3;
ddX2dq4dq2 = ddX2dq2dq4;
ddX2dq4dq3 = ddX2dq3dq4;

ddX3dq2dq1 = ddX3dq1dq2;
ddX3dq3dq1 = ddX3dq1dq3;
ddX3dq4dq1 = ddX3dq1dq4;
ddX3dq3dq2 = ddX3dq2dq3;
ddX3dq4dq2 = ddX3dq2dq4;
ddX3dq4dq3 = ddX3dq3dq4;


ddY1dq2dq1 = ddY1dq1dq2;
ddY1dq3dq1 = ddY1dq1dq3;
ddY1dq4dq1 = ddY1dq1dq4;
ddY1dq3dq2 = ddY1dq2dq3;
ddY1dq4dq2 = ddY1dq2dq4;
ddY1dq4dq3 = ddY1dq3dq4;

ddY2dq2dq1 = ddY2dq1dq2;
ddY2dq3dq1 = ddY2dq1dq3;
ddY2dq4dq1 = ddY2dq1dq4;
ddY2dq3dq2 = ddY2dq2dq3;
ddY2dq4dq2 = ddY2dq2dq4;
ddY2dq4dq3 = ddY2dq3dq4;

ddY3dq2dq1 = ddY3dq1dq2;
ddY3dq3dq1 = ddY3dq1dq3;
ddY3dq4dq1 = ddY3dq1dq4;
ddY3dq3dq2 = ddY3dq2dq3;
ddY3dq4dq2 = ddY3dq2dq4;
ddY3dq4dq3 = ddY3dq3dq4;


mH=m3;

% dX1dq2 = 0;
% dY1dq2 = 0;
% 
% dX1dq4 = 0;
% dY1dq4 = 0;
g = 9.81;

%% Ths - K = 1/2*Th_aa*q1dot^2 + Th ab *q1dot *th3dot etc.

Theta11 = m1*(dX1dq1^2 + dY1dq1^2) + m2*(dX2dq1^2 + dY2dq1^2) + mH*(dX3dq1^2 + dY3dq1^2);
Theta22 = m1*(dX1dq2^2 + dY1dq2^2) + m2*(dX2dq2^2 + dY2dq2^2) + mH*(dX3dq2^2 + dY3dq2^2);
Theta33 = m1*(dX1dq3^2 + dY1dq3^2) + m2*(dX2dq3^2 + dY2dq3^2) + mH*(dX3dq3^2 + dY3dq3^2);
Theta44 = m1*(dX1dq4^2 + dY1dq4^2) + m2*(dX2dq4^2 + dY2dq4^2) + mH*(dX3dq4^2 + dY3dq4^2);

Theta12 = m1*(dX1dq1*dX1dq2 + dY1dq1*dY1dq2) + m2*(dX2dq1*dX2dq2 + dY2dq1*dY2dq2) + mH*(dX3dq1*dX3dq2 + dY3dq1*dY3dq2);
Theta13 = m1*(dX1dq1*dX1dq3 + dY1dq1*dY1dq3) + m2*(dX2dq1*dX2dq3 + dY2dq1*dY2dq3) + mH*(dX3dq1*dX3dq3 + dY3dq1*dY3dq3);
Theta14 = m1*(dX1dq1*dX1dq4 + dY1dq1*dY1dq4) + m2*(dX2dq1*dX2dq4 + dY2dq1*dY2dq4) + mH*(dX3dq1*dX3dq4 + dY3dq1*dY3dq4);

Theta23 = m1*(dX1dq2*dX1dq3 + dY1dq2*dY1dq3) + m2*(dX2dq2*dX2dq3 + dY2dq2*dY2dq3) + mH*(dX3dq2*dX3dq3 + dY3dq2*dY3dq3);
Theta24 = m1*(dX1dq2*dX1dq4 + dY1dq2*dY1dq4) + m2*(dX2dq2*dX2dq4 + dY2dq2*dY2dq4) + mH*(dX3dq2*dX3dq4 + dY3dq2*dY3dq4);

Theta34 = m1*(dX1dq3*dX1dq4 + dY1dq3*dY1dq4) + m2*(dX2dq3*dX2dq4 + dY2dq3*dY2dq4) + mH*(dX3dq3*dX3dq4 + dY3dq3*dY3dq4);

Theta21 = Theta12;
Theta31 = Theta13;
Theta41 = Theta14;
Theta32 = Theta23;
Theta42 = Theta24;
Theta43 = Theta34;

T = 1/2*Theta11*q1dot^2 + 1/2*Theta22*q2dot^2 + 1/2*Theta33*q3dot^2 + ...
    1/2*Theta44*q4dot^2 + Theta12*q1dot*q2dot + Theta13*q1dot*q3dot + ...
    Theta14*q1dot*q4dot + Theta23*q2dot*q3dot + Theta24*q2dot*q4dot + ...
    Theta34*q3dot*q4dot;

% V = g*(m1*(X1*sin(alpha)+Y1*cos(alpha)) + m2*(X2*sin(alpha)+Y2*cos(alpha)) ...
%     + mH*(X3*sin(alpha)+Y3*cos(alpha))) + 1/2*k1*q3^2 + 1/2*k2*q4^2;

Xdot = dX3dq1*q1dot + dX3dq2*q2dot + dX3dq3*q3dot + dX3dq4*q4dot;

% ME = T + V;

dTdq1dot = Theta11*q1dot + Theta12*q2dot + Theta13*q3dot + ...
    Theta14*q3dot;

% Consider using dTh > dThdq1

% dT = 1/2*dTh_aa*q1dot^2 + 1/2*dTh_bb*q2dot^2 + 1/2*dTh_cc*q3dot^2 + ...
%     1/2*dTh_dd*q4dot^2 + dTh_ab*q1dot*q2dot + dTh_ac*q1dot*q3dot + ...
%     dTh_ad*q1dot*q3dot + dTh_bc*q2dot*q3dot + dTh_bd*q2dot*q4dot + ...
%     dTh_cd*q3dot*q4dot;

% V = g*(m1*dY1 + m2*dY2 + mH*dY3) + 1/2*k1*dq3^2 + 1/2*k2*dq4^2;

% Theta differentials

dTheta11dq1 = + m1*(ddX1dq1dq1*dX1dq1 + dX1dq1*ddX1dq1dq1 + ddY1dq1dq1*dY1dq1 + dY1dq1*ddY1dq1dq1)+ m2*(ddX2dq1dq1*dX2dq1 + dX2dq1*ddX2dq1dq1 + ddY2dq1dq1*dY2dq1 + dY2dq1*ddY2dq1dq1)+ m3*(ddX3dq1dq1*dX3dq1 + dX3dq1*ddX3dq1dq1 + ddY3dq1dq1*dY3dq1 + dY3dq1*ddY3dq1dq1);
dTheta12dq1 = + m1*(ddX1dq1dq1*dX1dq2 + dX1dq1*ddX1dq2dq1 + ddY1dq1dq1*dY1dq2 + dY1dq1*ddY1dq2dq1)+ m2*(ddX2dq1dq1*dX2dq2 + dX2dq1*ddX2dq2dq1 + ddY2dq1dq1*dY2dq2 + dY2dq1*ddY2dq2dq1)+ m3*(ddX3dq1dq1*dX3dq2 + dX3dq1*ddX3dq2dq1 + ddY3dq1dq1*dY3dq2 + dY3dq1*ddY3dq2dq1);
dTheta13dq1 = + m1*(ddX1dq1dq1*dX1dq3 + dX1dq1*ddX1dq3dq1 + ddY1dq1dq1*dY1dq3 + dY1dq1*ddY1dq3dq1)+ m2*(ddX2dq1dq1*dX2dq3 + dX2dq1*ddX2dq3dq1 + ddY2dq1dq1*dY2dq3 + dY2dq1*ddY2dq3dq1)+ m3*(ddX3dq1dq1*dX3dq3 + dX3dq1*ddX3dq3dq1 + ddY3dq1dq1*dY3dq3 + dY3dq1*ddY3dq3dq1);
dTheta14dq1 = + m1*(ddX1dq1dq1*dX1dq4 + dX1dq1*ddX1dq4dq1 + ddY1dq1dq1*dY1dq4 + dY1dq1*ddY1dq4dq1)+ m2*(ddX2dq1dq1*dX2dq4 + dX2dq1*ddX2dq4dq1 + ddY2dq1dq1*dY2dq4 + dY2dq1*ddY2dq4dq1)+ m3*(ddX3dq1dq1*dX3dq4 + dX3dq1*ddX3dq4dq1 + ddY3dq1dq1*dY3dq4 + dY3dq1*ddY3dq4dq1);

dTheta22dq1 = + m1*(ddX1dq2dq1*dX1dq2 + dX1dq2*ddX1dq2dq1 + ddY1dq2dq1*dY1dq2 + dY1dq2*ddY1dq2dq1)+ m2*(ddX2dq2dq1*dX2dq2 + dX2dq2*ddX2dq2dq1 + ddY2dq2dq1*dY2dq2 + dY2dq2*ddY2dq2dq1)+ m3*(ddX3dq2dq1*dX3dq2 + dX3dq2*ddX3dq2dq1 + ddY3dq2dq1*dY3dq2 + dY3dq2*ddY3dq2dq1);
dTheta23dq1 = + m1*(ddX1dq2dq1*dX1dq3 + dX1dq2*ddX1dq3dq1 + ddY1dq2dq1*dY1dq3 + dY1dq2*ddY1dq3dq1)+ m2*(ddX2dq2dq1*dX2dq3 + dX2dq2*ddX2dq3dq1 + ddY2dq2dq1*dY2dq3 + dY2dq2*ddY2dq3dq1)+ m3*(ddX3dq2dq1*dX3dq3 + dX3dq2*ddX3dq3dq1 + ddY3dq2dq1*dY3dq3 + dY3dq2*ddY3dq3dq1);
dTheta24dq1 = + m1*(ddX1dq2dq1*dX1dq4 + dX1dq2*ddX1dq4dq1 + ddY1dq2dq1*dY1dq4 + dY1dq2*ddY1dq4dq1)+ m2*(ddX2dq2dq1*dX2dq4 + dX2dq2*ddX2dq4dq1 + ddY2dq2dq1*dY2dq4 + dY2dq2*ddY2dq4dq1)+ m3*(ddX3dq2dq1*dX3dq4 + dX3dq2*ddX3dq4dq1 + ddY3dq2dq1*dY3dq4 + dY3dq2*ddY3dq4dq1);

dTheta33dq1 = + m1*(ddX1dq3dq1*dX1dq3 + dX1dq3*ddX1dq3dq1 + ddY1dq3dq1*dY1dq3 + dY1dq3*ddY1dq3dq1)+ m2*(ddX2dq3dq1*dX2dq3 + dX2dq3*ddX2dq3dq1 + ddY2dq3dq1*dY2dq3 + dY2dq3*ddY2dq3dq1)+ m3*(ddX3dq3dq1*dX3dq3 + dX3dq3*ddX3dq3dq1 + ddY3dq3dq1*dY3dq3 + dY3dq3*ddY3dq3dq1);
dTheta34dq1 = + m1*(ddX1dq3dq1*dX1dq4 + dX1dq3*ddX1dq4dq1 + ddY1dq3dq1*dY1dq4 + dY1dq3*ddY1dq4dq1)+ m2*(ddX2dq3dq1*dX2dq4 + dX2dq3*ddX2dq4dq1 + ddY2dq3dq1*dY2dq4 + dY2dq3*ddY2dq4dq1)+ m3*(ddX3dq3dq1*dX3dq4 + dX3dq3*ddX3dq4dq1 + ddY3dq3dq1*dY3dq4 + dY3dq3*ddY3dq4dq1);

dTheta44dq1 = + m1*(ddX1dq4dq1*dX1dq4 + dX1dq4*ddX1dq4dq1 + ddY1dq4dq1*dY1dq4 + dY1dq4*ddY1dq4dq1)+ m2*(ddX2dq4dq1*dX2dq4 + dX2dq4*ddX2dq4dq1 + ddY2dq4dq1*dY2dq4 + dY2dq4*ddY2dq4dq1)+ m3*(ddX3dq4dq1*dX3dq4 + dX3dq4*ddX3dq4dq1 + ddY3dq4dq1*dY3dq4 + dY3dq4*ddY3dq4dq1);


dTheta11dq2 = + m1*(ddX1dq1dq2*dX1dq1 + dX1dq1*ddX1dq1dq2 + ddY1dq1dq2*dY1dq1 + dY1dq1*ddY1dq1dq2)+ m2*(ddX2dq1dq2*dX2dq1 + dX2dq1*ddX2dq1dq2 + ddY2dq1dq2*dY2dq1 + dY2dq1*ddY2dq1dq2)+ m3*(ddX3dq1dq2*dX3dq1 + dX3dq1*ddX3dq1dq2 + ddY3dq1dq2*dY3dq1 + dY3dq1*ddY3dq1dq2);
dTheta12dq2 = + m1*(ddX1dq1dq2*dX1dq2 + dX1dq1*ddX1dq2dq2 + ddY1dq1dq2*dY1dq2 + dY1dq1*ddY1dq2dq2)+ m2*(ddX2dq1dq2*dX2dq2 + dX2dq1*ddX2dq2dq2 + ddY2dq1dq2*dY2dq2 + dY2dq1*ddY2dq2dq2)+ m3*(ddX3dq1dq2*dX3dq2 + dX3dq1*ddX3dq2dq2 + ddY3dq1dq2*dY3dq2 + dY3dq1*ddY3dq2dq2);
dTheta13dq2 = + m1*(ddX1dq1dq2*dX1dq3 + dX1dq1*ddX1dq3dq2 + ddY1dq1dq2*dY1dq3 + dY1dq1*ddY1dq3dq2)+ m2*(ddX2dq1dq2*dX2dq3 + dX2dq1*ddX2dq3dq2 + ddY2dq1dq2*dY2dq3 + dY2dq1*ddY2dq3dq2)+ m3*(ddX3dq1dq2*dX3dq3 + dX3dq1*ddX3dq3dq2 + ddY3dq1dq2*dY3dq3 + dY3dq1*ddY3dq3dq2);
dTheta14dq2 = + m1*(ddX1dq1dq2*dX1dq4 + dX1dq1*ddX1dq4dq2 + ddY1dq1dq2*dY1dq4 + dY1dq1*ddY1dq4dq2)+ m2*(ddX2dq1dq2*dX2dq4 + dX2dq1*ddX2dq4dq2 + ddY2dq1dq2*dY2dq4 + dY2dq1*ddY2dq4dq2)+ m3*(ddX3dq1dq2*dX3dq4 + dX3dq1*ddX3dq4dq2 + ddY3dq1dq2*dY3dq4 + dY3dq1*ddY3dq4dq2);

dTheta22dq2 = + m1*(ddX1dq2dq2*dX1dq2 + dX1dq2*ddX1dq2dq2 + ddY1dq2dq2*dY1dq2 + dY1dq2*ddY1dq2dq2)+ m2*(ddX2dq2dq2*dX2dq2 + dX2dq2*ddX2dq2dq2 + ddY2dq2dq2*dY2dq2 + dY2dq2*ddY2dq2dq2)+ m3*(ddX3dq2dq2*dX3dq2 + dX3dq2*ddX3dq2dq2 + ddY3dq2dq2*dY3dq2 + dY3dq2*ddY3dq2dq2);
dTheta23dq2 = + m1*(ddX1dq2dq2*dX1dq3 + dX1dq2*ddX1dq3dq2 + ddY1dq2dq2*dY1dq3 + dY1dq2*ddY1dq3dq2)+ m2*(ddX2dq2dq2*dX2dq3 + dX2dq2*ddX2dq3dq2 + ddY2dq2dq2*dY2dq3 + dY2dq2*ddY2dq3dq2)+ m3*(ddX3dq2dq2*dX3dq3 + dX3dq2*ddX3dq3dq2 + ddY3dq2dq2*dY3dq3 + dY3dq2*ddY3dq3dq2);
dTheta24dq2 = + m1*(ddX1dq2dq2*dX1dq4 + dX1dq2*ddX1dq4dq2 + ddY1dq2dq2*dY1dq4 + dY1dq2*ddY1dq4dq2)+ m2*(ddX2dq2dq2*dX2dq4 + dX2dq2*ddX2dq4dq2 + ddY2dq2dq2*dY2dq4 + dY2dq2*ddY2dq4dq2)+ m3*(ddX3dq2dq2*dX3dq4 + dX3dq2*ddX3dq4dq2 + ddY3dq2dq2*dY3dq4 + dY3dq2*ddY3dq4dq2);


dTheta33dq2 = + m1*(ddX1dq3dq2*dX1dq3 + dX1dq3*ddX1dq3dq2 + ddY1dq3dq2*dY1dq3 + dY1dq3*ddY1dq3dq2)+ m2*(ddX2dq3dq2*dX2dq3 + dX2dq3*ddX2dq3dq2 + ddY2dq3dq2*dY2dq3 + dY2dq3*ddY2dq3dq2)+ m3*(ddX3dq3dq2*dX3dq3 + dX3dq3*ddX3dq3dq2 + ddY3dq3dq2*dY3dq3 + dY3dq3*ddY3dq3dq2);
dTheta34dq2 = + m1*(ddX1dq3dq2*dX1dq4 + dX1dq3*ddX1dq4dq2 + ddY1dq3dq2*dY1dq4 + dY1dq3*ddY1dq4dq2)+ m2*(ddX2dq3dq2*dX2dq4 + dX2dq3*ddX2dq4dq2 + ddY2dq3dq2*dY2dq4 + dY2dq3*ddY2dq4dq2)+ m3*(ddX3dq3dq2*dX3dq4 + dX3dq3*ddX3dq4dq2 + ddY3dq3dq2*dY3dq4 + dY3dq3*ddY3dq4dq2);

dTheta44dq2 = + m1*(ddX1dq4dq2*dX1dq4 + dX1dq4*ddX1dq4dq2 + ddY1dq4dq2*dY1dq4 + dY1dq4*ddY1dq4dq2)+ m2*(ddX2dq4dq2*dX2dq4 + dX2dq4*ddX2dq4dq2 + ddY2dq4dq2*dY2dq4 + dY2dq4*ddY2dq4dq2)+ m3*(ddX3dq4dq2*dX3dq4 + dX3dq4*ddX3dq4dq2 + ddY3dq4dq2*dY3dq4 + dY3dq4*ddY3dq4dq2);


dTheta11dq3 = + m1*(ddX1dq1dq3*dX1dq1 + dX1dq1*ddX1dq1dq3 + ddY1dq1dq3*dY1dq1 + dY1dq1*ddY1dq1dq3)+ m2*(ddX2dq1dq3*dX2dq1 + dX2dq1*ddX2dq1dq3 + ddY2dq1dq3*dY2dq1 + dY2dq1*ddY2dq1dq3)+ m3*(ddX3dq1dq3*dX3dq1 + dX3dq1*ddX3dq1dq3 + ddY3dq1dq3*dY3dq1 + dY3dq1*ddY3dq1dq3);
dTheta12dq3 = + m1*(ddX1dq1dq3*dX1dq2 + dX1dq1*ddX1dq2dq3 + ddY1dq1dq3*dY1dq2 + dY1dq1*ddY1dq2dq3)+ m2*(ddX2dq1dq3*dX2dq2 + dX2dq1*ddX2dq2dq3 + ddY2dq1dq3*dY2dq2 + dY2dq1*ddY2dq2dq3)+ m3*(ddX3dq1dq3*dX3dq2 + dX3dq1*ddX3dq2dq3 + ddY3dq1dq3*dY3dq2 + dY3dq1*ddY3dq2dq3);
dTheta13dq3 = + m1*(ddX1dq1dq3*dX1dq3 + dX1dq1*ddX1dq3dq3 + ddY1dq1dq3*dY1dq3 + dY1dq1*ddY1dq3dq3)+ m2*(ddX2dq1dq3*dX2dq3 + dX2dq1*ddX2dq3dq3 + ddY2dq1dq3*dY2dq3 + dY2dq1*ddY2dq3dq3)+ m3*(ddX3dq1dq3*dX3dq3 + dX3dq1*ddX3dq3dq3 + ddY3dq1dq3*dY3dq3 + dY3dq1*ddY3dq3dq3);
dTheta14dq3 = + m1*(ddX1dq1dq3*dX1dq4 + dX1dq1*ddX1dq4dq3 + ddY1dq1dq3*dY1dq4 + dY1dq1*ddY1dq4dq3)+ m2*(ddX2dq1dq3*dX2dq4 + dX2dq1*ddX2dq4dq3 + ddY2dq1dq3*dY2dq4 + dY2dq1*ddY2dq4dq3)+ m3*(ddX3dq1dq3*dX3dq4 + dX3dq1*ddX3dq4dq3 + ddY3dq1dq3*dY3dq4 + dY3dq1*ddY3dq4dq3);

dTheta22dq3 = + m1*(ddX1dq2dq3*dX1dq2 + dX1dq2*ddX1dq2dq3 + ddY1dq2dq3*dY1dq2 + dY1dq2*ddY1dq2dq3)+ m2*(ddX2dq2dq3*dX2dq2 + dX2dq2*ddX2dq2dq3 + ddY2dq2dq3*dY2dq2 + dY2dq2*ddY2dq2dq3)+ m3*(ddX3dq2dq3*dX3dq2 + dX3dq2*ddX3dq2dq3 + ddY3dq2dq3*dY3dq2 + dY3dq2*ddY3dq2dq3);
dTheta23dq3 = + m1*(ddX1dq2dq3*dX1dq3 + dX1dq2*ddX1dq3dq3 + ddY1dq2dq3*dY1dq3 + dY1dq2*ddY1dq3dq3)+ m2*(ddX2dq2dq3*dX2dq3 + dX2dq2*ddX2dq3dq3 + ddY2dq2dq3*dY2dq3 + dY2dq2*ddY2dq3dq3)+ m3*(ddX3dq2dq3*dX3dq3 + dX3dq2*ddX3dq3dq3 + ddY3dq2dq3*dY3dq3 + dY3dq2*ddY3dq3dq3);
dTheta24dq3 = + m1*(ddX1dq2dq3*dX1dq4 + dX1dq2*ddX1dq4dq3 + ddY1dq2dq3*dY1dq4 + dY1dq2*ddY1dq4dq3)+ m2*(ddX2dq2dq3*dX2dq4 + dX2dq2*ddX2dq4dq3 + ddY2dq2dq3*dY2dq4 + dY2dq2*ddY2dq4dq3)+ m3*(ddX3dq2dq3*dX3dq4 + dX3dq2*ddX3dq4dq3 + ddY3dq2dq3*dY3dq4 + dY3dq2*ddY3dq4dq3);

dTheta33dq3 = + m1*(ddX1dq3dq3*dX1dq3 + dX1dq3*ddX1dq3dq3 + ddY1dq3dq3*dY1dq3 + dY1dq3*ddY1dq3dq3)+ m2*(ddX2dq3dq3*dX2dq3 + dX2dq3*ddX2dq3dq3 + ddY2dq3dq3*dY2dq3 + dY2dq3*ddY2dq3dq3)+ m3*(ddX3dq3dq3*dX3dq3 + dX3dq3*ddX3dq3dq3 + ddY3dq3dq3*dY3dq3 + dY3dq3*ddY3dq3dq3);
dTheta34dq3 = + m1*(ddX1dq3dq3*dX1dq4 + dX1dq3*ddX1dq4dq3 + ddY1dq3dq3*dY1dq4 + dY1dq3*ddY1dq4dq3)+ m2*(ddX2dq3dq3*dX2dq4 + dX2dq3*ddX2dq4dq3 + ddY2dq3dq3*dY2dq4 + dY2dq3*ddY2dq4dq3)+ m3*(ddX3dq3dq3*dX3dq4 + dX3dq3*ddX3dq4dq3 + ddY3dq3dq3*dY3dq4 + dY3dq3*ddY3dq4dq3);

dTheta44dq3 = + m1*(ddX1dq4dq3*dX1dq4 + dX1dq4*ddX1dq4dq3 + ddY1dq4dq3*dY1dq4 + dY1dq4*ddY1dq4dq3)+ m2*(ddX2dq4dq3*dX2dq4 + dX2dq4*ddX2dq4dq3 + ddY2dq4dq3*dY2dq4 + dY2dq4*ddY2dq4dq3)+ m3*(ddX3dq4dq3*dX3dq4 + dX3dq4*ddX3dq4dq3 + ddY3dq4dq3*dY3dq4 + dY3dq4*ddY3dq4dq3);


dTheta11dq4 = + m1*(ddX1dq1dq4*dX1dq1 + dX1dq1*ddX1dq1dq4 + ddY1dq1dq4*dY1dq1 + dY1dq1*ddY1dq1dq4)+ m2*(ddX2dq1dq4*dX2dq1 + dX2dq1*ddX2dq1dq4 + ddY2dq1dq4*dY2dq1 + dY2dq1*ddY2dq1dq4)+ m3*(ddX3dq1dq4*dX3dq1 + dX3dq1*ddX3dq1dq4 + ddY3dq1dq4*dY3dq1 + dY3dq1*ddY3dq1dq4);
dTheta12dq4 = + m1*(ddX1dq1dq4*dX1dq2 + dX1dq1*ddX1dq2dq4 + ddY1dq1dq4*dY1dq2 + dY1dq1*ddY1dq2dq4)+ m2*(ddX2dq1dq4*dX2dq2 + dX2dq1*ddX2dq2dq4 + ddY2dq1dq4*dY2dq2 + dY2dq1*ddY2dq2dq4)+ m3*(ddX3dq1dq4*dX3dq2 + dX3dq1*ddX3dq2dq4 + ddY3dq1dq4*dY3dq2 + dY3dq1*ddY3dq2dq4);
dTheta13dq4 = + m1*(ddX1dq1dq4*dX1dq3 + dX1dq1*ddX1dq3dq4 + ddY1dq1dq4*dY1dq3 + dY1dq1*ddY1dq3dq4)+ m2*(ddX2dq1dq4*dX2dq3 + dX2dq1*ddX2dq3dq4 + ddY2dq1dq4*dY2dq3 + dY2dq1*ddY2dq3dq4)+ m3*(ddX3dq1dq4*dX3dq3 + dX3dq1*ddX3dq3dq4 + ddY3dq1dq4*dY3dq3 + dY3dq1*ddY3dq3dq4);
dTheta14dq4 = + m1*(ddX1dq1dq4*dX1dq4 + dX1dq1*ddX1dq4dq4 + ddY1dq1dq4*dY1dq4 + dY1dq1*ddY1dq4dq4)+ m2*(ddX2dq1dq4*dX2dq4 + dX2dq1*ddX2dq4dq4 + ddY2dq1dq4*dY2dq4 + dY2dq1*ddY2dq4dq4)+ m3*(ddX3dq1dq4*dX3dq4 + dX3dq1*ddX3dq4dq4 + ddY3dq1dq4*dY3dq4 + dY3dq1*ddY3dq4dq4);

dTheta22dq4 = + m1*(ddX1dq2dq4*dX1dq2 + dX1dq2*ddX1dq2dq4 + ddY1dq2dq4*dY1dq2 + dY1dq2*ddY1dq2dq4)+ m2*(ddX2dq2dq4*dX2dq2 + dX2dq2*ddX2dq2dq4 + ddY2dq2dq4*dY2dq2 + dY2dq2*ddY2dq2dq4)+ m3*(ddX3dq2dq4*dX3dq2 + dX3dq2*ddX3dq2dq4 + ddY3dq2dq4*dY3dq2 + dY3dq2*ddY3dq2dq4);
dTheta23dq4 = + m1*(ddX1dq2dq4*dX1dq3 + dX1dq2*ddX1dq3dq4 + ddY1dq2dq4*dY1dq3 + dY1dq2*ddY1dq3dq4)+ m2*(ddX2dq2dq4*dX2dq3 + dX2dq2*ddX2dq3dq4 + ddY2dq2dq4*dY2dq3 + dY2dq2*ddY2dq3dq4)+ m3*(ddX3dq2dq4*dX3dq3 + dX3dq2*ddX3dq3dq4 + ddY3dq2dq4*dY3dq3 + dY3dq2*ddY3dq3dq4);
dTheta24dq4 = + m1*(ddX1dq2dq4*dX1dq4 + dX1dq2*ddX1dq4dq4 + ddY1dq2dq4*dY1dq4 + dY1dq2*ddY1dq4dq4)+ m2*(ddX2dq2dq4*dX2dq4 + dX2dq2*ddX2dq4dq4 + ddY2dq2dq4*dY2dq4 + dY2dq2*ddY2dq4dq4)+ m3*(ddX3dq2dq4*dX3dq4 + dX3dq2*ddX3dq4dq4 + ddY3dq2dq4*dY3dq4 + dY3dq2*ddY3dq4dq4);

dTheta33dq4 = + m1*(ddX1dq3dq4*dX1dq3 + dX1dq3*ddX1dq3dq4 + ddY1dq3dq4*dY1dq3 + dY1dq3*ddY1dq3dq4)+ m2*(ddX2dq3dq4*dX2dq3 + dX2dq3*ddX2dq3dq4 + ddY2dq3dq4*dY2dq3 + dY2dq3*ddY2dq3dq4)+ m3*(ddX3dq3dq4*dX3dq3 + dX3dq3*ddX3dq3dq4 + ddY3dq3dq4*dY3dq3 + dY3dq3*ddY3dq3dq4);
dTheta34dq4 = + m1*(ddX1dq3dq4*dX1dq4 + dX1dq3*ddX1dq4dq4 + ddY1dq3dq4*dY1dq4 + dY1dq3*ddY1dq4dq4)+ m2*(ddX2dq3dq4*dX2dq4 + dX2dq3*ddX2dq4dq4 + ddY2dq3dq4*dY2dq4 + dY2dq3*ddY2dq4dq4)+ m3*(ddX3dq3dq4*dX3dq4 + dX3dq3*ddX3dq4dq4 + ddY3dq3dq4*dY3dq4 + dY3dq3*ddY3dq4dq4);

dTheta44dq4 = + m1*(ddX1dq4dq4*dX1dq4 + dX1dq4*ddX1dq4dq4 + ddY1dq4dq4*dY1dq4 + dY1dq4*ddY1dq4dq4)+ m2*(ddX2dq4dq4*dX2dq4 + dX2dq4*ddX2dq4dq4 + ddY2dq4dq4*dY2dq4 + dY2dq4*ddY2dq4dq4)+ m3*(ddX3dq4dq4*dX3dq4 + dX3dq4*ddX3dq4dq4 + ddY3dq4dq4*dY3dq4 + dY3dq4*ddY3dq4dq4);

%% Repeats
dTheta21dq1 = dTheta12dq1;
dTheta31dq1 = dTheta13dq1;
dTheta41dq1 = dTheta14dq1;
dTheta32dq1 = dTheta23dq1;
dTheta42dq1 = dTheta24dq1;
dTheta43dq1 = dTheta34dq1;

dTheta21dq2 = dTheta12dq2;
dTheta31dq2 = dTheta13dq2;
dTheta41dq2 = dTheta14dq2;
dTheta32dq2 = dTheta23dq2;
dTheta42dq2 = dTheta24dq2;
dTheta43dq2 = dTheta34dq2;

dTheta21dq3 = dTheta12dq3;
dTheta31dq3 = dTheta13dq3;
dTheta41dq3 = dTheta14dq3;
dTheta32dq3 = dTheta23dq3;
dTheta42dq3 = dTheta24dq3;
dTheta43dq3 = dTheta34dq3;

dTheta21dq4 = dTheta12dq4;
dTheta31dq4 = dTheta13dq4;
dTheta41dq4 = dTheta14dq4;
dTheta32dq4 = dTheta23dq4;
dTheta42dq4 = dTheta24dq4;
dTheta43dq4 = dTheta34dq4;
%%

dVdq1 = g*(m1*(dX1dq1*sin(alpha)+dY1dq1*cos(alpha)) + ...
    m2*(dX2dq1*sin(alpha)+dY2dq1*cos(alpha)) + ...
    m3*(dX3dq1*sin(alpha)+dY3dq1*cos(alpha)));
dVdq2 = g*(m1*(dX1dq2*sin(alpha)+dY1dq2*cos(alpha)) + ...
    m2*(dX2dq2*sin(alpha)+dY2dq2*cos(alpha)) + ...
    m3*(dX3dq2*sin(alpha)+dY3dq2*cos(alpha)));
dVdq3 = g*(m1*(dX1dq3*sin(alpha)+dY1dq3*cos(alpha)) + ...
    m2*(dX2dq3*sin(alpha)+dY2dq3*cos(alpha)) + ...
    m3*(dX3dq3*sin(alpha)+dY3dq3*cos(alpha))) + k1*q3 + c1*q3dot;
dVdq4 = g*(m1*(dX1dq4*sin(alpha)+dY1dq4*cos(alpha)) + ...
    m2*(dX2dq4*sin(alpha)+dY2dq4*cos(alpha)) + ...
    m3*(dX3dq4*sin(alpha)+dY3dq4*cos(alpha))) + k2*q4 + c2*q4dot;

M = [   Theta11,Theta12,Theta13,Theta14;
        Theta21,Theta22,Theta23,Theta24;
        Theta31,Theta32,Theta33,Theta34;
        Theta41,Theta42,Theta43,Theta44     ];
                                                                                                                                                                                                    %  i,2                                                                                                                                                                                             %   i,3                                                                                                                                                                                            %   i,4                                                                                                                                                                                            %
N = [	 + dTheta11dq1*q1dot - 1/2*dTheta11dq1*q1dot	 + dTheta12dq1*q2dot - 1/2*dTheta12dq1*q2dot	 + dTheta13dq1*q3dot - 1/2*dTheta13dq1*q3dot	 + dTheta14dq1*q4dot - 1/2*dTheta14dq1*q4dot,	 + dTheta11dq2*q1dot - 1/2*dTheta21dq1*q1dot	 + dTheta12dq2*q2dot - 1/2*dTheta22dq1*q2dot	 + dTheta13dq2*q3dot - 1/2*dTheta23dq1*q3dot	 + dTheta14dq2*q4dot - 1/2*dTheta24dq1*q4dot,	 + dTheta11dq3*q1dot - 1/2*dTheta31dq1*q1dot	 + dTheta12dq3*q2dot - 1/2*dTheta32dq1*q2dot	 + dTheta13dq3*q3dot - 1/2*dTheta33dq1*q3dot	 + dTheta14dq3*q4dot - 1/2*dTheta34dq1*q4dot,	 + dTheta11dq4*q1dot - 1/2*dTheta41dq1*q1dot	 + dTheta12dq4*q2dot - 1/2*dTheta42dq1*q2dot	 + dTheta13dq4*q3dot - 1/2*dTheta43dq1*q3dot	 + dTheta14dq4*q4dot - 1/2*dTheta44dq1*q4dot,;
		 + dTheta21dq1*q1dot - 1/2*dTheta11dq2*q1dot	 + dTheta22dq1*q2dot - 1/2*dTheta12dq2*q2dot	 + dTheta23dq1*q3dot - 1/2*dTheta13dq2*q3dot	 + dTheta24dq1*q4dot - 1/2*dTheta14dq2*q4dot,	 + dTheta21dq2*q1dot - 1/2*dTheta21dq2*q1dot	 + dTheta22dq2*q2dot - 1/2*dTheta22dq2*q2dot	 + dTheta23dq2*q3dot - 1/2*dTheta23dq2*q3dot	 + dTheta24dq2*q4dot - 1/2*dTheta24dq2*q4dot,	 + dTheta21dq3*q1dot - 1/2*dTheta31dq2*q1dot	 + dTheta22dq3*q2dot - 1/2*dTheta32dq2*q2dot	 + dTheta23dq3*q3dot - 1/2*dTheta33dq2*q3dot	 + dTheta24dq3*q4dot - 1/2*dTheta34dq2*q4dot,	 + dTheta21dq4*q1dot - 1/2*dTheta41dq2*q1dot	 + dTheta22dq4*q2dot - 1/2*dTheta42dq2*q2dot	 + dTheta23dq4*q3dot - 1/2*dTheta43dq2*q3dot	 + dTheta24dq4*q4dot - 1/2*dTheta44dq2*q4dot,;
		 + dTheta31dq1*q1dot - 1/2*dTheta11dq3*q1dot	 + dTheta32dq1*q2dot - 1/2*dTheta12dq3*q2dot	 + dTheta33dq1*q3dot - 1/2*dTheta13dq3*q3dot	 + dTheta34dq1*q4dot - 1/2*dTheta14dq3*q4dot,	 + dTheta31dq2*q1dot - 1/2*dTheta21dq3*q1dot	 + dTheta32dq2*q2dot - 1/2*dTheta22dq3*q2dot	 + dTheta33dq2*q3dot - 1/2*dTheta23dq3*q3dot	 + dTheta34dq2*q4dot - 1/2*dTheta24dq3*q4dot,	 + dTheta31dq3*q1dot - 1/2*dTheta31dq3*q1dot	 + dTheta32dq3*q2dot - 1/2*dTheta32dq3*q2dot	 + dTheta33dq3*q3dot - 1/2*dTheta33dq3*q3dot	 + dTheta34dq3*q4dot - 1/2*dTheta34dq3*q4dot,	 + dTheta31dq4*q1dot - 1/2*dTheta41dq3*q1dot	 + dTheta32dq4*q2dot - 1/2*dTheta42dq3*q2dot	 + dTheta33dq4*q3dot - 1/2*dTheta43dq3*q3dot	 + dTheta34dq4*q4dot - 1/2*dTheta44dq3*q4dot,;
		 + dTheta41dq1*q1dot - 1/2*dTheta11dq4*q1dot	 + dTheta42dq1*q2dot - 1/2*dTheta12dq4*q2dot	 + dTheta43dq1*q3dot - 1/2*dTheta13dq4*q3dot	 + dTheta44dq1*q4dot - 1/2*dTheta14dq4*q4dot,	 + dTheta41dq2*q1dot - 1/2*dTheta21dq4*q1dot	 + dTheta42dq2*q2dot - 1/2*dTheta22dq4*q2dot	 + dTheta43dq2*q3dot - 1/2*dTheta23dq4*q3dot	 + dTheta44dq2*q4dot - 1/2*dTheta24dq4*q4dot,	 + dTheta41dq3*q1dot - 1/2*dTheta31dq4*q1dot	 + dTheta42dq3*q2dot - 1/2*dTheta32dq4*q2dot	 + dTheta43dq3*q3dot - 1/2*dTheta33dq4*q3dot	 + dTheta44dq3*q4dot - 1/2*dTheta34dq4*q4dot,	 + dTheta41dq4*q1dot - 1/2*dTheta41dq4*q1dot	 + dTheta42dq4*q2dot - 1/2*dTheta42dq4*q2dot	 + dTheta43dq4*q3dot - 1/2*dTheta43dq4*q3dot	 + dTheta44dq4*q4dot - 1/2*dTheta44dq4*q4dot,; ];

G = [ dVdq1;  dVdq2;  dVdq3;  dVdq4 ];  
     
 

% qddottest = M(1:3,1:3)\(-N(1:3,1:3)*qdot(1:3) - G(1:3) );

qddottest = linsolve(M(1:3,1:3),(-N(1:3,1:3)*qdot(1:3)-G(1:3) ));

% qddot = M\(N*qdot + G);

du = [ q1dot;   q2dot;  q3dot; 
        qddottest];

end


