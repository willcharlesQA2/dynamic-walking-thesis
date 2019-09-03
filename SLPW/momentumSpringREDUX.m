function [ Qpost ] = momentumSpringREDUX( q1,q2,q3,q4,q1dot,q2dot,q3dot,q4dot,step,Param )
%momentumSpring Calculates angular momentum after the single support phase
%but before the double support phase
%   Uses 4 equations:   Angular momentum about the rear leg
%                       Kinetic Energy before = Kinetic Energy after
%                       Xcdot = 0
%                       Ycdot = 0
[ Param ] = switchLeg( Param,step,1 );

m1 = Param.m1;
m2 = Param.m2;
m3 = Param.m3;

% Calculate first derivatives for single-support Pre-imapct and
% double-support post impact

% 
[~,~,dxth,dyth,ddxth,ddyth]=xth_yth(q1,1,0,Param);

[ sth, dsth, ddsth ] = arcLength(q1,1,1,Param);
% % 
% %% Roll-over shape for front leg
[xth2,yth2,dxth2,dyth2,ddxth2,ddyth2]=xth_yth(q1+q2,2,0,Param);
% 
% % sth = 0;
[ s2, dsth2, ddsth2 ] = arcLength((q1+q2),1,2,Param);


[  X1, Y1, X2, Y2, X3, Y3,...
    dX1dq1, dY1dq1, dX1dq2, dY1dq2, dX1dq3, dY1dq3, dX1dq4, dY1dq4,...
    dX2dq1, dY2dq1, dX2dq2, dY2dq2, dX2dq3, dY2dq3, dX2dq4, dY2dq4,...
    dX3dq1, dY3dq1, dX3dq2, dY3dq2, dX3dq3, dY3dq3, dX3dq4, dY3dq4]...
    = derivativesSS([q1,q2,q3],Param,'V');


 [  Xc, Yc,...
    dXcdq1, dYcdq1, dXcdq2, dYcdq2, dXcdq3, dYcdq3, dXcdq4, dYcdq4 ]...
    = derivativesPc([q1,q2,q3,q4,q1dot,q2dot,q3dot,q4dot],Param);
%% Find Centre of Mass
P1 = [X1,Y1];
P2 = [X2,Y2];
P3 = [X3,Y3];
m = [m1,m2,m3];
r = [P1,P2,P3];

mr = [m1*X1+m2*X2+m3*X3,m1*Y1+m2*Y2+m3*Y3];
M = sum(m);

% coordinates of CoM
R = 1/M*mr;
Rx = R(1);
Ry = R(2);
P1 = 0;



%% Angular momentum of P1 and P3 about rear foot
Swoma2 = [m1*(dY1dq1*(X1-sth)-dX1dq1*(Y1-0)) + m3*(dY3dq1*(X3-sth)-dX3dq1*(Y3-0)),...
          m1*(dY1dq2*(X1-sth)-dX1dq2*(Y1-0)) + m3*(dY3dq2*(X3-sth)-dX3dq2*(Y3-0)),...
          m1*(dY1dq3*(X1-sth)-dX1dq3*(Y1-0)) + m3*(dY3dq3*(X3-sth)-dX3dq3*(Y3-0)),...
          m1*(dY1dq4*(X1-sth)-dX1dq4*(Y1-0)) + m3*(dY3dq4*(X3-sth)-dX3dq4*(Y3-0))];

 
 XdotM = [ dXcdq1-dsth2, dXcdq2-dsth2, dXcdq3, dXcdq4 ];
 YdotM = [ dYcdq1, dYcdq2, dYcdq3, dYcdq4 ];
 
 Qpre = [q1dot; q2dot; q3dot; q4dot];



MA = [Swoma2; XdotM; YdotM];
MB = [Swoma2; 0, 0, 0, 0; 0, 0, 0, 0];
%% Kinetic Energy7
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
    
    
%     KE_Before = 1/2*Theta11*q1dot^2 + 1/2*Theta22*q2dot^2 + 1/2*Theta33*q3dot^2 + ...
%         1/2*Theta44*q4dot^2 + Theta12*q1dot*q2dot + Theta13*q1dot*q3dot + ...
%         Theta14*q1dot*q4dot + Theta23*q2dot*q3dot + Theta24*q2dot*q4dot + ...
%         Theta34*q3dot*q4dot;
    
    MassMatrix = [  Theta11, Theta12, Theta13, Theta14;
                    Theta12, Theta22, Theta23, Theta24;
                    Theta13, Theta23, Theta33, Theta34;
                    Theta14, Theta24, Theta34, Theta44   ];
                
    KE_Before = 1/2*Qpre'*MassMatrix*Qpre; 

x0 = [1];
options = optimset('Display','off');
[x,fval] = fsolve(@(x)limitcycle(x,Param,MA,MB,MassMatrix,KE_Before,Qpre),x0,options);

Qpost = [MA; 0, 0, 1, 0]\[MB; 0, 0, x, 0]*Qpre;

% % Test
% Qpost2 = [  1,0,0,0;
%             0,1,0,0;
%             0,0,1,0;
%             XdotM]\...
%         [  1,0,0,0;
%             0,1,0,0;
%             0,0,1,0;
%             0,0,0,0]*Qpre;

end

function [Fout,Qpost] = limitcycle(x,Param,MA,MB,MassMatrix,KE_B,Qpre)


Qpost = [MA; 0, 0, 1, 0]\[MB; 0, 0, x, 0]*Qpre;

    KE_A = 1/2*Qpost'*MassMatrix*Qpost;
    % KE_After - KE_Before
   Fout = KE_A - KE_B;

end
