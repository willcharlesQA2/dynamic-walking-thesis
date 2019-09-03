function [ME,KE,PE] = Energy_test(t,u,Param,step)

if nargin < 4
    step = 1;
    warning('Step not defined')
end

% global Param
% 
%    
%     % % Parameters
% 
% m1 = Param.m1;
% m2 = Param.m2;
% m3 = Param.m3;
% 
% k1 = Param.k1;
% k2 = Param.k2;
% alpha = Param.alpha;

g = 9.81;


for i = 1:size(t)
    
    q1 = u(i,1);
    q2 = u(i,2);
    q3 = u(i,3);

if size(u,2) == 8
    symmetry = 0; % if double-support
    q4 = u(i,4);
    q1dot = u(i,5);
    q2dot = u(i,6);
    q3dot = u(i,7);
    q4dot = u(i,8);
else
    symmetry = 1; % if single-support
    q4 = 0;
    q1dot = u(i,4);
    q2dot = u(i,5);
    q3dot = u(i,6);
    q4dot = 0;
end    

[ Param ] = switchLeg( Param,step,symmetry );

% % Parameters

m1 = Param.m1;
m2 = Param.m2;
m3 = Param.m3;

k1 = Param.k1;
k2 = Param.k2;
alpha = Param.alpha;


[  X1, Y1, X2, Y2, X3, Y3,...
    dX1dq1, dY1dq1, dX1dq2, dY1dq2, dX1dq3, dY1dq3, dX1dq4, dY1dq4,...
    dX2dq1, dY2dq1, dX2dq2, dY2dq2, dX2dq3, dY2dq3, dX2dq4, dY2dq4,...
    dX3dq1, dY3dq1, dX3dq2, dY3dq2, dX3dq3, dY3dq3, dX3dq4, dY3dq4]...
    = derivativesSS(u(i,:),Param,'V');
    
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
    
    
    KE(i) = 1/2*Theta11*q1dot^2 + 1/2*Theta22*q2dot^2 + 1/2*Theta33*q3dot^2 + ...
        1/2*Theta44*q4dot^2 + Theta12*q1dot*q2dot + Theta13*q1dot*q3dot + ...
        Theta14*q1dot*q4dot + Theta23*q2dot*q3dot + Theta24*q2dot*q4dot + ...
        Theta34*q3dot*q4dot;
    
    % just single support
    SpringWork(i) = 1/2*k1*q3^2;
    
    PE(i) = g*(m1*(X1*sin(alpha)+Y1*cos(alpha)) + m2*(X2*sin(alpha)+Y2*cos(alpha)) ...
    + m3*(X3*sin(alpha)+Y3*cos(alpha))) + 1/2*k1*q3^2 + 1/2*k2*q4^2;
    
    ME(i) = KE(i) + PE(i);
    
end

% figure(4)
% plot(t,SpringWork)