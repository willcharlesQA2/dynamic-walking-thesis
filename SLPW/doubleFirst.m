function [ t1,u1,t2,u2,event] = doubleFirst( q,qdot,tend,step,Param )
%oneStep outlines one step of the walker
%   Detailed explanation goes here

% Accuracy of solver
accuracy = 10^(-8);

% Switches around legs A and B
[ Param ] = switchLeg( Param,step,0 );

Param.xd2 = 0;
Param.yd2 = Param.Lr2-Param.LH1;
Param.xd1 = 0;
Param.yd1 = Param.Lr1-Param.LH1;

T2 = 1.0;
T1 = 2.0;

%% Double-support phase
optionsDS=odeset('Events',@(t2,u2)foot_strikeDS(t2,u2,Param),'abstol',accuracy,'reltol',accuracy);
 [t2,u2,~,~,~] = ode45(@(t,y)f_right_LM_compact(t,y,Param), [tend tend+T1], [q,qdot],optionsDS);


% If the take-off event from 'foot_strikeDS' has been met:
if t2(end) == tend+T1
    event = 'noTakeOff';
        t1 = NaN;
        u1 = NaN;
elseif u2(end,1) >= (pi/2) || u2(end,1) <= (-pi/2)
    event = 'fallen';
        t1 = NaN;
        u1 = NaN;
else
%%  Changes stance leg and swing leg.
        q1 = [u2(end,1)+u2(end,2),-u2(end,2),u2(end,4)];
        qdot1 = [u2(end,5)+u2(end,6),-u2(end,6),u2(end,8)];
        
[ Param ] = switchLeg( Param,step,1 );

%%    if not changing: 
%     q = [u2(end,1),u2(end,2),u2(end,3)];
%     qdot = [u2(end,5),u2(end,6),u2(end,7)];

%% CHECK PE GAINED FROM SWITCHING SUPPORT LEG FROM BACK TO FRONT
%     [~,~,~,~,~,~,~,~,IC,sth] = locations(t2(end),u2(end,:),Param);
%     dPE = (Param.m1+Param.m2+Param.m3)*9.81*IC*sin(Param.alpha);
%     
    optionsSS=odeset('Events',@(t1,u1)foot_strikeSS(t1,u1,Param.LH1,Param),'abstol',accuracy,'reltol',accuracy);
    [t1,u1,~,~,~] = ode45(@(t,y)f_right_compact(t,y,Param), [t2(end) t2(end)+T2], [q1,qdot1],optionsSS);
    if t1(end) ~= t2(end)+T2 && u1(end,1) < (pi/2) && u1(end,1) > (-pi/2) && u1(end,3) < 0
        event = 'none';
    elseif  u1(end,3) > 0
        event = 'prematureLiftOff';
    else
        event = 'noTouchDown';
    end
    
    % [Xc,Yc] = Spring_Plots(t1,u1);
    
end
end


function [value,isterminal,direction] = foot_strikeSS(t,u,LH1,Param)
q1 = u(1);
q2 = u(2);
q3 = u(3);
q4 = 0;

[xth,yth,~,~,~,~]=xth_yth(q1,1,0,Param);
[xth2,yth2,~,~,~,~]=xth_yth(q1+q2,2,0,Param);

xH2 = 0;
yH1 = LH1 + q3;
yH2 = LH1 + q4;

xc =  (xth2-xH2)*cos(q2) + (yth2-yH2)*sin(q2);
yc = -(xth2-xH2)*sin(q2) + (yth2-yH2)*cos(q2) + yH1;

% Collision point at the end of the foot for the swing leg
Yc = -(xc-xth)*sin(q1) + (yc-yth)*cos(q1);
Y3 = xth*sin(q1) + (yH1 - yth)*cos(q1);

value = Yc;  % value of the ith event function

if u(4)+u(5) >= 0  && u(1) > 0 && u(2) < 0
    isterminal = 1;
elseif u(1) >= (pi/2) || u(1)<= (-pi/2) || u(3) > 0 % if the hip mass touches the floor, walker falls over.
    isterminal = 1;
    value = 0;
    direction = 0;
else
    isterminal = 0;
end
direction = -1;
end
% 
% function [value,isterminal,direction] = foot_strikeDS_NEW(t,u,Param)
% 
% % [~,q4,~,~] = constained_DS(u(1),u(2),u(3),u(4),Param.Xst);
% 
% [~, ~, ~, Y1ddot, ~, ~,~,U] = accelerationsDS_NEW(t2,u2,Param);
% 
% q3ddot = U(:,3);
% % Acceleration or rear foot
% Tacc = Y1ddot + q3ddot.*cos(u2(:,1));
% 
% % if loop here if there are any more conditions
% isterminal = 1;
% direction = 1;
% 
% if u(1) >= (pi/2) || u(1) <= (-pi/2)
%     value = 0;
%     direction = 0;
% end
% end

%% End condition - q4 = 0
function [value,isterminal,direction] = foot_strikeDS(t,u,Param)

% [~,q4,~,~] = constained_DS(u(1),u(2),u(3),u(4),Param.Xst);

q3 = u(3);
q3dot = u(7);

c = Param.c1;
k = Param.k1;

value = k*(q3 + Param.load) + c*q3dot;  % value of the ith event function

% if loop here if there are any more conditions
isterminal = 1;

direction = 1;
if u(1) >= (pi/2) || u(1) <= (-pi/2)
    value = 0;
    direction = 0;
end
end
