function [t1,u1,t2,u2,qnew,event,IC] = oneStep(th_0,r_0,thdot_0,rdot_0,Param)
%One Step
% %   Detailed explanation goes here
% Param.tl = 0.2; % toe length
% Param.hl = 0; % heel length
% Param.fr = 1;   % foot radius
% Param.L0 = 1;
% Param.m = 80;
% Param.k = 12e3;
% Param.Aa = -20*(pi/180);
% % Minimum Theta, when leg becomes locked
% Param.Th = Param.hl/Param.fr;
% % Maximum Theta, when leg becomes locked
% Param.Tt = Param.tl/Param.fr;

t0 = 0;
T1 = 1;
T2 = 1;

% Initial Conditions

% Initial Conditions
% th_0 = -0.2806;
% r_0 = -0.01635;
%
% thdot_0 = 1.17;
% rdot_0 = -0.1743;

%% Single support phase
q1 = [th_0,r_0];
q1dot = [thdot_0,rdot_0];

optionsSS=odeset('Events',@(t1,u1)touchDown(t1,u1,Param),'abstol',10^(-8),'reltol',10^(-8));
optionsErr=odeset('abstol',10^(-12),'reltol',10^(-12));
[t1,u1] = ode45(@(t,y)SSdynamics(t,y,Param), [t0 t0+T1], [q1,q1dot],optionsSS);

if t1(end) == T1 || u1(end,1) >= (pi/2) || u1(end,1) <= (-pi/2) || u1(end,2) > 0
    event = 'noTouchDown';
    t2 = NaN;
    u2 = NaN;
    qnew = NaN;
    IC = NaN;
else
    
    q2 = [u1(end,1),u1(end,2)];
    q2dot = [u1(end,3),u1(end,4)];
    
    %% Calculate initial contact (IC) of touch down leg
    theta = u1(end,1);
    [ CoP,thC,~,~,~,~,~,~ ] = findCoP( theta,Param,1);
    
    xc =  Param.fr1*sin(thC);
    yc = -Param.fr1*cos(thC) + Param.fr1;
    L = Param.L0 + u1(end,2);
    X =  (-xc)*cos(theta) + (L-yc)*sin(theta) + CoP;
    Y = -(-xc)*sin(theta) + (L-yc)*cos(theta);
    
    [ CoP2,thC2,dCoP2,ddCoP2,dxc2,dyc2,ddxc2,ddyc2 ] = findCoP( Param.alpha,Param,2 );
    xc2 =  Param.fr2*sin(thC2);
    yc2 = -Param.fr2*cos(thC2) + Param.fr2;
    alpha = Param.alpha;
    L0 = Param.L0;
    Param.IC =  (xc2).*cos(alpha) + (yc2-L0).*sin(alpha) + X - CoP2; 
%     Yic = -(xc2).*sin(alpha) + (yc2-L0).*cos(alpha) + Y
    IC = Param.IC;
    %% Double - support phase
    
    optionsDS=odeset('Events',@(t2,u2)takeOff(t2,u2,Param),'abstol',10^(-8),'reltol',10^(-8));
    optionsErr=odeset('abstol',10^(-12),'reltol',10^(-12));
    [t2,u2] = ode45(@(t,y)DSdynamics(t,y,Param), [t1(end) t1(end)+T2], [q2,q2dot],optionsDS);
    
    if t2(end) == t1(end) + T2
        event = 'noTakeOff';
        qnew = NaN;
    else
        event = 'none';
        %% Now find initial conditions for next step (keeping Xdot and Ydot
        % constant)
        [ phi2,r2,th2dot,r2dot ] = newConditions( u2,Param );
        qnew = [ phi2,r2,th2dot,r2dot ];
    end
end
end

function [value,isterminal,direction] = takeOff(t,u,Param)
%% When ankle r1 = 0
% Rear spring r1
value = u(2);

isterminal = 1;

direction = 1;

if isnan(u(1))
    value = 0;
end

end


function [value,isterminal,direction] = touchDown(t,u,Param)
%% When front leg touches down


% if loop here if there are any more conditions
isterminal = 1;

% If walker falls over or there is early take-off
if u(1) >= (pi/2) || u(1) <= (-pi/2) || u(2) > 0
    value = 0;
    direction = 0;
elseif u(1) > 0
    L0 = Param.L0;
    fr1 = Param.fr1;
    theta = u(1);
    L = L0 + u(2);
    % angle of attack
    alpha = Param.alpha;
    
    [ CoP,thC,~,~,~,~,~,~ ] = findCoP( theta,Param,1 );
    
    xc =  fr1*sin(thC);
    yc = -fr1*cos(thC) + fr1;
    
%     X =  (-xc)*cos(theta) + (L-yc)*sin(theta) + CoP;
    Y = -(-xc)*sin(theta) + (L-yc)*cos(theta);
    % FrontLegExtension = (Y1^2-(Param.IC-X1)^2)^0.5;
    
    value = Y-L0*cos(alpha);  % value of the ith event function
    
    [ CoP2,thC2,dCoP2,ddCoP2,dxc2,dyc2,ddxc2,ddyc2 ] = findCoP( alpha,Param,2 );
    xc2 =  Param.fr2*sin(thC2);
    yc2 = -Param.fr2*cos(thC2) + Param.fr2;
    
%     Xic =  (xc2).*cos(alpha) + (yc2-L0).*sin(alpha) + X; 
    value = -(xc2).*sin(alpha) + (yc2-L0).*cos(alpha) + Y;

    direction = -1;
else
    isterminal = 0;
    value = -1;
    direction = 0;
end


end

