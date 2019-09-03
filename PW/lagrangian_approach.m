function [t,u,th1sw,th2sw,th1dotsw,th2dotsw,event,steptime,Dist]=lagrangian_approach(th1_0,th2_0,th1dot_0,th2dot_0,step)
%% LAGRANGIAN_APPROACH to gait modelling
% Calculates each step from heel strike to heel strike and double support
% transition phase. With legs of different masses, lengths etc, initial
% conditions need to be set from the parameter 'ini'.

%% INITIAL CONDITIONS
if nargin<5
    step=35;
end

if th2_0 < 0
    % This has happened when the walker has done the splits!
    fprintf('\n\color{red}theta_2 should be positive; initial conditions are wrong\n')
end

global x1 y1 m1 xr2 yr2 m2 xc yc mc alpha xh yh rh xm ym rm xf yf rf S1th S2th S3th S4th ini h

%     ini=[x_A y_A m_A x_B y_B m_B xc yc mc alpha...
%         rh_A rm_A rf_A h_A xh_A yh_A xm_A ym_A xf_A yf_A S1th_A S2th_A S3th_A S4th_A...
%         rh_B rm_B rf_B h_B xh_B yh_B xm_B ym_B xf_B yf_B S1th_B S2th_B S3th_B S4th_B];
%          11   12   13  14   15   16   17   18   19   20    21     22     23     24
%          25   26   27  28   29   30   31   32   33   34   35      36     37     38

x1=ini(1);
y1=ini(2);
m1=ini(3);
xr2=ini(4);          % x'2
yr2=ini(5);          % y'2
m2=ini(6);
xc=ini(7);
yc=ini(8);
mc=ini(9);

alpha=ini(10);     %slope angle

%%FOOT PARAMETERS

% s1=-0.03;      % hindfoot length
% s2=-0.01;       %
% s3=0.011;
% s4=0.14;       % forefoot length
rh=ini(11);      % hindfoot gain
rm=ini(12);      % midfoot gain
rf=ini(13);      % forefoot gain
h=ini(14);            % horizontal distance

%Calcuates theta(S) values
S2th=ini(22);
S3th=ini(23);

xm=ini(17);
ym=ini(18);
xh=ini(15);
yh=ini(16);
xf=ini(19);
yf=ini(20);

S1th=ini(21);
S4th=ini(24);


%% LEG POSITIONS AND VELOCITITES WITH NO INPUT ARGUMENTS
if nargin==0
    th1_0=-0.331065;
    th2_0=0.688414;
    th1dot_0=1.22061;
    th2dot_0=-0.00099235;
end
u0=[th1_0; th2_0; th1dot_0; th2dot_0];

T=1.2;            %FINAL TIME FOR SIMULATION %consider making adaptive to reduce simulation time

%% COMPUTATION
% options=odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5 1e-5 1e-5]);
options=odeset('Events',@(t,u)heel_strike2(t,u),'abstol',10^(-8),'reltol',10^(-8));
[t,u,TE,YE,IE] = ode45(@swing_stage, [0 T], u0, options);       % linspace(0,T,100)

% [energy]=check_energy(t,u);

% Calculates variables at heel strike using interpolation
%[th1_h,th2_h,th1dot_h,th2dot_h,heel,event,Xstrike,steptime]=heel_strike(t,u,step);  % NEEDS TO BE REPLACED!!
%%%%%%%%%%%%%%%%%%%%% FAILURE MODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xth,yth,~,~,~,~]=xth_yth(u(end,1),1);

Xcwithoutsth=(xc-xth)*cos(u(end,1))+(yc-yth)*sin(u(end,1));%+sth;
Yc=-(xc-xth)*sin(u(end,1))+(yc-yth)*cos(u(end,1));

% take out sth, as it is not needed and will take up lots of computing
% time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if time runss out or hip mass falls below floor (due to heel_srike2
% termination Yc is not always less than 0 but very close to it)
if t(end) == T || Yc <= 0.1
    if Xcwithoutsth < 0
        event = 'fallbackward';
    elseif Xcwithoutsth > 0
        event = 'fallforward';
    else
        event = 'none';
    end
elseif t(end) ~= T && Yc > 0
    event = 'hlstr';
else
    error('Problem with heel strike time')
end

switch event
    case {'none','fallforward','fallbackward'}
        %leave program
        %fprintf('No heel strike detected\n')
        th1sw=NaN;    th2sw=NaN;    th1dotsw=NaN; th2dotsw=NaN;
        steptime = NaN;
        Dist = NaN;
        %         plots(t,u,step);
    case 'hlstr'
        %% FIND STEP LENGTH
        
        % First find sth at start of stance phase (to find average velocity); then
        % find sth at end of stance phase (to find the point of contact at end of
        % stance phase).
        thi = 1;
        for theta_arc = [u(1,1),u(end,1)]
            if theta_arc<=S1th         %start point theta(S1)
                th_arc=S1th;
            elseif theta_arc>=S4th
                th_arc=S4th;      %end point theta(S4)
            else
                th_arc=theta_arc;
            end
            
            s0=h; %CHANGES WHEN h CHANGES
            if th_arc ~= 0
                [ths,s_arc]=ode45(@f_arc, [0 th_arc], s0);
            else
                s_arc=0;
            end
            sth(thi)=s_arc(end);
            
            [xth(thi),yth(thi),~,~,~,~]=xth_yth(theta_arc,1);      %STANCE LEG ROLLOVER SHAPE
            thi = thi+1;        % Should be 2 at arc length end
        end
        
        Xstart = sth(1);
        
        [xth2,yth2,~,~,~,~]=xth_yth(u(end,1)+u(end,2),2); %   SWING LEG ROLLOVER SHAPE
        xrp=(xth2-xc)*cos(u(end,2))+(yth2-yc)*sin(u(end,2))+xc;  % contact point in [x,y]
        yrp=-(xth2-xc)*sin(u(end,2))+(yth2-yc)*cos(u(end,2))+yc; % reference frames.
        
        Xstrike=(xrp-xth(2))*cos(u(end,1))+(yrp-yth(2))*sin(u(end,1))+sth(2);
        
        %% Not sure what this stuff is
        %     [xth(2),yth(2),~,~,~,~]=xth_yth(0.410086,1);      %STANCE LEG ROLLOVER SHAPE
        %
        %     [xth2,yth2,~,~,~,~]=xth_yth(0.410086-0.26841,2); %   SWING LEG ROLLOVER SHAPE
        %     xrp=(xth2-xc)*cos(-0.26841)+(yth2-yc)*sin(-0.26841)+xc;  % contact point in [x,y]
        %     yrp=-(xth2-xc)*sin(-0.26841)+(yth2-yc)*cos(-0.26841)+yc; % reference frames.
        %     X2strike=(xrp-xth(2))*cos(0.410086)+(yrp-yth(2))*sin(0.410086)+sth(2);
        
        %Distance travelled
        Dist = Xstrike - Xstart;
        
        
        %% Values at heel-strike
        th1_h = u(end,1);
        th2_h = u(end,2);
        th1dot_h = u(end,3);
        th2dot_h = u(end,4);
        steptime = t(end);
        heel = size(t);     % i at which heel strike occurs. Not necessary with 'Events'(heel_strike2)
        
        
        
        %% SWITCH CONDITIONS AT HEEL STRIKE
        
        Phi=momentum(th1_h,th2_h,th1dot_h,th2dot_h,Xstrike);
        %Phi2=momentum2(0.410086,-0.7134,1.5239,-1.15886,X2strike);
        
        th1sw=th1_h+th2_h;      %becomes th1 for next step
        th2sw=-th2_h;           %becomes th2
        th1dotsw=Phi(1);      %becomes th1dot
        th2dotsw=Phi(2);      %becomes th2dot
        
        %         fprintf('Conditions at heel strike, step %g:\nStep Period=%gs\nth1_0\t=\t%g;\nth2_0\t=\t%g;\nth1dot_0=\t%g;\nth2dot_0=\t%g;\n\n',step,steptime,th1sw,th2sw,th1dotsw,th2dotsw)
        
end





end


function [du]=swing_stage(t,u)

global x1 y1 m1 xr2 yr2 m2 xc yc mc alpha S1th S4th h

%Potential


%Pi=m1*g*(X1*sin(alpha)+Y1*cos(alpha))+mc*g*(Xc*sin(alpha)+Yc*cos(alpha))+m2*g*(X2*sin(alpha)+Y2*cos(alpha));
%coordinates xth,yth etc.
%INITIAL VARIABLES

g=9.81;

th1=u(1);
th2=u(2);

[xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(th1,1);

% xth=tan(th1)*(r/2);
% yth=r/4*tan(th1)^2;
% dxth=sec(th1)^2*(r/2);
% dyth=sec(th1)^2*xth;
% ddxth=(r/2)*2*sec(th1)^2*tan(th1);
% ddyth=dxth*sec(th1)^2+2*sec(th1)^2*tan(th1)*(xth);

%%CALCULATE ARC LENGTH
if th1<=S1th          %start point theta(S1)
    th_arc=S1th;
elseif th1>=S4th
    th_arc=S4th;      %end point theta(S4)
else
    th_arc=th1;
end

%s0=h;        %s0=(0+xm^2)^0.5;
% if th_arc ~=0
%     [ths,s_arc]=ode45(@f_arc, [0 th_arc], s0);
% else
%     s_arc=0;
% end
%
% sth=s_arc(end);
if th_arc == S1th || th_arc==S4th
    dsth=0; ddsth=0;
else
    dsth=(dxth^2+dyth^2)^0.5;
    ddsth=(dxth*ddxth+dyth*ddyth)/(dxth^2+dyth^2)^0.5;
end

% FOR POINT FOOT:
% xth=0;  yth=0;  dxth=0; dyth=0; ddxth=0;    ddyth=0;

% X1=(x1-xth)*cos(th1)+(y1-yth)*sin(th1)+sth;
% Y1=-(x1-xth)*sin(th1)+(y1-yth)*cos(th1);
% %
% Xcwithoutsth=(xc-xth)*cos(th1)+(yc-yth)*sin(th1);%+sth;
% Yc=-(xc-xth)*sin(th1)+(yc-yth)*cos(th1);
%
%  %% USED FOR FAILURE MODES
%  % take out sth, as it is not needed and will take up lots of computing
%  % time
%  if Yc < 0 && Xcwithoutsth < 0
%      fprintf('\nWALKER FALLS BACKWARDS')
%  elseif Yc < 0 && Xcwithoutsth > 0
%         fprintf('\nWALKER FALLS FORWARDS')
%  end
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
x2=(xr2-xc)*cos(th2)+(yr2-yc)*sin(th2)+xc;
y2=-(xr2-xc)*sin(th2)+(yr2-yc)*cos(th2)+yc;
%
% X2=(x2-xth)*cos(th1)+(y2-yth)*sin(th1)+sth;
% Y2=-(x2-xth)*sin(th1)+(y2-yth)*cos(th1);


%Kinetic
dX1dth1=-dxth*cos(th1)-(x1-xth)*sin(th1)-dyth*sin(th1)+(y1-yth)*cos(th1)+dsth;
dY1dth1=+dxth*sin(th1)-(x1-xth)*cos(th1)-dyth*cos(th1)-(y1-yth)*sin(th1);%changed

dXcdth1=-dxth*cos(th1)-(xc-xth)*sin(th1)-dyth*sin(th1)+(yc-yth)*cos(th1)+dsth;
dYcdth1=+dxth*sin(th1)-(xc-xth)*cos(th1)-dyth*cos(th1)-(yc-yth)*sin(th1);%changed

% dX2dth1=-dXcdth1*cos(th2)-dYcdth1*sin(th2)+dXcdth1;
% dY2dth1=+dXcdth1*sin(th2)-dYcdth1*cos(th2)+dYcdth1;
%
% dX2dth2=-(x2-Xc)*sin(th2)+(y2-Yc)*cos(th2);
% dY2dth2=-(x2-Xc)*cos(th2)-(y2-Yc)*sin(th2);

dX1dth2=0;
dY1dth2=0;
dXcdth2=0;
dYcdth2=0;

ddX1ddth1 = -ddxth*cos(th1)+2*dxth*sin(th1)-(x1-xth)*cos(th1)-...
    ddyth*sin(th1)-2*dyth*cos(th1)-(y1-yth)*sin(th1)+...
    ddsth;
ddY1ddth1 = +2*dxth*cos(th1)+ddxth*sin(th1)+(x1-xth)*sin(th1)+...
    -ddyth*cos(th1)+2*dyth*sin(th1)-(y1-yth)*cos(th1); %changed

ddXcddth1 = -ddxth*cos(th1)+2*dxth*sin(th1)-(xc-xth)*cos(th1)-...
    ddyth*sin(th1)-2*dyth*cos(th1)-(yc-yth)*sin(th1)+...
    ddsth;
ddYcddth1 = +2*dxth*cos(th1)+ddxth*sin(th1)+(xc-xth)*sin(th1)+...
    -ddyth*cos(th1)+2*dyth*sin(th1)-(yc-yth)*cos(th1); %changed


%%
% x2 =    (xr2-xc)*cos(th2)+(yr2-yc)*sin(th2)+xc;
% y2 =    -(xr2-xc)*sin(th2)+(yr2-yc)*cos(th2)+yc;

% figure(5)
% hold on
% plot(x2,y2,'*')
%
% dx2dth1=0;
% dy2dth1=0;
%
% ddx2ddth1=0;
% ddy2ddth1=0;

dx2dth2 =   -(xr2-xc)*sin(th2)+(yr2-yc)*cos(th2);
dy2dth2 =   -(xr2-xc)*cos(th2)-(yr2-yc)*sin(th2);

ddx2ddth2 = -(xr2-xc)*cos(th2)-(yr2-yc)*sin(th2);
ddy2ddth2 = +(xr2-xc)*sin(th2)-(yr2-yc)*cos(th2);

% ddx2dth1dth2 =  0;
% ddy2dth1dth2 =  0;

% X2=(x2-xth)*cos(th1)+(y2-yth)*sin(th1)+sth;
% Y2=-(x2-xth)*sin(th1)+(y2-yth)*cos(th1);

dX2dth1 =   -(x2-xth)*sin(th1)-dxth*cos(th1)+...
    (y2-yth)*cos(th1)-dyth*sin(th1)+dsth;
dY2dth1 =   -(x2-xth)*cos(th1)+dxth*sin(th1)-...
    (y2-yth)*sin(th1)-dyth*cos(th1);

ddX2ddth1 = -(x2-xth)*cos(th1)+2*dxth*sin(th1)-ddxth*cos(th1)-...
    (y2-yth)*sin(th1)-2*dyth*cos(th1)-ddyth*sin(th1)+ddsth;

ddY2ddth1 = (x2-xth)*sin(th1)+2*dxth*cos(th1)+ddxth*sin(th1)-...
    (y2-yth)*cos(th1)+2*dyth*sin(th1)-ddyth*cos(th1);

dX2dth2 =   dx2dth2*cos(th1)+dy2dth2*sin(th1);
dY2dth2 =   -dx2dth2*sin(th1)+dy2dth2*cos(th1);

ddX2ddth2 = ddx2ddth2*cos(th1)+ddy2ddth2*sin(th1);
ddY2ddth2 = -ddx2ddth2*sin(th1)+ddy2ddth2*cos(th1);

ddX2dth1dth2 =  -dx2dth2*sin(th1)+...
    +dy2dth2*cos(th1);
ddY2dth1dth2 =  -dx2dth2*cos(th1)+...
    -dy2dth2*sin(th1);

% ddX2dth1dth2 =  ddx2dth1dth2*cos(th1)-dx2dth2*sin(th1)+...
%                 ddy2dth1dth2*sin(th1)+dy2dth2*cos(th1);
% ddY2dth1dth2 =  -ddx2dth1dth2*sin(th1)-dx2dth2*cos(th1)+...
%                 ddy2dth1dth2*cos(th1)-dy2dth2*sin(th1);

%%

dPidth1 =   m1*g*(dX1dth1*sin(alpha)+dY1dth1*cos(alpha))+...
    mc*g*(dXcdth1*sin(alpha)+dYcdth1*cos(alpha))+...        %POTENTIAL ENERGY
    m2*g*(dX2dth1*sin(alpha)+dY2dth1*cos(alpha));
dPidth2 =   m1*g*(dX1dth2*sin(alpha)+dY1dth2*cos(alpha))+...
    mc*g*(dXcdth2*sin(alpha)+dYcdth2*cos(alpha))+...
    m2*g*(dX2dth2*sin(alpha)+dY2dth2*cos(alpha));

dTh1dth1 =  m1*2*(ddX1ddth1*dX1dth1+ddY1ddth1*dY1dth1)+...
    mc*2*(ddXcddth1*dXcdth1+ddYcddth1*dYcdth1)+...
    m2*2*(ddX2ddth1*dX2dth1+ddY2ddth1*dY2dth1);
dTh1dth2 =  m2*2*(ddX2dth1dth2*dX2dth1+ddY2dth1dth2*dY2dth1);   %this

dTh12dth1 = m2*(ddX2ddth1*dX2dth2+dX2dth1*ddX2dth1dth2+...
    ddY2ddth1*dY2dth2+dY2dth1*ddY2dth1dth2);
dTh12dth2 = m2*(ddX2dth1dth2*dX2dth2+dX2dth1*ddX2ddth2+...      %this
    ddY2dth1dth2*dY2dth2+dY2dth1*ddY2ddth2);

dTh2dth1 =  m2*2*(dX2dth2*ddX2dth1dth2 + dY2dth2*ddY2dth1dth2);
dTh2dth2 =  m2*2*(dX2dth2*ddX2ddth2 + dY2dth2*ddY2ddth2);

Theta1 =    m1*(dX1dth1^2 + dY1dth1^2) + mc*(dXcdth1^2 + dYcdth1^2) + m2*(dX2dth1^2 + dY2dth1^2);
Theta12 =   m2*(dX2dth1*dX2dth2 + dY2dth1*dY2dth2);
Theta2 =    m2*(dX2dth2^2 + dY2dth2^2);

% dTh12dth2=0;
th1dot=u(3);
th2dot=u(4);

delta = Theta1*Theta2-Theta12^2;
% f1 =    +1/2*dTh1dth1*u(3)^2+dTh12dth1*u(3)*u(4)+1/2*dTh2dth1*u(4)^2-dPidth1;
% f2 =    +1/2*dTh1dth2*u(3)^2+dTh12dth2*u(3)*u(4)+1/2*dTh2dth2*u(4)^2-dPidth2;

% f1 =    +1/2*dTh1dth1*u(3)^2+dTh12dth1*u(3)*u(4)+1/2*dTh2dth1*u(4)^2-dPidth1...
%     -dTh1dth1*th1dot^2-(dTh12dth1+dTh1dth2)*th1dot*th2dot-dTh12dth2*th2dot^2;
% f2 =    +1/2*dTh1dth2*u(3)^2+dTh12dth2*u(3)*u(4)+1/2*dTh2dth2*u(4)^2-dPidth2...
%     -dTh12dth1*th1dot^2-(dTh2dth1+dTh12dth2)*th1dot*th2dot-dTh2dth2*th2dot^2;

f1 =    -1/2*dTh1dth1*u(3)^2-dTh1dth2*u(3)*u(4)+(-dTh12dth2+1/2*dTh2dth1)*u(4)^2-dPidth1;

f2 =    (-dTh12dth1+1/2*dTh1dth2)*u(3)^2-dTh2dth1*u(3)*u(4)-1/2*dTh2dth2*u(4)^2-dPidth2;


% % Compare energy in system
% K   =   0.5*Theta1*u(3)^2+Theta12*u(3)*u(4)+0.5*Theta2*u(4)^2;
% Pi  =   m1*g*(X1*sin(alpha)+Y1*cos(alpha))+...
%         mc*g*(Xc*sin(alpha)+Yc*cos(alpha))+...
%         m2*g*(X2*sin(alpha)+Y2*cos(alpha));
%
% energy=K+Pi

%% Integrates the following equation
% fprintf('\n%g\n%g\n%g


du = [ u(3);
    u(4);
    1/delta*(f1*Theta2-f2*Theta12);
    1/delta*(-f1*Theta12+f2*Theta1) ];

end

function [value,isterminal,direction] = heel_strike2(t,u)

global ini xc yc
h=ini(14);
S1th=ini(21);
S4th=ini(24);

th1=u(1);
th2=u(2);

%     if th1<S1th         %start point theta(S1)
%         th_arc=S1th;
%     elseif th1>=S4th
%         th_arc=S4th;      %end point theta(S4)
%     else
%         th_arc=th1;
%     end
%
%     s0=h; %CHANGES WHEN h CHANGES
%     if th_arc ~= 0
%         [ths,s_arc]=ode45(@f_arc, [0 th_arc], s0);
%     else
%         s_arc=0;
%     end
%     sth=s_arc(end);

[xth,yth,~,~,~,~]=xth_yth(th1,1);
[xth2,yth2,~,~,~,~]=xth_yth(th1+th2,2); %   SWING LEG ROLLOVER SHAPE
xrp=(xth2-xc)*cos(th2)+(yth2-yc)*sin(th2)+xc;  % contact point in [x,y]
yrp=-(xth2-xc)*sin(th2)+(yth2-yc)*cos(th2)+yc; % reference frames.

%Xswp=(xrp-xth)*cos(th1)+(yrp-yth)*sin(th1)+sth; %X location of swing leg at contact point %CHANGE STH
Yswp=-(xrp-xth)*sin(th1)+(yrp-yth)*cos(th1);      %Y location of swing leg at contact point


% value(i) is the value of the function.
% isterminal(i) = 1, if the integration is to terminate at a zero of th1s
%   event function and 0 otherwise.
% direction(i) = 0 if all zeros are to be computed (the default), +1 if
%   only the zeros where the event function increases, and -1 if only the
%   zeros where the event function decreases.
value = Yswp;

if u(3)+u(4) >= 0 && u(1) > 0 && u(1)+u(2) < 0 % Heel strike occurs at these conditions.
    isterminal = 1; % Change termination condition
else
    isterminal = 0;
end

%% USED FOR FAILURE MODES
Yc=-(xc-xth)*sin(th1)+(yc-yth)*cos(th1);
if Yc < 0
    value = 0;
    isterminal =1;
end

direction = -1;
end



