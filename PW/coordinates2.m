function [X1,Y1,Xc,Yc,X2,Y2,dX1dth1,dXcdth1,dX2dth1,dX2dth2,...
    dY1dth1,dYcdth1,dY2dth1,dY2dth2,sth] = coordinates2(th1,th2,leg)
%COORDINATES2 like 'coordinates' but with sth and ddsth removed, in case
%the pivot point is moving
%   Converts th1 & th2 into X,Y coordinates with the origin at point of 
%   contact with the floor.
%   Is also used for the momentum equations

global ini

[xth,yth,dxth,dyth,~,~]=xth_yth(th1,leg);   % if leg==1: STANCE LEG
                                                    %    leg==2: SWING LEG
if leg==1               % STANCE LEG
    x1=ini(1);
    y1=ini(2);
    xr2=ini(4);
    yr2=ini(5);
    xc=ini(7);
    yc=ini(8);
    h=ini(14);
    S1th=ini(21);
    S4th=ini(24); 
elseif leg==2           % SWING LEG
    x1=ini(4);
    y1=ini(5);
    xr2=ini(1);
    yr2=ini(2);
    xc=ini(7);
    yc=ini(8);
    h=ini(28);
    S1th=ini(35);
    S4th=ini(38); 
end



%%CALCULATE ARC LENGTH                              
if th1<=S1th          %start point theta(S1)
    th_arc=S1th;
elseif th1>=S4th
    th_arc=S4th;      %end point theta(S4)
else
    th_arc=th1;
end

s0=h;        %s0=(0+xm^2)^0.5;
if th_arc ~=0
    [ths,s_arc]=ode45(@f_arc, [0 th_arc], s0); 
else
    s_arc=0;
end

sth=s_arc(end);         % arc length

X1=(x1-xth)*cos(th1)+(y1-yth)*sin(th1)+sth; 
Y1=-(x1-xth)*sin(th1)+(y1-yth)*cos(th1);
% 
Xc=(xc-xth)*cos(th1)+(yc-yth)*sin(th1)+sth; 
Yc=-(xc-xth)*sin(th1)+(yc-yth)*cos(th1);
% 
x2=(xr2-xc)*cos(th2)+(yr2-yc)*sin(th2)+xc;
y2=-(xr2-xc)*sin(th2)+(yr2-yc)*cos(th2)+yc;
%
X2=(x2-xth)*cos(th1)+(y2-yth)*sin(th1)+sth;
Y2=-(x2-xth)*sin(th1)+(y2-yth)*cos(th1);

dX1dth1=-dxth*cos(th1)-(x1-xth)*sin(th1)-dyth*sin(th1)+(y1-yth)*cos(th1);
dY1dth1=+dxth*sin(th1)-(x1-xth)*cos(th1)-dyth*cos(th1)-(y1-yth)*sin(th1);

dXcdth1=-dxth*cos(th1)-(xc-xth)*sin(th1)-dyth*sin(th1)+(yc-yth)*cos(th1);
dYcdth1=+dxth*sin(th1)-(xc-xth)*cos(th1)-dyth*cos(th1)-(yc-yth)*sin(th1);

dx2dth2 =   -(xr2-xc)*sin(th2)+(yr2-yc)*cos(th2);
dy2dth2 =   -(xr2-xc)*cos(th2)-(yr2-yc)*sin(th2);

% ddx2ddth2 = -(xr2-xc)*cos(th2)-(yr2-yc)*sin(th2); 
% ddy2ddth2 = +(xr2-xc)*sin(th2)-(yr2-yc)*cos(th2);

dX2dth1 =   -(x2-xth)*sin(th1)-dxth*cos(th1)+...
            (y2-yth)*cos(th1)-dyth*sin(th1);
dY2dth1 =   -(x2-xth)*cos(th1)+dxth*sin(th1)-...
            (y2-yth)*sin(th1)-dyth*cos(th1);

dX2dth2 =   dx2dth2*cos(th1)+dy2dth2*sin(th1); 
dY2dth2 =   -dx2dth2*sin(th1)+dy2dth2*cos(th1);

end

