function [K,Pi]=check_energy(t,u)

global x1 y1 m1 xr2 yr2 m2 xc yc mc alpha S1th S4th h

for i = 1: size(t)

    g=9.81;

    th1=u(i,1);
    th2=u(i,2);

    [xth,yth,dxth,dyth,~,~]=xth_yth(th1,1);

    %%CALCULATE ARC LENGTH
    if th1<=S1th          %start point theta(S1)
        th_arc=S1th;
    elseif th1>=S4th
        th_arc=S4th;      %end point theta(S4)
    else
        th_arc=th1;
    end

    s0=h;        %s0=(0+xm^2)^0.5;
    if th_arc ~= 0
        [~,s_arc]=ode45(@f_arc, [0 th_arc], s0); 
    else
        s_arc=0;
    end

    sth=s_arc(end);
    dsth=(dxth^2+dyth^2)^0.5;

    %%LOCATIONS
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

    %% DIFFERENTIATE X1, Xc
    dX1dth1=-dxth*cos(th1)-(x1-xth)*sin(th1)-dyth*sin(th1)+(y1-yth)*cos(th1)+dsth;
    dY1dth1=+dxth*sin(th1)-(x1-xth)*cos(th1)-dyth*cos(th1)-(y1-yth)*sin(th1);%changed

    dXcdth1=-dxth*cos(th1)-(xc-xth)*sin(th1)-dyth*sin(th1)+(yc-yth)*cos(th1)+dsth;
    dYcdth1=+dxth*sin(th1)-(xc-xth)*cos(th1)-dyth*cos(th1)-(yc-yth)*sin(th1);%changed

    % ddX1ddth1 = -ddxth*cos(th1)+2*dxth*sin(th1)-(x1-xth)*cos(th1)-...
    %             ddyth*sin(th1)-2*dyth*cos(th1)-(y1-yth)*sin(th1)+...
    %             ddsth; 
    % ddY1ddth1 = +2*dxth*cos(th1)+ddxth*sin(th1)+(x1-xth)*sin(th1)+...
    %             -ddyth*cos(th1)+2*dyth*sin(th1)-(y1-yth)*cos(th1); %changed
    % 
    % ddXcddth1 = -ddxth*cos(th1)+2*dxth*sin(th1)-(xc-xth)*cos(th1)-...
    %             ddyth*sin(th1)-2*dyth*cos(th1)-(yc-yth)*sin(th1)+...
    %             ddsth; 
    % ddYcddth1 = +2*dxth*cos(th1)+ddxth*sin(th1)+(xc-xth)*sin(th1)+...
    %             -ddyth*cos(th1)+2*dyth*sin(th1)-(yc-yth)*cos(th1); %changed

    %% DIFFERENTIATE X2

    dx2dth2 =   -(xr2-xc)*sin(th2)+(yr2-yc)*cos(th2);
    dy2dth2 =   -(xr2-xc)*cos(th2)-(yr2-yc)*sin(th2);

    % ddx2ddth2 = -(xr2-xc)*cos(th2)-(yr2-yc)*sin(th2);
    % ddy2ddth2 = +(xr2-xc)*sin(th2)-(yr2-yc)*cos(th2);

    dX2dth1 =   -(x2-xth)*sin(th1)-dxth*cos(th1)+...        %CHECK SIGNS!!!
                (y2-yth)*cos(th1)-dyth*sin(th1)+dsth;
    dY2dth1 =   -(x2-xth)*cos(th1)+dxth*sin(th1)-...
                (y2-yth)*sin(th1)-dyth*cos(th1);

    % ddX2ddth1 = -(x2-xth)*cos(th1)+2*dxth*sin(th1)-ddxth*cos(th1)-...
    %     (y2-yth)*sin(th1)-2*dyth*cos(th1)-ddyth*sin(th1)+ddsth;
    % 
    % ddY2ddth1 = (x2-xth)*sin(th1)+2*dxth*cos(th1)+ddxth*sin(th1)-...
    %     (y2-yth)*cos(th1)+2*dyth*sin(th1)-ddyth*cos(th1);

    dX2dth2 =   dx2dth2*cos(th1)+dy2dth2*sin(th1); 
    dY2dth2 =   -dx2dth2*sin(th1)+dy2dth2*cos(th1);

    % ddX2ddth2 = ddx2ddth2*cos(th1)+ddy2ddth2*sin(th1); 
    % ddY2ddth2 = -ddx2ddth2*sin(th1)+ddy2ddth2*cos(th1);
    %  
    % ddX2dth1dth2 =  ddx2dth1dth2*cos(th1)-dx2dth2*sin(th1)+...
    %                 ddy2dth1dth2*sin(th1)+dy2dth2*cos(th1); 
    % ddY2dth1dth2 =  -ddx2dth1dth2*sin(th1)-dx2dth2*cos(th1)+...
    %                 ddy2dth1dth2*cos(th1)-dy2dth2*sin(th1);

    %%CALCULATE THETAS            
    Theta1 =    m1*(dX1dth1^2 + dY1dth1^2) + mc*(dXcdth1^2 + dYcdth1^2) + m2*(dX2dth1^2 + dY2dth1^2);
    Theta12 =   m2*(dX2dth1*dX2dth2 + dY2dth1*dY2dth2);
    Theta2 =    m2*(dX2dth2^2 + dY2dth2^2);    

    K(i)  =   0.5*Theta1*u(i,3)^2+Theta12*u(i,3)*u(i,4)+0.5*Theta2*u(i,4)^2;
    Pi(i) =   m1*g*(X1*sin(alpha)+Y1*cos(alpha))+...
              mc*g*(Xc*sin(alpha)+Yc*cos(alpha))+...
              m2*g*(X2*sin(alpha)+Y2*cos(alpha));

% L(i)=K(i)-Pi(i)
    
    
    
end
