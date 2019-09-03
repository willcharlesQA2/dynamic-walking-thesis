function [Xc,Yc] = Spring_Plots(t,u)

%Plots trajectories of masses - Make this into a function
% Obtain lengths from model parameters
global Param
Lr1 = Param.Lr1;
Lr2 = Param.Lr2;
LH1 = Param.LH1;
xd2 = Param.xd2;
yd2 = Param.yd2;

S1th = Param.S1th;
S4th = Param.S4th;
% h = Param.h;
h = Param.h;

figure
stepsize = 1;

% initialise Xc and Yc
Xc = zeros(1,size(1:stepsize:size(t,1),2));
Yc = zeros(1,size(1:stepsize:size(t,1),2));

% Used to plot total roll-over shape
count = 1;
for i = linspace(S1th,S4th,50)
    [xthroll_total(count),ythroll_total(count),~,~,~,~]=xth_yth(i,1,0,Param);      % TOTAL rollover shape of stance foot
    [xth2roll_total(count),yth2roll_total(count),~,~,~,~]=xth_yth(i,2,0,Param);    % TOTAL rollover shape of swing foot
    count = count + 1;
end

for i = 1:stepsize:size(t)
    clf
    hold on
        
    q1 = u(i,1);
    q2 = u(i,2);
    q3 = u(i,3);
    
    if size(u,2) == 8
        % double support
        q4 = u(i,4);
    else
        %single support
        q4 = 0;
    end
        
%% Roll-over shape
[xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(q1,1,0,Param);
[xth2,yth2,~,~,~,~]=xth_yth(q1+q2,1,0,Param);

%% arc length
[ sth(i) ,~ ,~ ] = arcLength(q1,dxth,dyth,ddxth,ddyth,S1th,S4th,h,1);

    %% Lengths of each leg
    yr1 = Lr1 + q3;
    yH1 = LH1 + q3;
    yH2 = LH1 + q4;
    xH2 = 0;
    
    % Hip mass locations
    X3 = -xth*cos(q1) + (yH1 - yth)*sin(q1) + sth(i);
    Y3 = xth*sin(q1) + (yH1 - yth)*cos(q1);
    
    % Stance leg in SS. Trailing leg in DS.
    X1 = -xth*cos(q1) + (yr1 - yth)*sin(q1) + sth(i);
    Y1 = xth*sin(q1) + (yr1 - yth)*cos(q1);
    
    % Ankle joints
    X1_0=(-xth)*cos(q1)+(-yth)*sin(q1)+sth(i); % [x,y] at [0,0]
    Y1_0=-(-xth)*sin(q1)+(-yth)*cos(q1);
    
    xr0=(-yH2)*sin(q2); %[x',y'] at [0,0]
    yr0=(-yH2)*cos(q2)+yH1;
    
    Xsw0=(xr0-xth)*cos(q1)+(yr0-yth)*sin(q1)+sth(i); % [X,Y] locations at end of swing leg %
    Ysw0=-(xr0-xth)*sin(q1)+(yr0-yth)*cos(q1);     % USED TO PLOT LEGS

    % Swing leg in [x1,y1] coordinates
    xp2 =    (xd2)*cos(q2)+(yd2)*sin(q2);
    yp2 =    -(xd2)*sin(q2)+(yd2)*cos(q2)+yH1;
    
    % Swing leg in SS. Leading leg in DS.
    X2 = (xp2-xth)*cos(q1) + (yp2-yth)*sin(q1) + sth(i);
    Y2 = -(xp2-xth)*sin(q1) + (yp2-yth)*cos(q1);
    
    % End of swing foot in [x1,y1] coordinates
    xc =  (xth2-xH2)*cos(q2) + (yth2-yH2)*sin(q2);
    yc = -(xth2-xH2)*sin(q2) + (yth2-yH2)*cos(q2) + yH1;
    
    % End of swing leg in SS. Foot contact of lead leg in DS.
    Xc(i) = (xc-xth)*cos(q1) + (yc-yth)*sin(q1) + sth(i);
    Yc(i) = -(xc-xth)*sin(q1) + (yc-yth)*cos(q1);
    
    % Plot individual masses
    plot(X1,Y1,'g*')
    plot(X3,Y3,'b*')
    plot(X2,Y2,'r*')
    
    % Plot leg line from foot contact to hip mass
    plot([X1_0,X3],[Y1_0,Y3],'k')
    % Hip mass to end of the other foot
    plot([X3,Xsw0],[Y3,Ysw0],'r')
    
    %% PLOT TOTAL ROLLOVER
    % Plots rollover shape in [X,Y] coordinates
    Xroll_total=(xthroll_total-xth)*cos(q1)+(ythroll_total-yth)*sin(q1)+sth(i);
    Yroll_total=-(xthroll_total-xth)*sin(q1)+(ythroll_total-yth)*cos(q1);
    plot(Xroll_total,Yroll_total)
    
    % rollover shape of swing leg in [x,y] coordinates
    xth2roll_total_local=(xth2roll_total)*cos(q2)+(yth2roll_total-yH2).*sin(q2);        %same rollover shape as swing leg
    yth2roll_total_local=-(xth2roll_total)*sin(q2)+(yth2roll_total-yH2).*cos(q2)+yH1;
    
    % rollover shape of swing leg in [X,Y] coordinates
    Xrlsw_total=(xth2roll_total_local-xth)*cos(q1)+(yth2roll_total_local-yth)*sin(q1)+sth(i);
    Yrlsw_total=-(xth2roll_total_local-xth)*sin(q1)+(yth2roll_total_local-yth)*cos(q1);
    plot(Xrlsw_total,Yrlsw_total)
    
    %% Change axis. axisequal
    axis([-1,2,-1,2])
    
    % Controls live animation speed
    if i+10 <= size(t,1)
        pause((t(i+stepsize)-t(i))*1)
    end
end

fprintf('Xc(1) \t\t= %g\t\tYc(1) \t\t= %g\n',Xc(1),Yc(1))
fprintf('Xc(end) \t= %g\t\tYc(end) \t= %g\n\n',Xc(end),Yc(end))

totalChange = ( (Xc(end)-Xc(1))^2 + (Yc(end)-Yc(1))^2 )^0.5;

fprintf('Total change is %gmm\n\n',totalChange*1000)

end