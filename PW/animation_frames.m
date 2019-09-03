function [M,leftover_time] = animation_frames(t,u,step,tstart)
%% ANIMATES WALKER. % animates a few steps at a time
% Consider using getframe and movie functions in order to ger a smoother
% animation
% xlabel('$\dot{\theta_{1}}$','interpreter','latex')

% Figure1=figure(1);clf;
% set(Figure1,'defaulttextinterpreter','latex');
% plot(...);
% xlabel('$\delta$');
% ylabel('$\epsilon$');

% For movie try using interp1(t,th1,snap) where snap is the time at which
% to take a screenshot.

global r m1 m2  mc alpha S1th S4th h  x1 xr2 xc y1 yr2 yc
slope = -alpha;
if nargin < 4
    tstart = 0;
end

% xth=tan(th1).*(r/2);
% yth=r/4.*tan(th1).^2;
% dxth=sec(th1).^2.*(r/2);
% dyth=(sec(th1).^2).*xth;
% xth=xth';
% yth=yth';

for i = 1:size(t,1)
    [xthroll(i),ythroll(i),~,~,~,~]=xth_yth(u(i,1),1);      % USED rollover shape of stance foot
    [xth2roll(i),yth2roll(i),~,~,~,~]=xth_yth(u(i,1),2);    % USED rollover shape of swing foot
end

count = 1;
for i = linspace(S1th,S4th,50)
    [xthroll_total(count),ythroll_total(count),~,~,~,~]=xth_yth(i,1);      % TOTAL rollover shape of stance foot
    [xth2roll_total(count),yth2roll_total(count),~,~,~,~]=xth_yth(i,2);    % TOTAL rollover shape of swing foot
    count = count + 1;
end


%% CALCULATE ARC LENGTH
%tstart = 0;     %changes with more steps
figure(1)
set(gcf, 'Renderer', 'zbuffer');
axis([-0.7 0.7 -0.2 1.2])
frame=1;        % used for movie loop
fps = 20;       % frames per second (hz)
delta_t = 1 / fps; 
 
for time=(delta_t-tstart):delta_t:t(end) % (delta_t-tstart):delta_t:t(end)
    
th1 = interp1(t,u(:,1),time);
th2 = interp1(t,u(:,2),time);
th1dot = interp1(t,u(:,3),time);
th2dot = interp1(t,u(:,4),time);    
    
[xth,yth,~,~,~,~]=xth_yth(th1,1);
[xth2,yth2,~,~,~,~]=xth_yth(th1,2);      %leg2    
    
    
    thi=th1;
    if thi<=S1th         %start point theta(S1)
        th_arc=S1th;
    elseif thi>=S4th
        th_arc=S4th;      %end point theta(S4)
    else
        th_arc=thi;
    end

    s0=h; %CHANGES WHEN h CHANGES
    if th_arc ~= 0
        [ths,s_arc]=ode45(@f_arc, [0 th_arc], s0); 
    else
        s_arc=0;
    end
    
    sth=s_arc(end);
    
    % DISTANCE TRAVELLED = -sth(1) + Xsw

 

%%
% xth=0;
% yth=0;
% Location of masses in [X,Y] coordinates
% X1=(x1-xth).*cos(th1)+(y1-yth).*sin(th1)+sth'; 
% Y1=-(x1-xth).*sin(th1)+(y1-yth).*cos(th1);
% 
% Xc=(xc-xth).*cos(th1)+(yc-yth).*sin(th1)+sth'; 
% Yc=-(xc-xth).*sin(th1)+(yc-yth).*cos(th1);
% 
% x2=(xr2-xc).*cos(th2)+(yr2-yc).*sin(th2)+xc;
% y2=-(xr2-xc).*sin(th2)+(yr2-yc).*cos(th2)+yc;
% 
% X2=(x2-xth).*cos(th1)+(y2-yth).*sin(th1)+sth';
% Y2=-(x2-xth).*sin(th1)+(y2-yth).*cos(th1);

% X_G = (m1*X1 + m2*X2 + mc*Xc)/(m1+m2+mc);
% Y_G = (m1*Y1 + m2*Y2 + mc*Yc)/(m1+m2+mc);

% PLOT TRAJECTORY OF MASSES

%title(step)
%hold on
%axis equal
%grid on
% plot(X1,Y1)            % trajectory of stance leg mass
% plot(Xc,Yc)            % trajectory of hip mass
% plot(X2,Y2,'--')       % trajectory of swing leg mass
% plot(X_G,Y_G,'k')
% plot(X1(1),Y1(1),'g*') % start points
% plot(Xc(1),Yc(1),'*')  % '' '' 
% plot(X2(1),Y2(1),'r*') % '' ''
%plot(sth,0,'g.')       % contact point with the floor
%plot([X1(1),Xc(1)],[Y1(1),Yc(1)])

%%PLOT LEG LINES




%%


% tic
%delta = 10;      % animation framerate

% for i=1:delta:size(t,1)
    

        
    if mod(step,2)              % Physiological
        swing = 'r';
        stance = 'g';
    else                        % Prosthetic
        swing = 'g';
        stance = 'r';
    end

X1=(x1-xth)*cos(th1)+(y1-yth)*sin(th1)+sth; 
Y1=-(x1-xth)*sin(th1)+(y1-yth)*cos(th1);

Xc=(xc-xth)*cos(th1)+(yc-yth)*sin(th1)+sth; 
Yc=-(xc-xth)*sin(th1)+(yc-yth)*cos(th1);

x2=(xr2-xc)*cos(th2)+(yr2-yc)*sin(th2)+xc;
y2=-(xr2-xc)*sin(th2)+(yr2-yc)*cos(th2)+yc;

X2=(x2-xth).*cos(th1)+(y2-yth).*sin(th1)+sth;
Y2=-(x2-xth).*sin(th1)+(y2-yth).*cos(th1);

xr0=(-xc)*cos(th2)+(-yc)*sin(th2)+xc; %[x',y'] at [0,0]
yr0=-(-xc)*sin(th2)+(-yc)*cos(th2)+yc;
%%
Xsw0=(xr0-xth)*cos(th1)+(yr0-yth)*sin(th1)+sth; % [X,Y] locations at end of swing leg %
Ysw0=-(xr0-xth)*sin(th1)+(yr0-yth)*cos(th1);     % USED TO PLOT LEGS
%%
X1_0=(-xth)*cos(th1)+(-yth)*sin(th1)+sth; % [x,y] at [0,0]
Y1_0=-(-xth)*sin(th1)+(-yth)*cos(th1);

centre = Xc*cos(alpha) - Yc*sin(alpha);

    %disp(t(i))
    plot([Xsw0-centre,Xc-centre]*cos(slope)+[Ysw0,Yc]*sin(slope),[Xsw0-centre,Xc-centre]*-sin(slope)+[Ysw0,Yc]*cos(slope),swing,'LineWidth',2)       % swing leg
    hold on
    %grid on
    axis([-0.7 0.7 -0.2 1.2])
    plot([X1_0-centre,Xc-centre]*cos(slope)+[Y1_0,Yc]*sin(slope),[X1_0-centre,Xc-centre]*-sin(slope)+[Y1_0,Yc]*cos(slope),stance,'LineWidth',2)       % stance leg
    
    plot((X1-centre)*cos(slope)+Y1*sin(slope),(X1-centre)*-sin(slope)+Y1*cos(slope),'ko','MarkerSize',10,'MarkerFaceColor',stance) % location of masses
    plot((Xc-centre)*cos(slope)+Yc*sin(slope),(Xc-centre)*-sin(slope)+Yc*cos(slope),'ko','MarkerSize',10,'MarkerFaceColor','b')  % '' '' 
    plot((X2-centre)*cos(slope)+Y2*sin(slope),(X2-centre)*-sin(slope)+Y2*cos(slope),'ko','MarkerSize',10,'MarkerFaceColor',swing) % '' ''
%     plot(Xsw0(i),Ysw0(i),'r*')
%% PLOT FEET CONTACT USED    
    % Plots rollover shape in [X,Y] coordinates
    Xroll=(xthroll-xth)*cos(th1)+(ythroll-yth)*sin(th1)+sth;
    Yroll=-(xthroll-xth)*sin(th1)+(ythroll-yth)*cos(th1);
    plot((Xroll-centre)*cos(slope)+Yroll*sin(slope),(Xroll-centre)*-sin(slope)+Yroll*cos(slope))
    
    % rollover shape of swing leg in [x,y] coordinates
    xth2roll_local=(xth2roll-xc)*cos(th2)+(yth2roll-yc).*sin(th2)+xc;        %same rollover shape as swing leg
    yth2roll_local=-(xth2roll-xc)*sin(th2)+(yth2roll-yc).*cos(th2)+yc;
    
    % rollover shape of swing leg in [X,Y] coordinates
    Xrlsw=(xth2roll_local-xth)*cos(th1)+(yth2roll_local-yth)*sin(th1)+sth;
    Yrlsw=-(xth2roll_local-xth)*sin(th1)+(yth2roll_local-yth)*cos(th1);
    plot((Xrlsw-centre)*cos(slope)+Yrlsw*sin(slope),(Xrlsw-centre)*-sin(slope)+Yrlsw*cos(slope))       % NEEDS TO BE CHANGED WITH DIFFERENT SHAPES
%% PLOT TOTAL ROLLOVER
    % Plots rollover shape in [X,Y] coordinates
    Xroll_total=(xthroll_total-xth)*cos(th1)+(ythroll_total-yth)*sin(th1)+sth;
    Yroll_total=-(xthroll_total-xth)*sin(th1)+(ythroll_total-yth)*cos(th1);
    plot((Xroll_total-centre)*cos(slope)+Yroll_total*sin(slope),(Xroll_total-centre)*-sin(slope)+Yroll_total*cos(slope),'b','LineWidth',1)
    
    % rollover shape of swing leg in [x,y] coordinates
    xth2roll_total_local=(xth2roll_total-xc)*cos(th2)+(yth2roll_total-yc).*sin(th2)+xc;        %same rollover shape as swing leg
    yth2roll_total_local=-(xth2roll_total-xc)*sin(th2)+(yth2roll_total-yc).*cos(th2)+yc;
    
    % rollover shape of swing leg in [X,Y] coordinates
    Xrlsw_total=(xth2roll_total_local-xth)*cos(th1)+(yth2roll_total_local-yth)*sin(th1)+sth;
    Yrlsw_total=-(xth2roll_total_local-xth)*sin(th1)+(yth2roll_total_local-yth)*cos(th1);
    plot((Xrlsw_total-centre)*cos(slope)+Yrlsw_total*sin(slope),(Xrlsw_total-centre)*-sin(slope)+Yrlsw_total*cos(slope),'b','LineWidth',1)
%% PLOT FLOOR
    plot((sth-centre)*cos(slope),(sth-centre)*-sin(slope),[stance,'.'])       % Plot contact point with the floor
    plot(linspace(-centre-1,-centre+1.5,100)*cos(slope),linspace(-centre-1,-centre+1.5,100)*-sin(slope),'k--')         % Plot floor    
    % Point of origin [0,0] is [-Xc(i),0]
    
    set(gca,'Xtick',[],'Ytick',[]); %get rid of axis ticks
    drawnow         % Updates matlab plot
    %plot(X1_0(i)-Xc(i),Y1_0(i),'ch')   % Plots ankle of stance leg
   
    %pause(0.05)
     M(frame) = getframe;
     frame=frame+1;
    hold off
%     figure(17)
%     plot(xth2,yth2)
end
% figure(2)
% movie(M,5,10)
% toc

leftover_time = t(end) - time;
%movie(M)
end
