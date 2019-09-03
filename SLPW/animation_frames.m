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

global Param
alpha= Param.alpha;
S1th = Param.S1th;
S4th = Param.S4th;
Lr1 = Param.Lr1;
LH1 = Param.LH1;
h = Param.h;

xd2 = Param.xd2;
yd2 = Param.yd2;

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
    [xthroll(i),ythroll(i),~,~,~,~]=xth_yth(u(i,1),1,0,Param);      % USED rollover shape of stance foot
    [xth2roll(i),yth2roll(i),~,~,~,~]=xth_yth(u(i,1),2,0,Param);    % USED rollover shape of swing foot
end

count = 1;
for i = linspace(S1th,S4th,50)
    [xthroll_total(count),ythroll_total(count),~,~,~,~]=xth_yth(i,1,0,Param);      % TOTAL rollover shape of stance foot
    [xth2roll_total(count),yth2roll_total(count),~,~,~,~]=xth_yth(i,2,0,Param);    % TOTAL rollover shape of swing foot
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
 
for time=(delta_t-tstart)+t(1):delta_t:t(end) % (delta_t-tstart):delta_t:t(end)
    
th1 = interp1(t,u(:,1),time);
th2 = interp1(t,u(:,2),time);
r1 = interp1(t,u(:,3),time);
    if size(u,2) == 8
        % double support
        r2 = interp1(t,u(:,4),time);
    else
        %single support
        r2 = 0;
    end
    
%     q1 = u(i,1);
%     q2 = u(i,2);
%     q3 = u(i,3);
[xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(th1,1,0,Param);
[xth2,yth2,~,~,~,~]=xth_yth(th1,2,0,Param);      %leg2    
    
    
%% arc length
[ sth ,~ ,~ ] = arcLength(th1,dxth,dyth,ddxth,ddyth,S1th,S4th,h,1,Param);
    

    if mod(step,2)              % Physiological
        swing = 'r';
        stance = 'g';
    else                        % Prosthetic
        swing = 'g';
        stance = 'r';
    end

    %% Lengths of each leg
    yr1 = Lr1 + r1;
    yH1 = LH1 + r1;
    yH2 = LH1 + r2;
    xH2 = 0;
    
    % Hip mass locations
    Xc = -xth*cos(th1) + (yH1 - yth)*sin(th1) + sth;
    Yc = xth*sin(th1) + (yH1 - yth)*cos(th1);
    
    % Stance leg in SS. Trailing leg in DS.
    X1 = -xth*cos(th1) + (yr1 - yth)*sin(th1) + sth;
    Y1 = xth*sin(th1) + (yr1 - yth)*cos(th1);

    % Swing leg in [x1,y1] coordinates
    xp2 =    (xd2)*cos(th2)+(yd2)*sin(th2);
    yp2 =    -(xd2)*sin(th2)+(yd2)*cos(th2)+yH1;
    
    % Swing leg in SS. Leading leg in DS.
    X2 = (xp2-xth)*cos(th1) + (yp2-yth)*sin(th1) + sth;
    Y2 = -(xp2-xth)*sin(th1) + (yp2-yth)*cos(th1);

    % Ankle joints
    X1_0=(-xth)*cos(th1)+(-yth)*sin(th1)+sth; % [x,y] at [0,0]
    Y1_0=-(-xth)*sin(th1)+(-yth)*cos(th1);
    
    xr0=(-yH2)*sin(th2); %[x',y'] at [0,0]
    yr0=(-yH2)*cos(th2)+yH1;
    
    Xsw0=(xr0-xth)*cos(th1)+(yr0-yth)*sin(th1)+sth; % [X,Y] locations at end of swing leg %
    Ysw0=-(xr0-xth)*sin(th1)+(yr0-yth)*cos(th1);     % USED TO PLOT LEGS

centre = Xc*cos(alpha) - Yc*sin(alpha);

    %disp(t(i))
    plot([Xsw0-centre,Xc-centre]*cos(slope)+[Ysw0,Yc]*sin(slope),[Xsw0-centre,Xc-centre]*-sin(slope)+[Ysw0,Yc]*cos(slope),swing,'LineWidth',2)       % swing leg
    hold on
    %grid on
    axis([-0.7 0.7 -0.2 1.2])
    
    %% old 
        plot([X1_0-centre,Xc-centre]*cos(slope)+[Y1_0,Yc]*sin(slope),[X1_0-centre,Xc-centre]*-sin(slope)+[Y1_0,Yc]*cos(slope),stance,'LineWidth',2)       % stance leg
%     %%new Spring starts here
%     plot([X1_0-centre,X1-centre]*cos(slope)+[Y1_0,Y1]*sin(slope),[X1_0-centre,X1-centre]*-sin(slope)+[Y1_0,Y1]*cos(slope),stance,'LineWidth',2)       % stance leg
%     
%     % Plot spring for leg
%     X0s = [Xc,X1];
%     Y0s = [Yc,Y1];
%     L0 = Param.LH1 - Param.Lr1;
%     [ Xs,Ys ] = plotSpring( X0s,Y0s,L0 );
%     plot((Xs-centre)*cos(slope)+Ys*sin(slope),(Xs-centre)*-sin(slope)+Ys*cos(slope),stance,'LineWidth',2)       % stance leg
%     
%     %% Spring ends here
    
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
    xth2roll_local=(xth2roll)*cos(th2)+(yth2roll-yH2).*sin(th2);        %same rollover shape as swing leg
    yth2roll_local=-(xth2roll)*sin(th2)+(yth2roll-yH2).*cos(th2)+yH1;
    
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
    xth2roll_total_local=(xth2roll_total)*cos(th2)+(yth2roll_total-yH2).*sin(th2);        %same rollover shape as swing leg
    yth2roll_total_local=-(xth2roll_total)*sin(th2)+(yth2roll_total-yH2).*cos(th2)+yH1;
    
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
