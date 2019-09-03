function [M,leftover_time] = animation_frames(t,u,Param,Cond,step,tstart)
% Plot walker as a diagram

if nargin == 0
    % Original leg length
    L0 = 1;
    % Foot length
    fl = 0.2;
    fr = 0.3; % foot radius
    theta = -0.5;
    r = 0;
else
    L0 = Param.L0;  % Natural leg length
    fr = Param.fr1;  % foot radius
    Aa = Param.alpha;
end

% fl = L0*thetaEnd;
% Maximum Theta, when leg becomes locked
Th = Param.Th;
Tt = Param.Tt;

%% Plot walker
% Linewidth
LW = 2.5;
% Scaling for axis
scaling = 0.3;

figure(1)
set(gcf, 'Renderer', 'zbuffer');
frame=1;        % used for movie loop
fps = 20;       % frames per second (hz)
delta_t = 1 / fps;

if mod(step,2)              % Physiological
    swing = 'r';
    stance = 'b';
else                        % Prosthetic
    swing = 'b';
    stance = 'r';
end
    
for time=(delta_t-tstart)+t(1):delta_t:t(end)
    
    theta   = interp1(t,u(:,1),time);
    r       = interp1(t,u(:,2),time);
    thdot   = interp1(t,u(:,3),time);
    rdot    = interp1(t,u(:,4),time);
    

    
    L = L0 + r;

    
    % Whole foot
    xall =  fr*sin(Th:0.01:Tt);
    yall = -fr*cos(Th:0.01:Tt)+fr;
    
    % Locked rocker
    [ CoP,thC,~,~,~,~,~,~ ] = findCoP( theta,Param,1 );
    
    % CoP in local coordinates
    xc =  fr*sin(thC);
    yc = -fr*cos(thC) + fr;
    
    X =  (-xc)*cos(theta) + (L-yc)*sin(theta) + CoP;
    Y = -(-xc)*sin(theta) + (L-yc)*cos(theta);
    
    
    %Ankle
    Xank =  (-xc)*cos(theta) + (-yc)*sin(theta)+CoP;
    Yank = -(-xc)*sin(theta) + (-yc)*cos(theta);
    
    % Foot plot
    Xft =  (xall-xc)*cos(theta) + (yall-yc)*sin(theta)+CoP;
    Yft = -(xall-xc)*sin(theta) + (yall-yc)*cos(theta);
    



    
    clf
    hold on
        
    % Position of mass
    plot(X,Y,'ko','MarkerSize',15,'MarkerFaceColor','k')
    % Plot ankle
%     plot(Xank,Yank,'ko','MarkerSize',5,'MarkerFaceColor','k')
    % Plot virtual leg
%     plot([CoP,X],[0,Y],'r')
    % Plot leg
    plot([Xank,X],[Yank,Y],stance)
    % Plot foot
    plot(Xft,Yft,stance)
    % Plot COP
    plot(CoP,0,'k.','MarkerSize',2,'MarkerFaceColor','k')
    % Plot touch down leg
%     plot([X,Xtd],[Y,Ytd])

switch Cond
    case {'SS'}
        % Touch down leg
        Xtd = -L0*sin(Aa) + X;
        Ytd = -L0*cos(Aa) + Y;
        TD = [swing,'--'];
        plot([X,Xtd],[Y,Ytd],TD)
    case {'DS'}
%         Xtd = Param.IC;
%         Ytd = 0;
        TD = swing;
        
        [ th2,r2,th2dot,r2dot ] = newConditions( [theta,r,thdot,rdot],Param );
        
        phi = -theta + th2;
        
        [ CoP2,thC2,dCoP2,ddCoP2,dxc2,dyc2,ddxc2,ddyc2 ] = findCoP( th2,Param,2 );
        
        L2 = L0 + r2;
        
        xc2 =  Param.fr*sin(thC2);
        yc2 = -Param.fr*cos(thC2) + Param.fr;

            %% Swing leg
            % Mass position from theta_2
            X2 =  (-xc2)*cos(th2) + (L2-yc2)*sin(th2) + CoP2 + Param.IC;
            Y2 = -(-xc2)*sin(th2) + (L2-yc2)*cos(th2);

            xr0=(-L2)*sin(phi); %[x',y'] at [0,0]
            yr0=(-L2)*cos(phi)+L;

            Xsw0=(xr0-xc)*cos(theta)+(yr0-yc)*sin(theta)+CoP; % [X,Y] locations at end of swing leg %
            Ysw0=-(xr0-xc)*sin(theta)+(yr0-yc)*cos(theta);     % USED TO PLOT LEGS

        %% PLot leg
        plot([X2,Xsw0],[Y2,Ysw0],TD)
        
        % Plot roll-over for swing leg
        Xft2 =  (xall-xc2)*cos(th2) + (yall-yc2)*sin(th2)+CoP2 + Param.IC;
        Yft2 = -(xall-xc2)*sin(th2) + (yall-yc2)*cos(th2);
        
        plot(Xft2,Yft2,swing)
    otherwise
        error('Case not recognised')
end


   
    

    
    axis([X-L0/2-scaling,X+L0/2+scaling,-scaling,L0+scaling])
    axis square
    
    % Plot floor
    plot(xlim,[0 0],'k--')
    
%     set(gca,'Xtick',[],'Ytick',[]); %get rid of axis ticks
    drawnow         % Updates matlab plot
    %plot(X1_0(i)-Xc(i),Y1_0(i),'ch')   % Plots ankle of stance leg

    %pause(0.05)
    M(frame) = getframe;
    frame=frame+1;
    hold off
    
end

leftover_time = t(end) - time;

end

