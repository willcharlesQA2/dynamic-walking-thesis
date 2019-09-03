function [X1,Y1,X2,Y2,X3,Y3,Xc,Yc,Xst,sth] = locations(t,u,Param,step,phase)

%Plots trajectories of masses - Make this into a function
% Obtain lengths from model parameters
if strcmp(phase,'SS')
    phaseMod = 1; % odd
elseif strcmp(phase,'DS')
    phaseMod = 0; % even
end

[ Param ] = switchLeg( Param,step,phaseMod );

Lr1 = Param.Lr1;
Lr2 = Param.Lr2;
LH1 = Param.LH1;
xd2 = 0;
yd2 = Param.Lr2-Param.LH1;

S1th = Param.S1th;
S4th = Param.S4th;
h = Param.h;

% initialise Xc and Yc
X1 = zeros(size(u,1),1);
Y1 = zeros(size(u,1),1);

X2 = zeros(size(u,1),1);
Y2 = zeros(size(u,1),1);

X3 = zeros(size(u,1),1);
Y3 = zeros(size(u,1),1);

Xc = zeros(size(u,1),1);
Yc = zeros(size(u,1),1);

sth = zeros(size(u,1),1);

for i = 1:size(t)
    
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
[xth2,yth2,dxth2,dyth2,ddxth2,ddyth2]=xth_yth(q1+q2,2,0,Param);

[ sth(i), ~,~ ] = arcLength(q1,1,1,Param);
    
    %% Lengths of each leg
    yr1 = Lr1 + q3;
    yH1 = LH1 + q3;
    yH2 = LH1 + q4;
    xH2 = 0;
    
    % Hip mass locations
    X3(i) = -xth*cos(q1) + (yH1 - yth)*sin(q1) + sth(i);
    Y3(i) = xth*sin(q1) + (yH1 - yth)*cos(q1);
    
    % Stance leg in SS. Trailing leg in DS.
    X1(i) = -xth*cos(q1) + (yr1 - yth)*sin(q1) + sth(i);
    Y1(i) = xth*sin(q1) + (yr1 - yth)*cos(q1);
    
    % Swing leg in [x1,y1] coordinates
    xp2 =    (xd2)*cos(q2)+(yd2)*sin(q2);
    yp2 =    -(xd2)*sin(q2)+(yd2)*cos(q2)+yH1;
    
    % Swing leg in SS. Leading leg in DS.
    X2(i) = (xp2-xth)*cos(q1) + (yp2-yth)*sin(q1) + sth(i);
    Y2(i) = -(xp2-xth)*sin(q1) + (yp2-yth)*cos(q1);
    
    % End of swing foot in [x1,y1] coordinates
    xc =  (xth2-xH2)*cos(q2) + (yth2-yH2)*sin(q2);
    yc = -(xth2-xH2)*sin(q2) + (yth2-yH2)*cos(q2) + yH1;
    
    % End of swing leg in SS. Foot contact of lead leg in DS.
    Xc(i) = (xc-xth)*cos(q1) + (yc-yth)*sin(q1) + sth(i);
    Yc(i) = -(xc-xth)*sin(q1) + (yc-yth)*cos(q1);

end

% fprintf('Xc(1) \t\t= %g\t\tYc(1) \t\t= %g\n',Xc(1),Yc(1))
% fprintf('Xc(end) \t= %g\t\tYc(end) \t= %g\n\n',Xc(end),Yc(end))
% 
% totalChange = ( (Xc(end)-Xc(1))^2 + (Yc(end)-Yc(1))^2 )^0.5;
% 
% fprintf('Total change is %gmm\n\n',totalChange*1000)

[ sth2(i), ~,~ ] = arcLength(q1+q2,1,2,Param);
Xst = Xc(i) - sth2(i);
% fprintf('Function ''locations'' used here. Check Xc is right\n')
end