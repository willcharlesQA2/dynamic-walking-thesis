function [  ] = limitCycles( filename,slope,c1,c2 )
%Takes data stored from birfurcations and plot a series of points
% Angle:          4   5.4   5.7   5.78    5.81
% Bifurcations:   1   2     4     8       16
% filename - Name of the file eg. 'SymmSlope.mat'
% slope - slope angle of required phase-plane doagram and sequence graph
% n - number of legs. 

if nargin < 3 
    c1 = 'k.';
    c2 = 'k.';
elseif nargin < 4
    c2 = c1;
end

load(filename)

% AllVar(bifurcation,:) = {alpha*-180/pi, uAll, DistAll, tf, inter, angvel};

x_axis = cell2mat(AllVar(:,1));

% The whole thing
allDist = cell2mat(AllVar(:,3));     % Step length
alltf = cell2mat(AllVar(:,4));          % step period
allinter = cell2mat(AllVar(:,5));       % inter-leg angle
allSpeed = allDist./alltf;

alluAll = AllVar(:,2);

%% Bifurcations spotted in 'symmetric_point2.mat' file are
% Angle:          4   5.4   5.7   5.78    5.81
% Bifurcations:   1   2     4     8       16
if nargin == 1
    slope = 4; % in degrees
end
bifurcations = 32;
% Find what position in the array we can find the required slope angle
indice = find(abs(x_axis-slope) < 0.001);
slopeinter = allinter(indice,:);

%% Could put these as an input to the function
Xlabel = 'Slope angle (deg)';

%% Plot bifurcation diagram
% Choose angles of 4, 5.4, 5.7, 5.78, 5.81
figure
hold on
h(:,1) = plot(x_axis,alltf(:,end-30:2:end),c1); % intact even
h(:,2) = plot(x_axis,alltf(:,end-31:2:end-1),c2); % prosthetic odd
xlabel(Xlabel,'FontSize',15)
ylabel('Stance time (s)','FontSize',15)
set(gca,'FontSize',12)
box on


%% Plot step diagram
figure
plot(slopeinter,'-*','MarkerEdgeColor','k','MarkerSize',3)
title(['Stable sequence graph at ',num2str(slope),' deg'])
xlabel('Step Number','FontSize',15)
ylabel('Inter-leg angle (rad)','FontSize',15)
set(gca,'FontSize',12)

%% Phase plane diagram
% u for the last step of the specified angle
maxsteps = size(alltf,2);
bif = 0;
for i = maxsteps+1-bifurcations:maxsteps
    bif = bif + 1;
    uSlope = cell2mat(alluAll{indice}(i,2));

    th1 = uSlope(:,1);
    phi = uSlope(:,2) + th1;
    th1dot = uSlope(:,3);
    phidot = uSlope(:,4) + th1dot;
    
    if mod(i,2) %odd
        if exist('leg1')
            leg1 = [leg1;th1];
            leg1dot = [leg1dot;th1dot];
            leg2 = [leg2;phi];
            leg2dot = [leg2dot;phidot];
        else
            leg1 = [th1];
            leg1dot = [th1dot];
            leg2 = [phi];
            leg2dot = [phidot];            
        end 
    else        %even
        if exist('leg1')
            leg1 = [leg1;phi];
            leg1dot = [leg1dot;phidot];  
            leg2 = [leg2;th1];
            leg2dot = [leg2dot;th1dot]; 
        else
            leg1 = [phi];
            leg1dot = [phidot];  
            leg2 = [th1];
            leg2dot = [th1dot]; 
        end
    end
end

%% Plot phase plane to show limit cycles
figure
hold on
plot(leg1,leg1dot)
plot(leg2,leg2dot,'r')
title(['Phase-plane diagram of point foot model at ',num2str(slope),' deg'])
xlabel('Leg angle (rad)','FontSize',15)
ylabel('Leg angular velocity (rad/s)','FontSize',15)
set(gca,'FontSize',12)
box on

%% Labels and legends
hl(1) = plot(NaN,NaN,c1);
hl(2) = plot(NaN,NaN,c2);

asymmetric = 0;     % if there are 2 different legs, set this as 1
if asymmetric == 1
    legend(hl(1,:),'unaffected leg','modified leg')
end
if mod(size(alltf,2),2)
    warning('There are an odd number of ''maxsteps''Check that even steps are prosthetic.')
end



end

