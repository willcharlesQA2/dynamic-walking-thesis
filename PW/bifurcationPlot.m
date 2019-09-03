function [  ] = bifurcationPlotNouAll( filename,c1,c2 )
%Takes data stored from birfurcations and plot a series of points
% filename - Name of the file eg. 'SymmSlope.mat'
% n - number of legs. 

if nargin == 1
    c1 = 'r.';
    c2 = 'b.';
elseif nargin < 3
    c2 = c1;
end

load(filename)

% AllVar(bifurcation,:) = {alpha*-180/pi, DistAll, tf, inter, angvel};
% or AllVar(bifurcation,:) = {r, DistAll, tf, inter, angvel};

x_axis = cell2mat(AllVar(:,1));

% The whole thing
allDist = cell2mat(AllVar(:,2));     % Step length
alltf = cell2mat(AllVar(:,3));          % step period
allinter = cell2mat(AllVar(:,4));       % inter-leg angle
allthdotS = cell2mat(AllVar(:,5));      % support leg velocity beginning of step
MEall = cell2mat(AllVar(:,6));          % Mechanical energy during step
diffE = cell2mat(AllVar(:,7));          % Energy lost at end of step
allthdotS = cell2mat(AllVar(:,5)); 
allSpeed = allDist./alltf;

ALLplots = {alltf;allinter;allSpeed;allthdotS;MEall;-diffE./MEall*100};


%% Could put these as an input to the function
% Xlabel = 'Length ratio of leg B, \beta';
% Xlabel = 'Mass ratio of leg B, \mu';
% Xlabel = 'Radius of curvature of the foot (\times leg length)';
% Xlabel = 'Roller radius of leg B, R_B (m)';
% Xlabel = 'Slope angle (deg)';
asymmetric = 0;     % if there are 2 different legs, set this as 1
average = 0;

Ylabels = {'Stance time (s)';'Inter-leg angle (rad)';'Walking speed (m/s)';...
    'Angular velocity of support leg (rad/s)';'Mechanical energy during step (J)';'% energy lost at collision'};

% Sym = 3.6;

%% Labels and legends
% hl = plot(NaN,NaN,'b-',NaN,NaN,'r-');

bifurcations = 16;

for i = 1:size(Ylabels,1)
    figure(i)
    hold on
        hl(:,1) = plot(NaN,NaN,'b-','linewidth',2);
        hl(:,2) = plot(NaN,NaN,'r-','linewidth',2);
        hl(:,3) = plot(NaN,NaN,'k--');
    plot(x_axis,ALLplots{i}(:,end-(bifurcations-2):2:end),c1); % B even
    plot(x_axis,ALLplots{i}(:,end-(bifurcations-1):2:end-1),c2); % A odd
    xlabel(Xlabel,'FontSize',18)
    ylabel(Ylabels(i),'FontSize',18)
    if asymmetric
        legend(hl(1,:),{'Leg A - unaffected','Leg B - modified','Mean of both'},'FontSize',14)
    end
    if average
        plot(x_axis,mean(ALLplots{i}(:,end-31:1:end-1),2),'k--')
%         legend(hl(1,:),'unaffected leg','modified leg','mean of both')
    end
    set(gca,'FontSize',13)
    box on
%     axis manual
    plot([Sym,Sym],ylim,'color',[0,0,0]+0.5)
end

if mod(size(alltf,2),2)
    warning('There are an odd number of ''maxsteps''Check that even steps are prosthetic.')
end

%% Point foot
% hl(1) = plot(NaN,NaN,'k-','linewidth',2);
% hl(2) = plot(NaN,NaN,'b-','linewidth',2);
% legend(hl(1,:),'Point foot','Curved foot')

end

