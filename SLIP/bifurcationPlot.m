function [  ] = bifurcationPlot( filename,c1,c2 )
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



% tf
% ws
% inter-leg
% GRFPeak1
% GRFPeak2
BW = 9.81*Param.m;
ALLplots = {StanceT;-interleg;ws;XNetImpulse;GRFPeak1/BW;GRFPeak2/BW};


%% Could put these as an input to the function
% Xlabel = 'Mechanical Energy (J)';

% Xlabel = 'Angle of attack of leg B, \alpha_B (deg)';
% Xlabel = 'Radius of curvature of foot B (\times leg length)';
% x_axis = alpha_total*-180/pi;

% Xlabel = 'Effective leg stiffness of leg B (kN/m)';
% x_axis = k_total/1000;
% Xlabel = 'Radius of curvature of foot B (\times leg length)';

% Sym = 0.6;

asymmetric = 0;     % if there are 2 different legs, set this as 1
average = 1;

% Ylabels = {'Stance time (s)';'Inter-leg angle (rad)';'Walking speed (m/s)';...
%     'Angular velocity of support leg (rad/s)';'Mechanical energy during step (J)';'% energy lost at collision'};
Ylabels = {'Stance time (s)';'Inter-leg angle (rad)';'Walking speed (m/s)';'Net X-impulse (Ns)';'FY1 / Body Weight';'FY2 / Body Weight'};

%% Labels and legends
% hl = plot(NaN,NaN,'b-',NaN,NaN,'r-');

FS = 18;
for i = 1:size(Ylabels,1)
    figure(i)
    hold on
        hl(:,1) = plot(NaN,NaN,'b-','linewidth',2);
        hl(:,2) = plot(NaN,NaN,'r-','linewidth',2);
        hl(:,3) = plot(NaN,NaN,'k--');
    plot(x_axis,ALLplots{i}(:,2:2:end),c1); % B even
    plot(x_axis,ALLplots{i}(:,1:2:end),c2); % A odd
    xlabel(Xlabel,'FontSize',FS)
    ylabel(Ylabels(i),'FontSize',FS)
    if asymmetric
        legend(hl(1,:),'unaffected leg','modified leg')
    end
    if average
        plot(x_axis,mean(ALLplots{i},2),'k--')
%         legend(hl(1,:),'unaffected leg','modified leg','mean of both')
    end
    set(gca,'FontSize',13)
    box on
    axis manual
    plot([Sym,Sym],ylim,'color',[0,0,0]+0.5)
end

% if mod(size(alltf,2),2)
%     warning('There are an odd number of ''maxsteps''Check that even steps are prosthetic.')
% end

%% Point foot
% hl(1) = plot(NaN,NaN,'k-','linewidth',2);
% hl(2) = plot(NaN,NaN,'b-','linewidth',2);
% legend(hl(1,:),'Point foot','Curved foot')

end