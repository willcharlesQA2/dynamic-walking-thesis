function [  ] = bifurcationPlotPerc( filename,percPoint )
%Takes data stored from birfurcations and plot a series of points
% filename - Name of the file eg. 'SymmSlope.mat'
% n - number of legs. 


c1 = 'r.';
c2 = 'b.';


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
% Xlabel = 'Length ratio of leg B, \beta_B';
% Xlabel = 'Mass ratio of leg B, \mu_B';
% Xlabel = 'Radius of curvature of the foot (\times leg length)';
% Xlabel = 'Roller radius of leg B, R_B (m)';
% Xlabel = 'Slope angle (deg)';
asymmetric = 0;     % if there are 2 different legs, set this as 1
average = 1;

Ylabels = {'Stance time (s)';'Inter-leg angle (rad)';'Walking speed (m/s)';...
    'Angular velocity of support leg (rad/s)';'Mechanical energy during step (J)';'% energy lost at collision'};

% Sym = 3.6;

%% Labels and legends
% hl = plot(NaN,NaN,'b-',NaN,NaN,'r-');


for i = 1:size(Ylabels,1)
    figure(i)
    hold on
        hl(:,1) = plot(NaN,NaN,'b-','linewidth',2);
        hl(:,2) = plot(NaN,NaN,'r-','linewidth',2);
        hl(:,3) = plot(NaN,NaN,'k--');
    plot(x_axis,ALLplots{i}(:,end-30:2:end),c1); % B even
    plot(x_axis,ALLplots{i}(:,end-31:2:end-1),c2); % A odd
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
    axis manual
    plot([Sym,Sym],ylim,'color',[0,0,0]+0.5)
    
        %% Find percentage difference for percPoint
    % Symmetric values
    vSym(i) = interp1(x_axis,ALLplots{i}(:,end),Sym);
    % For each step, interpolate to find percPoint
    for j = 1:16
%         vSym(j) = interp1(x_axisbif,ALLplots{i}(j,:),Sym);
        vq(i,j) = interp1(x_axis,ALLplots{i}(:,end-16+j),percPoint);
    end
    
    %% Percentage difference
    % Odd steps
    pDiffOdd(i,:) = (vq(i,1:2:end)-vSym(i))/vSym(i)*100;
    
    % Even stops
    pDiffEven(i,:) = (vq(i,2:2:end)-vSym(i))/vSym(i)*100;
    
    % Print results if no /0 error
    if abs(pDiffEven(i,1)) < 1e5
        fprintf('%s P change: \tred = %.3f%%/%.3f%%\tblue = %.3f%%/%.3f%%\n',Ylabels{i},pDiffEven(i,1),pDiffEven(i,2),pDiffOdd(i,1),pDiffOdd(i,2))
    end
    
end

if mod(size(alltf,2),2)
    warning('There are an odd number of ''maxsteps''Check that even steps are prosthetic.')
end

%% Point foot
% hl(1) = plot(NaN,NaN,'k-','linewidth',2);
% hl(2) = plot(NaN,NaN,'b-','linewidth',2);
% legend(hl(1,:),'Point foot','Curved foot')

end

