function [  ] = bifurcationPlotPerc( filename,percPoint)
%Takes data stored from birfurcations and plot a series of points
% filename - Name of the file eg. 'SymmSlope.mat'
% n - number of legs. 

c1 = 'r.';
c2 = 'b.';


load(filename)

% tf
% ws
% inter-leg
% GRFPeak1
% GRFPeak2

BW = ones(size(GRFPeak1,1),1)*AllM_total*9.81;
% BW = 9.81*(Param.mA+Param.mB+Param.m3);
StanceT = AllT2(1:end-1,:) + AllT1(1:end-1,:) + AllT2(2:end,:);

ALLplots = {StanceT;-Allinterleg;Speeds;AllGCrdot;GRFPeak1./BW;GRFPeak2./BW;XNetImpulse};


%% Could put these as an input to the function
% Xlabel = 'Angle of attack of leg B (deg)';
% x_axis = alpha_total*-180/pi;

% Xlabel = 'Mass ratio of leg B, \mu_B';

% Xlabel = 'Length ratio of leg B, \beta_B';

% x_axis = k_total/1000;
% Xlabel = 'Radius of curvature of the foot (\times leg length)';

% Xlabel = 'Stiffness of leg B, k_B (kN/m)';

asymmetric = 0;     % if there are 2 different legs, set this as 1
average = 1;

% Ylabels = {'Stance time (s)';'Inter-leg angle (rad)';'Walking speed (m/s)';...
%     'Angular velocity of support leg (rad/s)';'Mechanical energy during step (J)';'% energy lost at collision'};
Ylabels = {'Stance time (s)';'Inter-leg angle (rad)';'Walking speed (m/s)';'Spring velocity at mid-stance (m/s)';'FY1 / Body Weight';'FY2 / Body Weight';'Net X-impulse (Ns)'};

%% Labels and legends
% hl = plot(NaN,NaN,'b-',NaN,NaN,'r-');
% Sym = 18;

for i = 1:size(Ylabels,1)
    figure(i)
    hold on
        hl(:,1) = plot(NaN,NaN,'b-','linewidth',2);
        hl(:,2) = plot(NaN,NaN,'r-','linewidth',2);
        hl(:,3) = plot(NaN,NaN,'k--');
    plot(x_axisbif,ALLplots{i}(2:2:end,:),c1); % B even red
    plot(x_axisbif,ALLplots{i}(1:2:end,:),c2); % A odd blue
    xlabel(Xlabel,'FontSize',17)
    ylabel(Ylabels(i),'FontSize',17)
    if asymmetric
        legend(hl(1,:),'unaffected leg','modified leg')
    end
    if average
        plot(x_axisbif,mean(ALLplots{i}),'k--')
%         legend(hl(1,:),'unaffected leg','modified leg','mean of both')
    end
    set(gca,'FontSize',13)
    box on
    axis manual
    plot([Sym,Sym],ylim,'color',[0,0,0]+0.5)
    
    %% Find percentage difference for percPoint
    % Symmetric values
    vSym(i) = interp1(x_axisbif,ALLplots{i}(1,:),Sym);
    % For each step, interpolate to find percPoint
    for j = 1:14
%         vSym(j) = interp1(x_axisbif,ALLplots{i}(j,:),Sym);
        vq(i,j) = interp1(x_axisbif,ALLplots{i}(j,:),percPoint);
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



end

