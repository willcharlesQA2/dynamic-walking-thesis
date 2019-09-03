function [ uStable, stepRfinal, speed, GRFs, Param, results ] = findStableRAND(Param,n)
%Used to find stable point for bifurcation analysis
% 1. Take a random selection
% 2. Run code for N times
% 3. Find highest number reached
% 4. Take smaller subsection
% ...
% Repeat until steps = maxsteps

%%
% Number of nests
if nargin == 1
    n = 200;
end
% Fraction taken
pf = 0.2;

% Maximum number of steps required
maxSteps = Param.maxsteps;

% Lower bounds and upper bounds
%   th1     r1          r1dot
Lb=[0.2     -0.05        0.2   ];% 7000  -0.3*pi/180];
Ub=[0.4     -0.02        0.5   ];% 18000 -3*pi/180  ];

% Do you want to plot the stable points taken?
plotting = 0;

if plotting == 1
    figure(1)
    clf
    subplot(2,2,1)
    set(gca,'XLim',[Lb(1) Ub(1)],'YLim',[Lb(2) Ub(2)])
    hold on
    subplot(2,2,2)

    set(gca,'XLim',[Lb(3) Ub(3)],'YLim',[Lb(4) Ub(4)])
    hold on
end

% Number of times the function is repeated
for time = 1:1
for i=1:n,
    nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
    
    % Steps reached
    [stepR(i),speed,GRFs,Param, results] = Walk(nest(i,:),Param);
    if stepR(i) == maxSteps
        uStable = nest(i,:);
        stepRfinal = stepR(i);
        return
    end
    % Print out a message after each attempt
%     fprintf('n = %g\t%g steps complete\n',i,stepR(i))
    fprintf('-')
    
    if plotting == 1
            subplot(2,2,1)
            plot3(nest(i,1),nest(i,2),stepR(i),'o')
    %         set(gca,'XLim',[Lb(1) Ub(1)],'YLim',[Lb(2) Ub(2)])
            subplot(2,2,2)
            plot3(nest(i,3),nest(i,4),stepR(i),'o')
    %         set(gca,'XLim',[Lb(3) Ub(3)],'YLim',[Lb(4) Ub(4)])
    %         subplot(2,2,[3,4])
    %         plot(FPlot(1:G,1),'-+r')
            %set(gca,'YScale','log')
            refreshdata
            drawnow
    end
end

fprintf('\nPass %g complete.\n',i)
% %% Breed
% [newNests] = Breed (nest,stepR,3);
% clear stepR
% for i=1:size(newNests,1), 
%     % Steps reached
%     [stepR(i),speed] = Walk([newNests(i,:),input]);
%     if stepR(i) == maxSteps
%         uStable = nest(i,:);
%         return
%     end
%     % Print message after each breeding attempt
% %     fprintf('b = %g\t%g steps complete\n',i,stepR(i))
%     fprintf('-')
% end
% 
% % End of breeding
% fprintf('\nBreed %g complete.\n',i)

% Maximum steps reached
[maxStep(time),K] = max(stepR);

%% Plot steps taken
if plotting == 1
    averageStep(time) = mean(stepR);
    subplot(2,2,[3,4])
    % clf
    hold on
    plot(1:time,averageStep,'r+-')
    plot(1:time,maxStep,'k+-')
end
% Take a subsection
% diff = (Ub-Lb);
% Lb = newNests(K,:)-(pf*diff/2);
% Ub = newNests(K,:)+(pf*diff/2);

end
 
fprintf('No stable solution found.\n')
uStable = NaN;
stepRfinal = NaN;
speed = NaN;
GRFs = NaN;
results = NaN;
end

