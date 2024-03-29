clear all

% Parameters
%VALUES FOUND FROM GOSWAMI PG. 31
Param.Beta=0.6;         % length ratio (b/a)        % default Beta=1,   Mew=2
Param.Mu=3.6;          % mass ratio (mh/m)
M_total=80;
Param.L=1;

Beta = Param.Beta;
Mu = Param.Mu;

L = Param.L;

Param.mA = M_total/(2+Mu);
Param.mB = M_total/(2+Mu);
Param.m3 = Param.mA*Mu;

Param.LrA = L/(1+Beta);
Param.LrB = L/(1+Beta);
Param.LH1 = L;

Param.slope = 0;

Param.ME_0 = 710;
stiffness = 18e3;

Param.kA = stiffness; %9000 is original for both
Param.kB = stiffness;
Param.c1 = 0;
Param.c2 = 0;

% Foot parameters
%% Foot geometry
s1 = -0.4712;      % hindfoot length
s2 = 0;          %
s3 = 0;
s4 = 0.4712;       % forefoot length
%
% %% Point foot
% s1 = 0;
% s4 = 0;

% Foot gain
Param.frA = 0.3*Param.L;
Param.frB = 0.3*Param.L;

Param.h = 0; % Location of stance
% 
% [Param.S1th,Param.S2th,Param.S3th,Param.S4th,Param.xh,Param.yh,...
%     Param.xm,Param.ym,Param.xf,Param.yf]=rollshape(s1,s2,s3,s4,r,r,r,Param.h);

% arc length s = r*ang
Param.S1th = -pi/2;
Param.S2th = 0;
Param.S3th = 0;
Param.S4th =  pi/2;

Param.S1th2 = Param.S1th;
Param.S4th2 = Param.S4th;

%% Bifurcations
sizeA =  1;
sizeB = sizeA;
sizeC = sizeA;
sizeD = sizeA;
slope = 0*-pi/180;
slopes = linspace(0,0.1,sizeB)*-pi/180; % 0 -> 0.8


stiffness = 12e3; % try with 18 and time = 3
stiffnesses = linspace(15e3,22e3,sizeB); % 5 -> 25
% stiffnesses = [11.4,12,12.6]*1e3;

damping = 0;
dampings = linspace(0,5000,sizeB);

Param.maxsteps = 300;
maxsteps = Param.maxsteps;

speed = nan(11,sizeB);

%minimum ME (I think) is 667.08. Based on th1 = 0, th2 = 0 and KE = 0
% ME = 710;
MEs = linspace(680,710,sizeB);

%% Prelocate for speed
allSteps = nan(sizeA,sizeB);
uStableNew = NaN;

% Lower bounds and upper bounds
%   th1     r1      th1dot  r1dot
% Lb=[0.25     -0.11    1       0.1   ];
% Ub=[0.45     0        1.9       1   ];
Lb=[0  -0*pi/180];
Ub=[200 -5*pi/180  ];
bifurcation = 1;
tic
% Initial conditions for walking
%           q1        q3    q3dot
uStable = [ 0.302082,   -0.036527 ,   0.393525];
% uIni = [0.2767, -0.0587, 1.3, 0.34];
for aa = 1:1
    for bb = 1:1
        %                 uStable = [th1,r1,th1dot,r1dot];
        Param.frB = 0.26;
%         Param.kB = 20e3;
%         Beta = 0.57;
%         Mu = 3.8;


        Param.LrB = L/(1+Beta);
        Param.mB = Param.m3/Mu;
        
        Param.M_total = Param.mA + Param.mB + Param.m3;
        
        [ step, speed, GRFs, Param, results] = Walk(uStable,Param);
        
        if step ~= maxsteps
            [ uStableNew, step, speed, GRFs, Param, results] = findStableRAND(Param,200);
            
            if isnan(uStableNew) == 0
                uStable = uStableNew;
%                 [ step, speed, GRFs, Param] = Walk([uStable,stiffnesses(bb),slope,damping,ME]);
                
            end
        end
        
        allSteps(aa) = step;
        
        if step == maxsteps
            % Find Stable solution
            StableSolution(bifurcation,:) = {uStable,speed,MEs(aa),stiffnesses(bb),Param};
            Speeds(bifurcation,:) = speed;
            Allinterleg(:,bifurcation) = results.interleg;
            AllGRFs(bifurcation,:) = GRFs;
            [TP,FP,pks,locs,FPx] = GRFplot2(GRFs,16);
            bifurcation = bifurcation + 1;
            
            % Gait time
            TPS1 = TP{1};
            % Force
            FPSy1 = FP{1};
            
            TPS2 = TP{2};
            FPSy2 = FP{2};
            
            % Percentage of gait
            PG1 = 100/TPS1(end);
            PG2 = 100/TPS2(end);
            
            % multiply by F to obtain percentage of body weight
            PBW = 1/(Param.M_total*9.81);

            % Line width
            LW = 2;
            % Font size
            FS = 17;
            
            load('normalGRFs.mat')
            figure
            hold on
            plot(TPS1*PG1,FPSy1*PBW,'b','linewidth',LW)
            plot(TPS2*PG2,FPSy2*PBW,'r','linewidth',LW)
            plot(TPN*PGN,FPSyN*PBW,'k--','linewidth',LW)
            box on
            xlabel('% gait cycle','FontSize',FS)
            ylabel('GRF / Body Weight','FontSize',FS)

            set(gca,'FontSize',13)
            box on
            
            %% X-GRF
            FPSx1 = FPx{1};
            FPSx2 = FPx{2};
            
%             plot(TPS1*PG1,FPSx1*PBW,'b--','linewidth',LW)
%             plot(TPS2*PG2,FPSx2*PBW,'r--','linewidth',LW)
        end
    end
    fprintf('\nEstimate time left is %g hours.\n',toc*(sizeA/aa-1)/3600)
end
fprintf('\nTime complete with %g simuations done is %g hours. ',sizeA*sizeB,toc/3600)
% fprintf('\nEstimate time left is %g hours. ',toc*(sizeA/aa-1)/3600)
toc



% save('resultsTest.mat','StableSolution','allSteps','AllGRFs','slopes','stiffnesses')
% save('resultsNEW2.mat','StableSolution','Speeds','allSteps')

% figure
% FS = 17;
% plot(1:maxsteps,results.GCrdot,'-*','MarkerEdgeColor','k','MarkerSize',3)
% title('Stable sequence graph at ?')
% xlabel('Step Number','FontSize',FS)
% ylabel('rdot at mid-stance','FontSize',FS)
% set(gca,'FontSize',13)