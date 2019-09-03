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

Param.ME_0 = 710;
Param.slope = 0;


Param.kA = 18000; %9000 is original for both
Param.kB = 18000;
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
sizeA =  301; % was 301
sizeB = sizeA;
sizeC = sizeA;
sizeD = sizeA;
slope = 0*-pi/180;
slopes = linspace(0,0.1,sizeB)*-pi/180; % 0 -> 0.8


stiffness = 12e3; % try with 18 and time = 3
stiffnesses = linspace(15e3,22e3,sizeB); % 5 -> 25
% stiffnesses = [11.4,12,12.6]*1e3;

Mus = linspace(1.4,4,sizeB);
% Betas = linspace(0.55,0.75,sizeB);
Betas = linspace(0.4,0.8,sizeB);

damping = 0;
dampings = linspace(0,5000,sizeB);

Param.maxsteps = 200;
maxsteps = Param.maxsteps;

speed = nan(11,sizeB);

%minimum ME (I think) is 667.08. Based on th1 = 0, th2 = 0 and KE = 0
ME = 710;
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
uStable = [ 0.3288,   -0.0126 ,   0.3891];
% uIni = [0.2767, -0.0587, 1.3, 0.34];
for aa = 1:1
    for bb = 1:sizeB
%         Param.kB = stiffnesses(bb);
%         Param.mB = Param.m3/Mus(bb);
        Param.LrA = L/(1+Betas(bb));
        Param.LrB = L/(1+Betas(bb));
        x_axis = Betas(bb);
        
        Param.M_total = Param.mA + Param.mB + Param.m3;
        
        [ step, speed, GRFs, Param, results] = Walk(uStable,Param);
        
        if step ~= maxsteps
            [ uStableNew, step, speed, GRFs, Param, results] = findStableRAND(Param,100);
            
            if isnan(uStableNew) == 0
                uStable = uStableNew;
%                 [ step, speed, GRFs, Param] = Walk([uStable,stiffnesses(bb),slope,damping,ME]);
                
            end
        end
        
        allSteps(aa) = step;
        
        if step == maxsteps
            %% Compile results
            % Found the Stable solution
            StableSolution(:,bifurcation) = {uStable,speed,MEs(aa),stiffnesses(bb),Param};
            Speeds(:,bifurcation) = speed;
            AllGRFs(:,bifurcation) = GRFs;
            
            % Change this from bifurcation in future. CHange to bb
            x_axisbif(:,bifurcation) = x_axis;
            AllT1(:,bifurcation) = results.T1;
            AllT2(:,bifurcation) = results.T2;
            Allth1dot(:,bifurcation) = results.th1dot;
            Allinterleg(:,bifurcation) = results.interleg;
            AllGCth1dot(:,bifurcation) = results.GCth1dot;
            AllGCr(:,bifurcation) = results.GCr;
            AllGCrdot(:,bifurcation) = results.GCrdot;
            AllM_total(:,bifurcation) = Param.M_total;
            AllAlpha(:,bifurcation) = results.touchDownAngle;
            
            %% GRFs
            [TP,FP,pks,locs,FPx] = GRFplot2(AllGRFs(:,bifurcation),16);
            for imp = 1:15
                XNetImpulse(imp,bifurcation) = trapz(TP{imp},FPx{imp});
            end
            %% To get pks, check they are all the same length
            Allpks = [pks{:}];
            size_pks(:,bifurcation) = size(Allpks,1);
            GRFPeak1(:,bifurcation) = Allpks(1,:);
            GRFPeak2(:,bifurcation) = Allpks(2,:);
            
            %% One value for each stable solution
            bifurcation = bifurcation + 1;
        end
    end
    fprintf('\nEstimate time left is %g hours.\n',toc*(sizeA/aa-1)/3600)
end
fprintf('\nTime complete with %g simuations done is %g hours. ',sizeA*sizeB,toc/3600)
% fprintf('\nEstimate time left is %g hours. ',toc*(sizeA/aa-1)/3600)
toc

% save('dataSA.mat','StableSolution','Speeds')
% figure
% for  i = 1:size(StableSolution,1)
%     plot3(StableSolution{i,4},StableSolution{i,3},StableSolution{i,2}(8:end),'bo')
%     hold on
% end
% grid on
% xlabel('Stiffness')
% ylabel('Energy')
% zlabel('Walking Speed (m/s)')

Xlabel = 'Length ratio of legs, \beta';
Sym = 0.6;
save('resultsBifBetaBothV1.mat')
% save('resultsNEW2.mat','StableSolution','Speeds','allSteps')

% figure
% FS = 17;
% plot(1:maxsteps,results.GCrdot,'-*','MarkerEdgeColor','k','MarkerSize',3)
% title('Stable sequence graph at ?')
% xlabel('Step Number','FontSize',FS)
% ylabel('rdot at mid-stance','FontSize',FS)
% set(gca,'FontSize',13)