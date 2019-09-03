%% Run the code and plot GRFs for one asymmetric value
for loopF = 3:Bifs
    tend = -t_all{loopF-1}.t2(1); % last t value
    TP1 = tend+t_all{loopF-1}.t2;
    tend = tend+t_all{loopF-1}.t2(end);
    T{loopF-2} = [TP1(1:end-1);tend+t_all{loopF}.t1(1:end-1);tend+t_all{loopF}.t2];
    FPY{loopF-2} = [GRF_all{loopF-1}.FY2L(1:end-1);GRF_all{loopF}.FY1(1:end-1);GRF_all{loopF}.FY2T];            
end

Start = 7;
% Gait time
TPS1 = T{Start};
% Force
FPSy1 = FPY{Start};

TPS2 = T{Start+1};
FPSy2 = FPY{Start+1};

% Percentage of gait
PG1 = 100/TPS1(end);
PG2 = 100/TPS2(end);

% multiply by F to obtain percentage of body weight
PBW = 1/(Param.m*9.81);

% Line width
LW = 2;
% Font size
FS = 18;

load('normalGRFs.mat')
figure
hold on
plot(TPS1*PG1,FPSy1*PBW,'b','linewidth',LW)
plot(TPS2*PG2,FPSy2*PBW,'r','linewidth',LW)
plot(TPN*PGN,FPNy*PBW,'k--','linewidth',LW)
box on
xlabel('% gait cycle','FontSize',FS)
ylabel('GRF / Body Weight','FontSize',FS)

set(gca,'FontSize',13)
box on