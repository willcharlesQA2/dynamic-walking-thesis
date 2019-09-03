function [ TPS,FPSy ] = SlowGRF( alpha,k,ME, Param )

%% Input paramaters to run code again and plot GRF


[ t_all,GRF_all ] = GRFfromPS( alpha, k, ME, Param );

% figure
% hold on
% tend = -t_all{2-1}.t2(1); % last t value
% for i = 2:3; % number of steps
%     plot(tend+t_all{i-1}.t2,GRF_all{i-1}.FY2L)
%     tend = tend+t_all{i-1}.t2(end);
%     plot(tend+t_all{i}.t1,GRF_all{i}.FY1)
%     plot(tend+t_all{i}.t2,GRF_all{i}.FY2T)
% end
% 
% xlabel('time')
% ylabel('GRF')
% title(['\alpha = ',alpha*-180/pi,'deg; k = ',k/1000,'kN/m; ME = ',ME,'J'])


%% Find GRF of one leg

for i = 3; % number of steps
    tend1 = -t_all{i-1}.t2(1); % last t value so graph starts at 0
    tend2 = tend1+t_all{i-1}.t2(end);
    TPDS1 = tend1+t_all{i-1}.t2;
    FPDS1y = GRF_all{i-1}.FY2T;
    FPDS1x = GRF_all{i-1}.FX2T;
    TPDS2 = tend2+t_all{i}.t2;
    FPDS2y = GRF_all{i}.FY2L;
    FPDS2x = GRF_all{i}.FX2L;

    % Time of Stance
    TPS = [TPDS1(1:end-1);tend2+t_all{i}.t1(1:end-1);TPDS2];
    % Force of Stance leg
    FPSy = [GRF_all{i-1}.FY2L(1:end-1);GRF_all{i}.FY1(1:end-1);GRF_all{i}.FY2T];
    FPSx = [GRF_all{i-1}.FX2L(1:end-1);GRF_all{i}.FX1(1:end-1);GRF_all{i}.FX2T];
end
figure
hold on

xlabel('% Gait Cycle','FontSize',15)
ylabel('GRF / Body Weight','FontSize',15)
% title(['\alpha = ',num2str(alpha*-180/pi),'deg; k = ',num2str(k/1000),'kN/m; ME = ',num2str(ME),'J'],'FontSize',15)
% multiply by T to obtain percentage gat=it cycle
PG = 100/TPS(end);

% multiply by F to obtain percentage of body weight
PBW = 1/(Param.m*9.81);

% Line width
LW = 2;

plot(TPS*PG,FPSy*PBW,'b','linewidth',LW)
plot(TPS*PG,FPSx*PBW,'b--','linewidth',LW)
plot(TPDS1*PG,FPDS1y*PBW,'r','linewidth',LW)
plot(TPDS1*PG,FPDS1x*PBW,'r--','linewidth',LW)
plot(TPDS2*PG,FPDS2y*PBW,'r','linewidth',LW)
plot(TPDS2*PG,FPDS2x*PBW,'r--','linewidth',LW)

set(gca,'FontSize',12)
box on

figure
hold on
LW = 5;
plot(TPS*PG,FPSy*PBW,'k','linewidth',LW)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'Color','y')
box on
end