load('resultskBv1.mat')

%% Parameters for graph
LW = 2;


%% Intact leg (odd step number)
figure; hold on

[TP,FP,pks,locs] = GRFplot2(AllGRFs(2,:));
plot(TP{2}/TP{2}(end)*100,FP{2}/(9.81*80),'LineWidth',LW)

[TP,FP,pks,locs] = GRFplot2(AllGRFs(1,:));
plot(TP{2}/TP{2}(end)*100,FP{2}/(9.81*80),'LineWidth',LW)

[TP,FP,pks,locs] = GRFplot2(AllGRFs(3,:));
plot(TP{2}/TP{2}(end)*100,FP{2}/(9.81*80),'LineWidth',LW)



title('Intact leg')
xlabel('Percentage Gait Cycle')
ylabel('Y Force / Body Weight')
legend('normal','5% compliant','5% stiffer')

%% Intact leg (even step number)
figure; hold on

[TP,FP,pks,locs] = GRFplot2(AllGRFs(2,:));
plot(TP{1}/TP{1}(end)*100,FP{1}/(9.81*80),'LineWidth',LW)

[TP,FP,pks,locs] = GRFplot2(AllGRFs(1,:));
plot(TP{1}/TP{1}(end)*100,FP{1}/(9.81*80),'LineWidth',LW)

[TP,FP,pks,locs] = GRFplot2(AllGRFs(3,:));
plot(TP{1}/TP{1}(end)*100,FP{1}/(9.81*80),'LineWidth',LW)

title('Residual leg')
xlabel('Percentage Gait Cycle')
ylabel('Y Force / Body Weight')
legend('normal','5% compliant','5% stiffer')