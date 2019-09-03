function GRFplot(GRFinput,i)
%                 GRFs = {GRFy1,GRFy2,GRFx1,GRFx2,t1,t2,u1,u3};
GRFy1   = GRFinput{i,1};
GRFy2   = GRFinput{i,2};
GRFx1   = GRFinput{i,3};
GRFx2   = GRFinput{i,4};
t1      = GRFinput{i,5};
t2      = GRFinput{i,6};
u1      = GRFinput{i,7};
u3      = GRFinput{i,8};
Param   = GRFinput{i,9};

totalM = Param.m1 + Param.m2 + Param.m3;

colourRSuY = 'b';   %Rear Support
colourFSwY = 'b--'; %Front Swing
colourRSuX = 'r';   %Rear Support
colourFSwX = 'r--'; %Front Swing 

%% MUST BE SYMMETRIC WALKING
LW = 3;
figure
hold on
tSS = [t2-t2(1);t1-t2(1)];
tDS = t2-t2(1);
plot(tSS,GRFy1/(totalM*9.81)*100,colourRSuY,'LineWidth',LW)
% plot(tDS,GRFy2/(totalM*9.81)*100,colourFSwY,'LineWidth',LW)
plot(tDS+tSS(end),GRFy2/(totalM*9.81)*100,colourRSuY,'LineWidth',LW)
% plot(tSS+tSS(end),GRFy1/(totalM*9.81)*100,colourFSwY,'LineWidth',LW)

plot(tSS,GRFx1/(totalM*9.81)*100,colourRSuX,'LineWidth',LW)
% plot(tDS,GRFx2/(totalM*9.81)*100,colourFSwX,'LineWidth',LW)
plot(tDS+tSS(end),GRFx2/(totalM*9.81)*100,colourRSuX,'LineWidth',LW)
% plot(tDS+tSS(end),GRFx1/(totalM*9.81)*100,colourFSwY,'LineWidth',LW)
% Original:
% plot(t1-t2(1),F1ySS/(totalM*9.81)*100,colourRSuY)
% plot(t2-t2(1),F1y/(totalM*9.81)*100,colourFSwY)
% plot(t2-t2(1),F2y/(totalM*9.81)*100,colourRSuY)
% 
% plot(t1-t2(1),F1xSS/(totalM*9.81)*100,colourRSuX)
% plot(t2-t2(1),F1x/(totalM*9.81)*100,colourFSwX)
% plot(t2-t2(1),F2x/(totalM*9.81)*100,colourRSuX)
%% Grey double-support phase
DSt1 = t2(1)-t2(1);
SSt1 = t2(end)-t2(1);
DSt2 = t1(end)-t2(1);
SSt2 = t1(end) + t2(end) - 2*t2(1);
fill([DSt1,DSt1,SSt1,SSt1],[ylim fliplr(ylim)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3)
fill([DSt2,DSt2,SSt2,SSt2],[ylim fliplr(ylim)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3)

xlabel('Time (s)')
ylabel('GRF %BW')

end