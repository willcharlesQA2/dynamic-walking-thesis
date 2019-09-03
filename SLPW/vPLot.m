function vPLot(GRFinput,i)
% Used to compare displacements and velocities of DOFs with experimental
% data
%  get(gca,'ColorOrder')
t1      = GRFinput{i,5};
t2      = GRFinput{i,6};
u1      = GRFinput{i,7};
u3      = GRFinput{i,8};
Param   = GRFinput{i,9};

% Default colors
col1 = [       0    0.4470    0.7410  ];
col2 = [  0.8500    0.3250    0.0980  ];
col3 = [  0.9290    0.6940    0.1250  ];
col4 = [  0.4940    0.1840    0.5560  ];

%Line width
LW = 3;

v2 = u2v(u3);
v1 = u2v(u1);

T = [t2;t1];
% v Switched
v1s = [v1(:,2),v1(:,1),v1(:,3)];
V = [v2(:,1:2),v2(:,4);v1s(:,1:3)];

% figure
% title('Diaplacements')
% plot(T,V)
% legend('th1','th2','r1')
% hold on

%% Requires stable, symmetric walking for it to look nice
tAll  = [t2-t2(1);t1-t2(1);t2-2*t2(1)+t1(end);t1-2*t2(1)+t1(end)];
q1 = [v2(:,1);v1(:,2);v2(:,2);v1(:,1)];
q2 = [v2(:,2);v1(:,1);v2(:,1);v1(:,2)];

%

q3d1 = [v2(:,3)];
q3dd = [zeros(size(t1))];
q3d2 = [ v2(:,4);v1(:,3)];
q4d1 = [v2(:,4);v1(:,3);v2(:,3)];
q4dd = [zeros(size(t1))];

figure
hold on
plot(tAll,q1,'Color',col1,'LineWidth',LW)
plot(tAll,q2,'Color',col2,'LineWidth',LW)
plot(t2-t2(1),q3d1,'Color',col3,'LineWidth',LW)
plot([t2-t2(1);t1-t2(1);t2-2*t2(1)+t1(end)],q4d1,'Color',col4,'LineWidth',LW)
legend('th1','th2','r1','r2','LineWidth',LW)
plot(t1-t2(1),q3dd,'--','Color',col3,'LineWidth',LW)
plot([t2-2*t2(1)+t1(end);t1-2*t2(1)+t1(end)],q3d2,'Color',col3,'LineWidth',LW)
plot(t1-2*t2(1)+t1(end),q4dd,'--','Color',col4,'LineWidth',LW)

% Plot grey double support space
DSt1 = t2(1)-t2(1);
SSt1 = t2(end)-t2(1);
DSt2 = t1(end)-t2(1);
SSt2 = t1(end) + t2(end) - 2*t2(1);
% Stops axis from messing up
YBound = ylim;
fill([DSt1,DSt1,SSt1,SSt1],[YBound fliplr(YBound)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3)
fill([DSt2,DSt2,SSt2,SSt2],[YBound fliplr(YBound)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3)
