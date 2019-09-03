% x = 0:.05:2*pi;
% y = sin(x);
% z = zeros(size(x));
% col = x;  % This is the color, vary with x in this case.
% surface([x;x],[y;y],[z;z],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);

% Support leg in Single stance
th1SS = results.v1(:,1);
th1dotSS = results.v1(:,4);
r1SS = results.v1(:,3);

%Rear leg in single-stance
th1DS = results.v2(:,1);
th1dotDS = results.v2(:,5);
r1DS = results.v2(:,3);

%Swing leg in swing stance
th2SS = results.v1(:,2);
th2dotSS = results.v1(:,5);
r2SS = zeros(size(results.v1(:,1)));

% Front leg in double-stance
th2DS = results.v2(:,2);
th2dotDS = results.v2(:,6);
r2DS = results.v2(:,4);

% %% 2D plot
% figure; hold on;
% plot(th1SS,th1dotSS)    % Support leg in Single stance
% plot(th1DS,th1dotDS)    %Rear leg in single-stance
% plot(th2SS,th2dotSS)    %Swing leg in swing stance
% plot(th2DS,th2dotDS)    % Front leg in double-stance

%% Surf plot
Allth = [th1SS;th1DS;th2SS;th2DS];
Allthdot = [th1dotSS;th1dotDS;th2dotSS;th2dotDS];
Allr = [r1SS;r1DS;r2SS;r2DS];
figure; hold on;
surface([Allth,Allth],[Allthdot,Allthdot],[-Allr*100,-Allr*100],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',3);
h = colorbar;
FS = 17;
ylabel(h, 'Compression of spring, r (cm)','FontSize',FS)
colormap(jet)
% title(['Phase-plane diagram of point foot model at ',num2str(slope),' deg'])
xlabel('Leg angle (rad)','FontSize',FS)
ylabel('Leg angular velocity (rad/s)','FontSize',FS)
set(gca,'FontSize',13)
box on

plot(th1SS(1),th1dotSS(1),'k+','MarkerSize',5)
plot(th1SS(end),th1dotSS(end),'k+','MarkerSize',5)
plot(th1DS(end),th1dotDS(end),'k+','MarkerSize',5)
plot(th2DS(1),th2dotDS(1),'k+','MarkerSize',5)
% surf(th1SS,th1dotSS,r1SS)    % Support leg in Single stance
% surf(th1DS,th1dotDS,r1DS)    %Rear leg in single-stance
% surf(th2SS,th2dotSS,r2SS)    %Swing leg in swing stance
% surf(th2DS,th2dotDS,r2DS)    % Front leg in double-stance
