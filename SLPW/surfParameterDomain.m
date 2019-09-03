% load('resultsParameterExperimentSwoma3.mat')
% load('resultsParameterAlpha.mat')
colourPlot = 1;

k_total = stiffnesses;

ME_total = MEs;

for i = 1:15
    rdotDiff(:,:,i) = AllrdotStart(:,:,i+1)-AllrdotEnd(:,:,i);
%     surf(k_total*1e3,ME_total,rdotDiff(:,:,i),'FaceColor',[0 1 0],'FaceAlpha',0.10)
end

Colours = {[1 0 0],[0 1 0],[0 0 1]};
zlabels = {'Alpha (deg)','Change in spring velocity rdot','Walking Speed (m/s)','Number of vGRF peaks'};
zParam = {AllAlpha*-180/pi,rdotDiff,AllSpeeds,BOApks};
iSize = {16,15,16,1};

for plotLoop = 1:4
    figure; hold on
    for i = 1:iSize{plotLoop}
        if colourPlot == 1
            surf(k_total*1e-3,ME_total,zParam{plotLoop}(:,:,i)) %,'FaceAlpha',1/iSize{plotLoop}
            colormap(jet)
        else
            surf(k_total*1e3,ME_total,zParam{plotLoop}(:,:,i),'FaceColor',Colours{plotLoop},'FaceAlpha',0.10)
        end
    end

    FS = 17;

    % xlabel('Angle of attack, \alpha, deg')
    xlabel('Leg stiffness, k (kN/m)','FontSize',FS)
    ylabel('Mechanical Energy (J)','FontSize',FS)
    zlabel(zlabels{plotLoop},'FontSize',FS)

    view(0,90)
    set(gca,'FontSize',13)
    box on
end

% %% Number of GRF peaks
%     figure; hold on
%     
%     surf(k_total*1e3,ME_total,BOApks)
% 
% 
%     xlabel('Leg stiffness, k (kN/m)','FontSize',FS)
%     ylabel('Mechanical Energy (J)','FontSize',FS)
%     zlabel('Number of vGRF peaks','FontSize',FS)
% 
%     view(90,0)
%     set(gca,'FontSize',13)
%     box on

%% Old code
%     figure; hold on
%     for i = 1:16
%         surf(k_total*1e3,ME_total,AllAlpha(:,:,i),'FaceColor',[1 0 0],'FaceAlpha',0.10)
%     end
% 
%     FS = 17;
% 
%     % xlabel('Angle of attack, \alpha, deg')
%     xlabel('Leg stiffness, k (kN/m)','FontSize',FS)
%     ylabel('Mechanical Energy (J)','FontSize',FS)
%     zlabel('Alpha (rad)','FontSize',FS)
% 
%     view(90,0)
%     set(gca,'FontSize',13)
%     box on
    
% figure; hold on
% for i = 1:15
%     rdotDiff(:,:,i) = AllrdotStart(:,:,i+1)-AllrdotEnd(:,:,i);
%     surf(k_total*1e3,ME_total,rdotDiff(:,:,i),'FaceColor',[0 1 0],'FaceAlpha',0.10)
% end
% 
% 
% % xlabel('Angle of attack, \alpha, deg')
% xlabel('Leg stiffness, k (kN/m)','FontSize',FS)
% ylabel('Mechanical Energy (J)','FontSize',FS)
% zlabel('Change in spring velocity rdot','FontSize',FS)
% 
% view(90,0)
% set(gca,'FontSize',13)
% box on

%% If we need asymmetries
% if mod(i,2)
%     % red - physiological
%     surf(x_axisall,z_axisall,alltf(:,:,i),'FaceColor',[1 0 0],'FaceAlpha',0.25) 
%     hold on
% else
%     % green - prosthetic
%     surf(x_axisall,z_axisall,alltf(:,:,i),'FaceColor',[0 1 0],'FaceAlpha',0.25) 
% %% to see bifurcation
% end