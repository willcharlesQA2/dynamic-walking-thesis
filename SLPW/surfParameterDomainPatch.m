% load('resultsParameterExperimentSwoma3.mat')
% load('resultsParameterAlpha.mat')
colourPlot = 1;

k_total = stiffnesses;

ME_total = MEs;

AA = sizeA;
BB = sizeB;

for i = 1:15
    rdotDiff(:,:,i) = AllrdotStart(:,:,i+1)-AllrdotEnd(:,:,i);
%     [totalMomDiff(:,:,i),~] = momentumDiffCollision( Allu1End{:}(i,:),Allu2Start{:}(i,:),Param );
    %This is the old one
%      [totalMomDiff(imp,bifurcation),~] = momentumDiffCollision( Allu1End{:,bifurcation}(imp,:),Allu2Start{:,bifurcation}(imp+1,:),Param );
%     surf(k_total*1e3,ME_total,rdotDiff(:,:,i),'FaceColor',[0 1 0],'FaceAlpha',0.10)
end

totalMomDiff = nan(AA,BB);
for aa = 1:AA
    for bb = 1:BB
        if isnan(BOA(aa,bb)) == 0
            for i = 1:15
                [totalMomDiff(aa,bb,i),~] = momentumDiffCollision( Allu1End{aa,bb}(i,:),Allu2Start{aa,bb}(i+1,:),Param );
            end
        end
    end
end
Colours = {[1 0 0],[0 1 0],[0 0 1]};
zlabels = {'Alpha (deg)','Change in spring velocity rdot','Walking Speed (m/s)','Number of vGRF peaks','Collision impulse (Ns)'};
zParam = {AllAlpha*-180/pi,rdotDiff,AllSpeeds,BOApks,totalMomDiff};
for Loop = 1:size(zParam,2)
    zParamDiff{Loop} = diff(zParam{Loop},1,3);
end
zParamDiff{4} = zeros(size(zParam{Loop}));

% Bifurcation Size
iSize = {16,15,16,1,15};

EdgeColor = 'none';%[0 0 0.5];
FaceColor = 'flat';%'blue';
Col = 'b';

x1 = (k_total(2)-k_total(1))/2 /1000;
x2 = (ME_total(2)-ME_total(1))/2;



for plotLoop = 1:5
    Max(plotLoop) = max(zParam{plotLoop}(:));
    Min(plotLoop) = min(zParam{plotLoop}(:));
    MaxDiff(plotLoop) = max(zParamDiff{plotLoop}(:));
    MinDiff(plotLoop) = min(zParamDiff{plotLoop}(:));
    figure; hold on
    for aa = 1:AA
        for bb = 1:BB
            if BOA(aa,bb) == 1
                ME = ME_total(aa);
                k = k_total(bb)/1000;

                vert = [k-x1, ME-x2; k+x1, ME-x2; ...
                        k+x1, ME+x2; k-x1, ME+x2     ];
                fac = [1 2 ...
                       3 4 ];
                   
                   % zParam{plotLoop}(:,:,i) for colour
                  %  zParamDiff for black bifurcation?
                  % C = (x - Min) / (Max - Min)
%                   C = (zParam{plotLoop}(aa,bb,1) - Min(plotLoop)) / (Max(plotLoop) - Min(plotLoop));
%                 percDiff = abs(zParamDiff{plotLoop}(aa,bb,:))./abs(zParam{plotLoop}(aa,bb,1:end-1));

                if abs(zParamDiff{plotLoop}(aa,bb,:)) < 1e-2
                    C = zParam{plotLoop}(aa,bb,1);
                    FaceColor = 'flat';
                    patch('Vertices',vert,'Faces',fac,'FaceVertexCData',C,'FaceColor',FaceColor,'FaceLighting','gouraud','EdgeColor',EdgeColor,'linewidth',0.1)
                else
                    FaceColor = 'k';
                    patch('Vertices',vert,'Faces',fac,'FaceColor',FaceColor,'FaceLighting','gouraud','EdgeColor',EdgeColor,'linewidth',0.1)
%                     C = [1,1,1];
%                     C = [0,0,0];
                end
%                 patch('Vertices',vert,'Faces',fac,'FaceVertexCData',C,'FaceColor',FaceColor,'FaceLighting','gouraud','EdgeColor',EdgeColor,'linewidth',0.1)
        % Bifurcation 1-16

            end
        end
    end
    FS = 17;

    h = colorbar;
    colormap(jet)
    
    % xlabel('Angle of attack, \alpha, deg')
    xlabel('Leg stiffness, k (kN/m)','FontSize',FS)
    ylabel('Mechanical Energy (J)','FontSize',FS)
%     zlabel(zlabels{plotLoop},'FontSize',FS)
    ylabel(h,zlabels{plotLoop},'FontSize',FS)

    
    view(0,90)
    set(gca,'FontSize',13)
    box on
end

%% Surf
% for i = 1:iSize{plotLoop}
%     if colourPlot == 1
%         surf(k_total*1e3,ME_total,zParam{plotLoop}(:,:,i)) %,'FaceAlpha',1/iSize{plotLoop}
%         colormap(jet)
%     else
%         surf(k_total*1e3,ME_total,zParam{plotLoop}(:,:,i),'FaceColor',Colours{plotLoop},'FaceAlpha',0.10)
%     end
% end
                
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