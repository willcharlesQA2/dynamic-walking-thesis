%% readAllResults
% Take data gained from previous studies and read results.
% Extra things to include:
% Inital slopes and stiffnesses
% Param

% Damping values used
Dmps = 50:50:300;

col1 = [0, 0, 1];   % blue
col2 = [1, 0, 0];   % red
length = size(Dmps,2);
colours = [linspace(col1(1),col2(1),length)', linspace(col1(2),col2(2),length)', linspace(col1(3),col2(3),length)'];
% Dmp is damping ratio
c = 0;
figure
for Dmp = Dmps
    c = c+1;
    colour = colours(c,:);
    load(['Compiled Results/resultsDamped',num2str(Dmp),'.mat'])
%     figure
    for  i = 1:size(StableSolution,1)
%         Value = StableSolution{i,2}(8:end);     % Walking Speeds
        Value = AllGRFs{i,7}(end,2)*-180/pi;
        h(c) = plot3(StableSolution{i,4},StableSolution{i,3}*180/pi,Value,'o','MarkerFaceColor',colour,'MarkerEdgeColor','none');
        hold on
    end
    
end

grid on
xlabel('Leg Stiffness (kN/m)')
ylabel('Slope angle (deg)')
% zlabel('Walking Speed (m/s)')
zlabel('Inter-leg angle (deg)')
% legend(h, '50','100','150','200','250','300'); % or num2str(Dmps)
legend(h, num2str(Dmps'));

%% Example of saved results
% save('resultsDampedSA.mat','StableSolution','Speeds','allSteps','AllGRFs')
% add 'Param' for each step as well as slopes and stiffnesses for all
% solutions. Param can be under StableSolution?
%% GRF saved
% GRFs = {GRFy1,GRFy2,GRFx1,GRFx2,t1,t2,u1,u3,Param};
%% Find GRFs
% Take individual results and analyse GRFs. Maybe even to re-run results to
% gain t1,v1 etc.
GRFplot(AllGRFs,1)

%% Annotate figure?
% When a GRFplot is taken, an arrow will be plotted to the workspace so
% that we will know where the plot is taken from? Possibility of having the
% GRF plot on the same figure or re-plotting the data?

% Is there a way to automatically determine if a GRF plot has more/less
% than two humps?