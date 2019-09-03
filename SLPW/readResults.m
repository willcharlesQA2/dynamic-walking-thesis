%% readResults 
% Take data gained from previous studies and read results.
% Extra things to include:
% Inital slopes and stiffnesses 
% Param

% load('resultsDamped100v2.mat')
load('resultsDamped150NEW2.mat')
figure
for  i = 1:size(StableSolution,1)
% plot3(StableSolution{i,4},StableSolution{i,3}*180/pi,Speeds(i,8:end),'bo')
plot3(StableSolution{i,4},StableSolution{i,3}*180/pi,StableSolution{i,2}(8:end),'bo')
hold on
end
grid on
xlabel('Leg Stiffness (kN/m)')
ylabel('Slope angle (deg)')
zlabel('Walking Speed (m/s)')

%% Example of saved results
% save('resultsDampedSA.mat','StableSolution','Speeds','allSteps','AllGRFs')
% add 'Param' for each step as well as slopes and stiffnesses for all
% solutions. Param can be under StableSolution?

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