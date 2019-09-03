% load('resultsParameterExperimentSwoma3.mat')
% If there is 
oneCol = 0;

k_total = stiffnesses;

ME_total = MEs;

EdgeColor = 'none';%[0 0 0.5];
FaceColor = 'flat';%'blue';
Colours = ['g';'b';'r';'c';'m';'y';'y';'y'];
%% Colour
fvcBlue2 = [52,106,157;
53,144,186;
79,152,226;
70,166,224;
129,167,218;
87,185,228]/255;

figure; hold on

% Represents half the width of each edge. Used later to plot each cube. 

x1 = 0;

x2 = (k_total(2)-k_total(1))/2 /1000;
x3 = (ME_total(2)-ME_total(1))/2;
% plot3(StableSolution{i,4},StableSolution{i,3},StableSolution{i,2}(8:end),'bo')
% StableSolution(bifurcation,:) = {uStable,speed,MEs(aa),stiffnesses(bb),Param};
% for h = 1
figure; hold on;
    for i = 1:size(ME_total,2)
        for j = 1:size(k_total,2)
            if BOA(i,j) == 1
                %               plot3(alpha_total(h)*-180/pi,k_total(i)/1000,ME_total(j),'b*')
                
                alpha = 0;
                ME = ME_total(i);
                k = k_total(j)/1000;
                
                % vertices of the cube
                vert = [alpha-x1, k-x2, ME-x3; alpha+x1, k-x2, ME-x3;...
                        alpha+x1, k+x2, ME-x3; alpha-x1, k+x2, ME-x3;...
                        alpha-x1, k-x2, ME+x3; alpha+x1, k-x2, ME+x3;...
                        alpha+x1, k+x2, ME+x3; alpha-x1, k+x2, ME+x3];
                fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
                if BOApks(i,j) == 0
                    warning('No peaks detected')
                    Col = 'k';
                elseif oneCol ==1
                    Col = 'b';
                else
                    Col = Colours(BOApks(i,j));
                end
%               Col = 'b';
%               Col = fvcBlue2
                
                % Fill in a blue cube for each point.
%                 patch('Vertices',vert,'Faces',fac,'FaceVertexCData',Col,'FaceColor',FaceColor,'FaceLighting','gouraud','EdgeColor',EdgeColor,'linewidth',0.1) %- not good for 3D
                patch('Vertices',vert,'Faces',fac,'FaceColor',Col,'FaceLighting','gouraud','EdgeColor',EdgeColor,'linewidth',0.1) %- not good for 3D

            end
        end
    end


% Font size
FS = 17;

% xlabel('Angle of attack, \alpha, deg')
ylabel('Leg stiffness, k (kN/m)','FontSize',FS)
zlabel('Mechanical Energy (J)','FontSize',FS)

view(90,0)
set(gca,'FontSize',13)
box on