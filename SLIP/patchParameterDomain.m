% figure; hold on
% 
% for h = 1:size(alpha_total,2)
%     for i = 1:size(k_total,2)
%         for j = 1:size(ME_total,2)
%           if BOA(h,i,j) == maxsteps
%               plot3(alpha_total(h)*-180/pi,k_total(i)/1000,ME_total(j),'b*')
%           end
%       end
%   end
% end
% 
% xlabel('Slope angle, \alpha, deg')
% ylabel('Stiffness, k, kN/m')
% zlabel('Mechanical Energy, J')

EdgeColor = 'none';%[0 0 0.5];
FaceColor = 'flat';%'blue';
FaceAlpha = 1;
Colours = ['b';'g';'r';'c';'m';'y'];
%% Colour
fvcBlue2 = [52,106,157;
53,144,186;
79,152,226;
70,166,224;
129,167,218;
87,185,228]/255;

figure; hold on

% Represents half the width of each edge. Used later to plot each cube.
if size(alpha_total,2) == 1
    x1 = 0;
else
    x1 = (alpha_total(2)-alpha_total(1))/2 * -180/pi;
end
x2 = (k_total(2)-k_total(1))/2 /1000;
x3 = (ME_total(2)-ME_total(1))/2;

for h = 1:size(alpha_total,2)
    for i = 1:size(k_total,2)
        for j = 1:size(ME_total,2)
            if BOA(h,i,j) == maxsteps
                %               plot3(alpha_total(h)*-180/pi,k_total(i)/1000,ME_total(j),'b*')
                
                alpha = alpha_total(h)*-180/pi;
                k = k_total(i)/1000;
                ME = ME_total(j);
                
                % vertices of the cube
                vert = [alpha-x1, k-x2, ME-x3; alpha+x1, k-x2, ME-x3;...
                        alpha+x1, k+x2, ME-x3; alpha-x1, k+x2, ME-x3;...
                        alpha-x1, k-x2, ME+x3; alpha+x1, k-x2, ME+x3;...
                        alpha+x1, k+x2, ME+x3; alpha-x1, k+x2, ME+x3];
                fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
%                 if NumPeaks(h,i,j) == 0
%                     warning('No peaks detected')
%                     Col = 'k';
%                 else
%                     Col = Colours(NumPeaks(h,i,j));
%                 end
              
                Col = fvcBlue2;
                
                % Fill in a blue cube for each point.
                patch('Vertices',vert,'Faces',fac,'FaceVertexCData',Col,'FaceColor',FaceColor,'FaceLighting','gouraud','EdgeColor',EdgeColor,'linewidth',0.1,'FaceAlpha',FaceAlpha) %- not good for 3D
%                 patch('Vertices',vert,'Faces',fac,'FaceColor',Col,'FaceLighting','gouraud','EdgeColor',EdgeColor,'linewidth',0.1) %- not good for 3D

            end
        end
    end
end



xlabel('Angle of attack, \alpha, deg')
ylabel('Leg stiffness, k, kN/m')
zlabel('Mechanical Energy, J')

view(125, 20); grid on