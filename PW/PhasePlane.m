function [  ] = PhasePlane(filename,indice,asymmetric)
% PhasePlane(filename,indice,asymmetric)
%Takes values from uAllBif and plots phase plane diagrams
% filname in '.mat' format
% indice, value of x_axis as 1, 2 or 3 etc.
% asymmetric as 1 if there are 2 different legs

load(filename)

if nargin == 1
    indice = 2;
    asymmetric = 0;     % if there are 2 different legs, set this as 1
end

% AllVar(bifurcation,:) = {alpha*-180/pi, uAll, DistAll, tf, inter, angvel};

x_axis = cell2mat(uAllbif(:,1));
uAll   = uAllbif(:,2);

alltf = cell2mat(AllVar(:,4)); 


%% Phase plane diagram
% u for the last step of the specified angle
maxsteps = size(alltf,2);
bif = 0;
% 2.5, 3, 3.5

if asymmetric
    bifurcations = 32;
else
    bifurcations = 1;
end
for i = maxsteps+1-bifurcations:maxsteps
    bif = bif + 1;
    uSlope = cell2mat(uAll{indice}(i,2));

    th1 = uSlope(:,1);
    phi = uSlope(:,2) + th1;
    th1dot = uSlope(:,3);
    phidot = uSlope(:,4) + th1dot;
    
    if mod(i,2) %odd
        if exist('leg1')
            leg1 = [leg1;th1];
            leg1dot = [leg1dot;th1dot];
            leg2 = [leg2;phi];
            leg2dot = [leg2dot;phidot];
        else
            leg1 = [th1];
            leg1dot = [th1dot];
            leg2 = [phi];
            leg2dot = [phidot];            
        end 
    else        %even
        if exist('leg1')
            leg1 = [leg1;phi];
            leg1dot = [leg1dot;phidot];  
            leg2 = [leg2;th1];
            leg2dot = [leg2dot;th1dot]; 
        else
            leg1 = [phi];
            leg1dot = [phidot];  
            leg2 = [th1];
            leg2dot = [th1dot]; 
        end
    end
end

%% Plot phase plane to show limit cycles
figure(1)
hold on
if asymmetric
    plot(leg1,leg1dot,'b')
    plot(leg2,leg2dot,'r')
    title(['Phase-plane with modified leg ROC at ',num2str(x_axis(indice)),'L'])
else
    plot([leg1;leg2;leg1(1)],[leg1dot;leg2dot;leg1dot(1)],'k--')
%     plot(leg1,leg1dot,'k--')
%     plot(leg2,leg2dot,'k--')
%     plot([leg1(end),leg2(1)],[leg1dot(end),leg2dot(1)],'k--')
%     plot([leg1(1),leg2(end)],[leg1dot(1),leg2dot(end)],'k--')
end
xlabel('Leg angle (rad)','FontSize',15)
ylabel('Leg angular velocity (rad/s)','FontSize',15)
set(gca,'FontSize',12)
box on

%% Labels and legends
hl(1) = plot(NaN,NaN,'b');
hl(2) = plot(NaN,NaN,'r');
hl(3) = plot(NaN,NaN,'k--');

if asymmetric
    legend(hl(1,:),'unaffected leg','modified leg','symmetric legs')
end
if mod(size(alltf,2),2)
    warning('There are an odd number of ''maxsteps''Check that even steps are prosthetic.')
end



end

