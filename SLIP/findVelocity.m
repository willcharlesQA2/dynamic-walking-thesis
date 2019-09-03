function findVelocity
%% Find the set of initial conditions (th1dot and r) to find the required 
% walking velocity
clear all
step = 1;
% Initial Conditions
Param.tl = 0.4712; % toe length
Param.hl = -0.4712; % heel length (keep it -ve)
Param.fr = 0.3;   % foot radius
Param.L0 = 1;
Param.m = 80;
Param.k = 16e3;
Param.Aa = -20*(pi/180);
% Minimum Theta, when leg becomes locked
Param.Th = Param.hl/Param.fr;
% Maximum Theta, when leg becomes locked
Param.Tt = Param.tl/Param.fr;

%% From data at start of step
% th_0 = -0.2179;
% rdot_0 = -0.2206;
% r_0 = -0.01858;
% thdot_0 = 1.068;

th_0 = 0;
rdot_0 = 0;
Param.ME = 852.2378;

% a = 1/2*Param.m*thdot_0^2+1/2*Param.k;
% b = thdot_0^2*Param.L0*Param.m + Param.m*9.81;
% c = Param.m*9.81*Param.L0 + 1/2*Param.m*Param.L0^2*thdot_0^2 - MEsimple;
% 
% x = roots([a b c]);    
% 
% r_0 = min(x);

% r_total = linspace(-0.0,-0.06,31);
% thdot_total = linspace(0.2,1.5,31);
x0 = [-0.036,-20*(pi/180)];
% Don't work!!!!!!
% optimset('TolFun',1e-6,'TolX',1e-6,'display','iter');
% options = optimset('TolFun',1e-12);
[x,fval] = fsolve(@(x)limitcycle(x,Param,th_0,rdot_0),x0)
% r_total = r_0;
% thdot_total = thdot_0;
% figure
% hold on
% xlabel('r')
% ylabel('thdot')
% for i = 1:size(r_total,2)
%     for j = 1:size(thdot_total,2)
%         r_0 = r_total(i);
%         thdot_0 = thdot_total(j);
%         maxsteps = 50;
%         th_0 = 0;
%         rdot_0 = 0;
% 
%         % Need to change oneStep to end at the next apex condition
%         [th1dot,t1,u1,t2,u2,t3,u3,qnew,event,IC] = Apex2Apex(th_0,r_0,thdot_0,rdot_0,Param);
%         %     [t1,u1,t2,u2,qnew,event,IC] = oneStep(th_0,r_0,thdot_0,rdot_0,Param);
%         if strcmp(event,'none')
%             fprintf('Finally at r = %g and thdot = %g\n',r_0,thdot_0)
%             plot(r_0,thdot_0,'b.')
% %              Param.IC = IC;
%         end
%        
%     end
% end
end

function [Fout] = limitcycle(x0,Param,th_0,rdot_0)
ME_0 = Param.ME;
r_0 = x0(1);
Param.Aa = x0(2);

        thdot_0 = sqrt( (2*ME_0-Param.k*r_0^2 - 2*Param.m*9.81*(Param.L0+r_0))/(Param.m*(Param.L0+r_0)^2) );


        [t1,u1,t2,u2,t3,u3,qnew,event,IC] = Apex2Apex(th_0,r_0,thdot_0,rdot_0,Param);
        if isnan(qnew)
            Fout1 = NaN;
            Fout2 = NaN;
        else 
            Fout2 = qnew(3) - thdot_0;
            Fout2 = qnew(4);
            Fout1 = qnew(2) - r_0;
        end
   Fout = [Fout1 Fout2];

end