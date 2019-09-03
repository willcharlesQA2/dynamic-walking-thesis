function [ walking_speed,GRFStance ] = gaitCharacteristics( t1,t2,u1,u2,IC,Param )
% Use to collect Gait characteristics such as walking speeds, peak GRFs, 
% step length, step period, inter-leg angle. (and thdot)

% IC - initial contact, X position of other leg at \theta_2 = 0

%%GRFs
GRFStance = [-Param.k1*u1(:,2).*cos(u1(:,1));-Param.k1*u2(:,2).*cos(u2(:,1))];

% FZ1 =
% hint: take the first half of gait and then find MAX. This ignores
% double-support phase at the start.

walking_speed = Param.IC/t2(end);


end

