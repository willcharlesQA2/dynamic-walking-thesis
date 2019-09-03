function [ eigenVals ] = linearised( x,step )
%Linearises the function and calculates a return map
% initial conditions
% x = th1_0,th2_0,th1dot_0,th2dot_0

%(gF)dx = F(x+dx)-x
% Need to find gF

% perturbations dx
p = rand(1,4)*1e-3;
tau = diag(p);

% x = th1_0,th2_0,th1dot_0,th2dot_0

for i = 1:4
    [t,u,th1sw,th2sw,th1dotsw,th2dotsw,event,steptime,Dist]=lagrangian_approach(x(1)+tau(i,1),x(2)+tau(i,2),x(3)+tau(i,3),x(4)+tau(i,4),step);

    psi(i,:) = [th1sw,th2sw,th1dotsw,th2dotsw] - x;
end

gF = tau\psi;

eigenVals = eig(gF);
end

