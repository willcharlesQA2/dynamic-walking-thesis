function [ eigenVals ] = localStability( x,Param)
%Linearises the function and calculates a return map
% initial conditions
% x = th1_0,th2_0,th1dot_0,th2dot_0

%(gF)dx = F(x+dx)-x
% Need to find gF

% perturbations dx
p = rand(1,3)*1e-6;
tau = diag(p);

% x = th1_0,th2_0,th1dot_0,th2dot_0

% number of steps after perturbation
pSteps = 1;
% For each perturbation
for i = 1:3
%                th_0,         r_0,        rdot_0,        thdot_0
    qStart = [x(1)+tau(i,1),x(2)+tau(i,2),x(4)+tau(i,3)]; %,x(4)+tau(i,4)
    [ thdot_0 ] = InitialConditionsV2(qStart(1), qStart(2), qStart(3),Param );
    for step = 1:pSteps
        [t1,u1,t2,u2,qNew,event,IC] = oneStep(qStart(1),qStart(2),thdot_0,qStart(3),Param);
        qStart = qNew;
        if isnan(qNew) == 1
            qNew = nan(1,4);
        end
    end
    psi(i,:) = [qNew(1:2),qNew(4)] - [x(1:2),x(4)];
end

gF = tau\psi;

eigenVals = eig(gF);
end

