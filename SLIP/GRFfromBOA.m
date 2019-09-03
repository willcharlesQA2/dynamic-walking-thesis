function [ t_all,GRF_all ] = GRFfromBOA( alpha,k,ME,Param )
% Runs the simulation again then plots the GRF plot

maxsteps = 20;
Param.Aa = alpha;
Param.k = k;
Param.ME_0 = ME;

%Initial conditions (need to add randomisation)
th_0 = 0;
r_0 = -Param.m*9.81/Param.k;
rdot_0 = 0;

%% Find theta from r_0 and ME_0
thdot_0 = sqrt( (2*Param.ME_0-Param.k*r_0^2 - 2*Param.m*9.81*(Param.L0+r_0))/(Param.m*(Param.L0+r_0)^2) );

for step = 1:maxsteps
    [t1,u1,t2,u2,qnew,event,IC] = oneStep(th_0,r_0,thdot_0,rdot_0,Param);

    Param.IC = IC;

    switch event
        case 'none'
            % successful step
            th_0    = qnew(1);
            r_0     = qnew(2);
            thdot_0 = qnew(3);
            rdot_0  = qnew(4);

            if step > maxsteps - 2
               [ ~,GRFStance ] = gaitCharacteristics( t1,t2,u1,u2,IC,Param );

%                 ws(h,i,j,step) = walking_speed;
                % to read do: C = permute(GC1,[1 3 2]);
                t_all{step-(maxsteps-2)}   = [t1;t2];
                GRF_all{step-(maxsteps-2)} = GRFStance;
                
            end

        case 'noTouchDown'
            %             fprintf('\nWalker has fallen over at step %g due to no touch down.\n',step)

            break
        case 'noTakeOff'
            %             fprintf('\nWalker has fallen over at step %g due to no take off.\n',step)

            break
    end
end

end

