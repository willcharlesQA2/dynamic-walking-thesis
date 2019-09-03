function [ t_all,GRF_all ] = GRFfromPS( alpha,k,ME,Param )
% Runs the simulation again then plots the GRF plot

maxsteps = 50;
% Param.Aa = alpha;
Param.alphaA = alpha;
Param.alphaB = alpha;
Param.k = k;
Param.kA = k;
Param.kB = k;
Param.ME_0 = ME;

%Initial conditions (need to add randomisation)
Bifs = 3;
th_0 = 0;
r_0 = -Param.m*9.81/Param.k/2;
r_0 = -0.015;
% r_0 = -0.034;
rdot_0 = 0;

%% Find theta from r_0 and ME_0
thdot_0 = sqrt( (2*Param.ME_0-Param.k*r_0^2 - 2*Param.m*9.81*(Param.L0+r_0))/(Param.m*(Param.L0+r_0)^2) );

for step = 1:maxsteps
                if mod(step,2) == 1 % odd step number
                    Param.k1 = Param.kA;
                    Param.k2 = Param.kB;
                    Param.alpha = Param.alphaB;
                else                % even step number
                    Param.k2 = Param.kA;
                    Param.k1 = Param.kB;
                    Param.alpha = Param.alphaA;
                end
    [t1,u1,t2,u2,qnew,event,IC] = oneStep(th_0,r_0,thdot_0,rdot_0,Param);

    Param.IC = IC;

    switch event
        case 'none'
            % successful step
            th_0    = qnew(1);
            r_0     = qnew(2);
            thdot_0 = qnew(3);
            rdot_0  = qnew(4);

            if step > maxsteps - Bifs
%                [ ~,GRFStance ] = gaitCharacteristics( t1,t2,u1,u2,IC,Param );
            [ FX,FY ] = GRFLOA2( t1,u1,t2,u2,Param );
            % F1 = Stance leg   SS
            % F2 = Trailing leg DS
            % F3 = Lead leg     DS
            
            GRF_all{step-(maxsteps-Bifs),:} = struct('FX1',FX{1},'FX2T',FX{2}...
                ,'FX2L',FX{3},'FY1',FY{1},'FY2T',FY{2},'FY2L',FY{3});
%             GRF_all{step-(maxsteps-Bifs),:} = struct('FY',FY);
            t_all{step-(maxsteps-Bifs),:}   = struct('t1',t1,'t2',t2);
%             t_all{step-(maxsteps-Bifs),:}   = struct('t2',t2);
%                 GRF_all{step-(maxsteps-Bifs)} = GRFStance;
                
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

