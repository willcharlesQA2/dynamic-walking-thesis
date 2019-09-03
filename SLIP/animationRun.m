
%% Basin of Attraction for spring one mass model
clear all

tic

step = 1;
% Initial Conditions
Param.tl = 0.1178; %0.4712; % toe length
Param.hl = -0.1178; %-0.4712; % heel length (keep it -ve)
Param.fr1 = 0.3;   % foot radius
Param.fr2 = 0.3;   % foot radius
Param.L0 = 1;
Param.m = 80;
Param.kA = 25e3;
Param.kB = 25e3;
Param.alphaA    = 20*-pi/180;
Param.alphaB    = 20*-pi/180;
Param.ME_0  = 860;
% Minimum Theta, when leg becomes locked
Param.Th = Param.hl/Param.fr1;
% Maximum Theta, when leg becomes locked
Param.Tt = Param.tl/Param.fr1;

maxsteps = 200;

%% From data at start of step
% th_0 = -0.2179;
% rdot_0 = -0.2206;
% r_0 = -0.01858;
% thdot_0 = 1.068;

th_0 = 0;
rdot_0 = 0;

% r_total = r_0;
% thdot_total = thdot_0;

% Save all stable points
stb = 1;
% ME_0 = 832.2378;

% Number of steps plotted at the end
Bifs = 16;

% Prelocate to increase speed
BOA = zeros(1,1);
GCr =   NaN(1,Bifs);
GCthdot = GCr;
GCrdot  = GCr;
ws = GCr;
ICs = GCr;
XdotIC = GCr;
YdotIC = GCr;
T1 = GCr;
T2 = GCr;
check = GCr;
GRFPeak1 = NaN(1,Bifs-2);
GRFPeak2 = GRFPeak1;
StanceT = GRFPeak1;






th_0 = 0;
r_0 = -0.01425;
rdot_0 = 0;

%% Find theta from r_0 and ME_0
kS = Param.kA;          % Initial k value
thdot_0 = sqrt( (2*Param.ME_0-kS*r_0^2 - 2*Param.m*9.81*(Param.L0+r_0))/(Param.m*(Param.L0+r_0)^2) );

% Used for animation
tstart = 0;
h = 1; % used if we had bifurcations
total_M = [];
for step = 1:maxsteps
    
    if mod(step,2) == 1 % odd step number
        Param.k1 = Param.kA;
        Param.k2 = Param.kB;
        Param.alpha = Param.alphaB;
        %                     legI = 'odd';
    else                % even step number
        Param.k2 = Param.kA;
        Param.k1 = Param.kB;
        Param.alpha = Param.alphaA;
        %                     legI = 'even';
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
            
            % Check ME for every step
            %                     [ME1,KE1,PE1] = energy(t1,u1,Param,'SS');
            %                     [ME2,KE2,PE2] = energy(t2,u2,Param,'DS');
            %
            %                     fprintf('i = %g\tj = %g\nME1 = %g\nME2 = %g\n\n',i,j,ME1(1),ME2(end))
            if step > maxsteps - Bifs
                indice = step-maxsteps+Bifs;
                %%Animation
                [M1,tstart] = animation_frames(t1,u1,Param,'SS',step,tstart);
                [M2,tstart] = animation_frames(t2,u2,Param,'DS',step,tstart);
                total_M = [total_M,M1,M2];
                %% Other Gait Characteristics
                GCthdot(h,indice) = interp1(u1(:,1),u1(:,3),0);
                GCr(h,indice)     = interp1(u1(:,1),u1(:,2),0);
                GCrdot(h,indice)  = interp1(u1(:,1),u1(:,4),0);
                [ walking_speed,GRFStance ] = gaitCharacteristics( t1,t2,u1,u2,IC,Param );
                %                         plot(t1s,GRFs);
                ws(h,indice) = walking_speed;
                ICs(h,indice) = Param.IC;
                T1(h,indice) = t1(end)-t1(1);
                T2(h,indice) = t2(end)-t2(1);
                check(h,indice) = step;
                %                             [ XdotIC(h,indice),YdotIC(h,indice) ] = velocties( u1(end,1),u1(end,2),u1(end,3),u1(end,4),Param );
                
                [ FX,FY ] = GRFLOA2( t1,u1,t2,u2,Param );
                % F1 = Stance leg   SS
                % F2 = Trailing leg DS
                % F3 = Lead leg     DS
                
                GRF_all{indice,:} = struct('FX1',FX{1},'FX2T',FX{2}...
                    ,'FX2L',FX{3},'FY1',FY{1},'FY2T',FY{2},'FY2L',FY{3});
                
                t_all{indice,:}   = struct('t1',t1,'t2',t2);
                
                % to read do: C = permute(GC1,[1 3 2]);
            end
            if step == maxsteps
                
                BOA(h) = step;
                %% Find GRF of each leg
                for loopF = 3:Bifs % number of steps. Stars at 3 so that odd/even bifurcations are not messed up.
                    tend = -t_all{loopF-1}.t2(1); % last t value
                    TP1 = tend+t_all{loopF-1}.t2;
                    tend = tend+t_all{loopF-1}.t2(end);
                    TP1 = [TP1(1:end-1);tend+t_all{loopF}.t1(1:end-1);tend+t_all{loopF}.t2];
                    FP1 = [GRF_all{loopF-1}.FY2L(1:end-1);GRF_all{loopF}.FY1(1:end-1);GRF_all{loopF}.FY2T];
                    %                                 plot(TP1,FP1)
                end
                %% Calculate Stance time (could just use TP1(end)-TP1(1)?
                for loopF = 3:Bifs
                    StanceT(h,loopF-2) = T2(h,loopF-1)+T1(h,loopF)+T2(h,loopF);
                end
            end
        case 'noTouchDown'
            %             fprintf('\nWalker has fallen over at step %g due to no touch down.\n',step)
            fprintf('fallen\n')
            break
        case 'noTakeOff'
            %             fprintf('\nWalker has fallen over at step %g due to no take off.\n',step)
            fprintf('fallen\n')
            break
    end
    
end
toc




