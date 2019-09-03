
%% Basin of Attraction for spring one mass model
clear all

tic

step = 1;
% Initial Conditions
Param.tl = 0.4712; % toe length
Param.hl = -0.4712; % heel length (keep it -ve)
Param.fr1 = 0.3;   % foot radius
Param.fr2 = Param.fr1;
Param.L0 = 1;
Param.m = 80;
Param.k = 20e3;
Param.Aa = -20*(pi/180);
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

ME_total = linspace(820,880,51);
alpha_total = linspace(10,40,51)*-pi/180;
k_total = linspace(5,30,51)*1000;

r_total = linspace(-0.0,-0.07,31);
rdot_total = linspace(-0.1,0.1,31);
thdot_total = linspace(1.2,1.5,21);
% r_total = r_0;
% thdot_total = thdot_0;

% Save all stable points
stb = 1;
% ME_0 = 832.2378;

Bifs = 4;

BOA = zeros(size(alpha_total,2),size(k_total,2),size(ME_total,2));
GCr =   NaN(size(alpha_total,2),size(k_total,2),size(ME_total,2),Bifs);
GCthdot = GCr;
GCrdot  = GCr;
ws = GCr;
ICs = GCr;
GRFPeak = BOA;
NumPeaks = BOA;
GRFPeak1 = NaN(size(alpha_total,2),size(k_total,2),size(ME_total,2),Bifs-2);
GRFPeak2 = GRFPeak1;
for h = 1:size(alpha_total,2)
    for i = 1:size(k_total,2)
        for j = 1:size(ME_total,2)
            Param.alpha = alpha_total(h);
            Param.k1 = k_total(i);
            Param.k2 = Param.k1;
            Param.ME_0 = ME_total(j);
%             r_0 = r_total(j);
    %         thdot_0 = thdot_total(j);
            th_0 = 0;
%             r_0 = -0.015;
            r_0 = -0.015;
            rdot_0 = 0;

            %% Find theta from r_0 and ME_0
            thdot_0 = sqrt( (2*Param.ME_0-Param.k*r_0^2 - 2*Param.m*9.81*(Param.L0+r_0))/(Param.m*(Param.L0+r_0)^2) );
            % To calculate the energy of the system assuming rdot = 0 at the
            % apex
    %         KEsimple = 1/2*Param.m*((Param.L0+r_0)*thdot_0)^2;
    %         PEsimple = 1/2*Param.k*r_0^2 + Param.m*9.81*(Param.L0+r_0);
    %         MEsimple = KEsimple + PEsimple;
    %         a = 1/2*Param.m*thdot_0^2+1/2*Param.k;
    %         b = thdot_0^2*Param.L0*Param.m + Param.m*9.81;
    %         c = Param.m*9.81*Param.L0 + 1/2*Param.m*Param.L0^2*thdot_0^2 - ME_0;
    %         x = roots([a b c]);
    %         r_0 = min(x);

    % Check ME initial
    %         fprintf('x = %g\tcheck = %g\n',thdot_0,MEsimple-ME_0)
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

                        % Check ME for every step
    %                     [ME1,KE1,PE1] = energy(t1,u1,Param,'SS');
    %                     [ME2,KE2,PE2] = energy(t2,u2,Param,'DS');
    % 
    %                     fprintf('i = %g\tj = %g\nME1 = %g\nME2 = %g\n\n',i,j,ME1(1),ME2(end))
                        if step > maxsteps - Bifs
                            indice = step-maxsteps+Bifs;
                            GCthdot(h,i,j,indice) = interp1(u1(:,1),u1(:,3),0);
                            GCr(h,i,j,indice)     = interp1(u1(:,1),u1(:,2),0);
                            GCrdot(h,i,j,indice)  = interp1(u1(:,1),u1(:,4),0);
                            [ walking_speed,GRFStance ] = gaitCharacteristics( t1,t2,u1,u2,IC,Param );
        %                         plot(t1s,GRFs);
                            ws(h,i,j,indice) = walking_speed;
                            ICs(h,i,j,indice) = Param.IC;

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
                            fprintf('x')
                            % GRFs
    %                         ME(end)
                            BOA(h,i,j) = step;
%                             GRFPeak(h,i,j) = max(GRFStance);
                            for loopF = 3:Bifs % number of steps. Stars at 3 so that odd/even bifurcations are not messed up.
                                tend = -t_all{loopF-1}.t2(1); % last t value
                                TP1 = tend+t_all{loopF-1}.t2;
                                tend = tend+t_all{loopF-1}.t2(end);
                                TP1 = [TP1(1:end-1);tend+t_all{loopF}.t1(1:end-1);tend+t_all{loopF}.t2];
                                FP1 = [GRF_all{loopF-1}.FY2L(1:end-1);GRF_all{loopF}.FY1(1:end-1);GRF_all{loopF}.FY2T];
%                                 plot(TP1,FP1)
                                [pks,locs] = findpeaks(FP1,TP1);
                                GRFPeak1(h,i,j,loopF-2) = pks(1);
%                                 GRFPeak2(h,i,j,loopF-2) = pks(2);
                            end
                            NumPeaks(h,i,j) = size(pks,1);
                        end
                    case 'noTouchDown'
                        %             fprintf('\nWalker has fallen over at step %g due to no touch down.\n',step)
                        BOA(h,i,j) = step;
                        fprintf('-')
                        break
                    case 'noTakeOff'
                        %             fprintf('\nWalker has fallen over at step %g due to no take off.\n',step)
                        BOA(h,i,j) = step;
                        fprintf('-')
                        break
                end

            end
        end
    end
    fprintf('\n')
end
% figure
% % surf(BOA)
% surf(r_total,alpha_total,BOA)
% xlabel('r')
% ylabel('alpha')
toc
save('PSZoomedIn.mat')



