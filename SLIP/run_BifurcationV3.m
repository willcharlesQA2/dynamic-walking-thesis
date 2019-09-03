
%% Basin of Attraction for spring one mass model
clear all

tic

step = 1;
% Initial Conditions
Param.tl = 0.4712; % toe length
Param.hl = -0.4712; % heel length (keep it -ve)
Param.frA = 0.3;   % foot radius
Param.frB = 0.3;   
Param.L0 = 1;
Param.m = 80;
% Minimum Theta, when leg becomes locked
% Param.Th = Param.hl/Param.fr;
Param.Th = -pi/2;
% Maximum Theta, when leg becomes locked
% Param.Tt = Param.tl/Param.fr;
Param.Tt = pi/2;

maxsteps = 400;

%% From data at start of step
% th_0 = -0.2179;
% rdot_0 = -0.2206;
% r_0 = -0.01858;
% thdot_0 = 1.068;

th_0 = 0;
rdot_0 = 0;

plotSize = 241;

ME_total = linspace(800,900,plotSize);
alpha_total = linspace(14,25,plotSize)*-pi/180;
k_total = linspace(15,30,plotSize)*1000;
fr_total = linspace(0,0.5,plotSize);
% k_total = linspace(20,26,21)*1000;

r_total = linspace(-0.0,-0.07,31);
rdot_total = linspace(-0.1,0.1,31);
thdot_total = linspace(1.2,1.5,21);

% r_total = r_0;
% thdot_total = thdot_0;

% Save all stable points
stb = 1;
% ME_0 = 832.2378;

% Number of steps plotted at the end
Bifs = 16;

% Prelocate to increase speed
maxEigenVals = NaN(size(k_total,2),1);
BOA = zeros(size(k_total,2),1);
GCr =   NaN(size(k_total,2),Bifs);
GCthdot = GCr;
GCrdot  = GCr;
ws = GCr;
ICs = GCr;
XdotIC = GCr;
YdotIC = GCr;
T1 = GCr;
T2 = GCr;
check = GCr;
interleg = GCr;
GRFPeak1 = NaN(size(k_total,2),Bifs-2);
GRFPeak2 = GRFPeak1;
StanceT = GRFPeak1;
XNetImpulse = GRFPeak1;

Param.alphaA    = 20*-pi/180;
Param.alphaB    = 20*-pi/180;
Param.kA     = 25*1e3;
Param.kB     = 25*1e3;
Param.ME_0  = 860;


% alpha_total = -0.38;
% x_axis = k_total;
x_axis = alpha_total*-180/pi;
% x_axis = fr_total;

eigenPlots = 0;

% y = 12.5/k*x -80*9.81/k;
for h = 1:size(k_total,2)
    for i = 1
        for j = 1
%             Param.frB = fr_total(h);
%             Param.alphaA = alpha_total(h);
            Param.alphaB = alpha_total(h);
%             Param.kA = k_total(h);
%             Param.kB = k_total(h);
%             Param.kB = k_total(h);
%             Param.ME_0 = ME_total(h);
%             r_0 = r_total(j);
    %         thdot_0 = thdot_total(j);
            th_0 = 0;
            r_0 = -0.01425;
            rdot_0 = 0;

            %% Find theta from r_0 and ME_0
            kS = Param.kA;          % Initial k value
            thdot_0 = sqrt( (2*Param.ME_0-kS*r_0^2 - 2*Param.m*9.81*(Param.L0+r_0))/(Param.m*(Param.L0+r_0)^2) );
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

                if mod(step,2) == 1 % odd step number
                    Param.k1 = Param.kA;
                    Param.k2 = Param.kB;
                    Param.alpha = Param.alphaB;
                    Param.fr1 = Param.frA;
                    Param.fr2 = Param.frB;
%                     legI = 'odd';
                else                % even step number
                    Param.k2 = Param.kA;
                    Param.k1 = Param.kB;
                    Param.alpha = Param.alphaA;
                    Param.fr2 = Param.frA;
                    Param.fr1 = Param.frB;
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
                            [ phi,~,~,~] = newConditions( u1(end,:),Param );
                            interleg(h,indice) = phi - u1(end,1);
                            
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
                            fprintf('x')
                            [ eigenVals(h,:) ] = localStability( qnew,Param);
                            maxEigenVals(h,:) = max(abs(eigenVals(h,:)));
                            % GRFs
    %                         ME(end)
                            BOA(h) = step;
                            %% Find GRF of each leg
                            for loopF = 3:Bifs % number of steps. Starts at 3 so that odd/even bifurcations are not messed up.
                                tend = -t_all{loopF-1}.t2(1); % last t value
                                TP1 = tend+t_all{loopF-1}.t2;
                                tend = tend+t_all{loopF-1}.t2(end);
                                TP1 = [TP1(1:end-1);tend+t_all{loopF}.t1(1:end-1);tend+t_all{loopF}.t2];
                                FP1 = [GRF_all{loopF-1}.FY2L(1:end-1);GRF_all{loopF}.FY1(1:end-1);GRF_all{loopF}.FY2T];
                                FPx = [GRF_all{loopF-1}.FX2L(1:end-1);GRF_all{loopF}.FX1(1:end-1);GRF_all{loopF}.FX2T];
%                                 plot(TP1,FP1)
                                [pks,locs] = findpeaks(FP1,TP1);
                                GRFPeak1(h,loopF-2) = pks(1);
                                GRFPeak2(h,loopF-2) = pks(2);
                                XNetImpulse(h,loopF-2) = trapz(TP1,FPx);
                            end
                            %% Calculate Stance time (could just use TP1(end)-TP1(1)?
                            for loopF = 3:Bifs
                                StanceT(h,loopF-2) = T2(h,loopF-1)+T1(h,loopF)+T2(h,loopF);
                            end
                        end
                    case 'noTouchDown'
                        %             fprintf('\nWalker has fallen over at step %g due to no touch down.\n',step)
                        BOA(h) = step;
                        fprintf('-')
                        break
                    case 'noTakeOff'
                        %             fprintf('\nWalker has fallen over at step %g due to no take off.\n',step)
                        BOA(h) = step;
                        fprintf('-')
                        break
                end

            end
        end
    end
%     fprintf('\n')
end
fprintf('\n')
% figure
% % surf(BOA)
% surf(r_total,alpha_total,BOA)
% xlabel('r')
% ylabel('alpha')
toc
Xlabel = 'Angle of attack of leg B, \alpha_B (deg)';
Sym = 20;
save('Bif_AlphaBV3.mat')
% save('Bif_RBV3.mat')

%% Set up eigenvalue figures
if eigenPlots == 1
    figure(1)
    hold on
    grid on
    viscircles([0,0],1);

    figure(2)
    hold on

    %% Plot eigenvalues
    figure(1)
    plot(eigenVals,'k.','MarkerSize',10)

    figure(2)
    plot(x_axis,maxEigenVals,'k.','MarkerSize',10)
end



