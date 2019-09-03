function [fitness,Speed,GRFs,Param,results] = Walk(u,Param)
% Input
if nargin ~= 0
    q1_0 = u(1);
    q3_0 = u(2);
    %     q1dot_0 = u(3);
    q3dot_0 = u(3);
%     stiffness = u(4);%u(5);
%     slope = u(5);%u(6);
%     damping = u(6);
%     ME_0 = u(7);
end
% bifurcations == 1 for plotting
bifurcations = 1;
animation = 0;
writeMovie = 0;
energyPlots = 0;
phasePlane = 0;

% Maximum number of steps taken until it is decided the walker is stable
% (~50)
maxsteps = Param.maxsteps;

if isfield(Param,'Bifs')
    Bifs = Param.Bifs;
else
    Bifs = 16;
end

% Model Parameters
ME_0 = Param.ME_0;
slope = Param.slope;

% Param.xd2 = 0;
% Param.yd2 = Param.Lr2-Param.LH1;
% Param.xd1 = 0;
% Param.yd1 = Param.Lr1-Param.LH1;

Param.push = 0e3;
% Param.load = -0.012;
Param.load = 0;

totalM = Param.mA + Param.mB + Param.m3;


% Param.kB = stiffness;
% Param.c1 = damping;
% Param.c2 = damping;

%slope angle
Param.alpha = slope;
% Param.alpha = -0.5*pi/180;


%Bifurcation
% For breaking outside of for loops
bb = 0;


%% While the walker is in double-support, finds th2, r2 for a given th1, r1.
% [q2_0,q4_0,q2dot_0,q4dot_0] = constrained_DS(q1_0,q3_0,q1dot_0,q3dot_0);
[q2_0,q1dot_0,q2dot_0,q4dot_0 ] = constrained_ME( q1_0,q3_0,q3dot_0,ME_0,Param );
q4_0 = 0;

q = [q1_0;q2_0;q3_0;q4_0];
qdot = [q1dot_0;q2dot_0;q3dot_0;q4dot_0];

% Set time = 0 at first step.
t2temp = 0;

for step = 1:maxsteps
    indice = step-maxsteps+Bifs;
    
    %Numerical Integration
    % options=odeset('Events',@(t,u)heel_strike2(t,u),'abstol',10^(-8),'reltol',10^(-8));
    [ t1,u1,t2,u2,event ] = doubleFirst( q,qdot,t2temp,step,Param );
    
    switch event
        case 'none'
            

            
            %             fprintf('Step %g complete\n',step)

            %% ANIMATION
            if step >= maxsteps - 11 && animation == 1;
                
                % Gets frames for movie file
                leftover_time = 0;
                
                [M2,leftover_time] = animation_frames(t2,u2,step+1,leftover_time);
                [M1,leftover_time] = animation_frames(t1,u1,step,leftover_time);
                if exist('total_M') == 1
                    total_M = [total_M,M2,M1];
                else
                    total_M = [M2,M1];
                end
                
                %movie(total_M,1,20)
                % movie(M2,1,20)
                %plots(t,u,step)
                
            end
            
            %% PHASE-PLANE DIAGRAM
            if phasePlane == 1
                if mod(step,2) == 1
                    leg = [u2(:,1),u2(:,5),u2(:,3),u2(:,7);...
                        u1(:,1)+u1(:,2),u1(:,4)+u1(:,5),zeros(size(t1)),zeros(size(t1))];
                else
                    leg = [u2(:,1)+u2(:,2),u2(:,5)+u2(:,6),u2(:,4),u2(:,8);...
                        u1(:,1),u1(:,4),u1(:,3),u1(:,6)];
                end
                
                if exist('PP') == 1
                    PP = [PP;leg];
                else
                    PP = leg;
                end
            end
            
            %% Collision phase outside function?
            [ Qpost ] = momentumSpringREDUX( u1(end,1),u1(end,2),u1(end,3),0,u1(end,4),u1(end,5),u1(end,6),0,step,Param);
            qdot = [Qpost'];
            q = [u1(end,1),u1(end,2),u1(end,3),0];
            
            %             fprintf('\nInitial conditions for step %g heel-strike are:\n',step+1)
            %             fprintf('q1    = %g\nq2    = %g\nq3    = %g\nq4    = %g\n',q(1),q(2),q(3),q(4))
            %             fprintf('q1dot = %g\nq2dot = %g\nq3dot = %g\nq4dot = %g\n',qdot(1),qdot(2),qdot(3),qdot(4))
            
            % end time used to plot more than one step
            t2temp = t1(end);
            
            %% Ground Reaction Forces
            if step > maxsteps - Bifs
                indice = step-maxsteps+Bifs;
                [ GRFs{indice} ] = GRFLOA2( t1,u1,t2,u2,Param,step );
                
            end
            
            if step > maxsteps - Bifs && bifurcations == 1;
                indice = step-maxsteps+Bifs;
                %% Ground Reaction Forces
                [ GRFs{indice} ] = GRFLOA2( t1,u1,t2,u2,Param,step );
                
%                 results.stancetime(indice) = t1(end)-t1(1)+t2(end)-t2(1); % Or t1(end)-t2(1)
                results.T1(indice,:) = t1(end)-t1(1);
                results.T2(indice,:) = t2(end)-t2(1);
                results.th1End(indice,:)   = u1(end,1);
                results.th1dotStart(indice,:) = u2(1,5);
                results.th1dot(indice,:) = u1(end,4);
                results.th2dot(indice,:) = u1(end,5);
                results.rdotStart(indice,:) = u2(1,7);
                results.rEnd(indice,:)   = u1(end,3);
                results.rdotEnd(indice,:)   = u1(end,6);
                results.interleg(indice,:) = u1(end,2);
                results.u1End(indice,:) = u1(end,:);
                results.u2Start(indice,:) = u2(1,:);
                results.GCth1dot(indice,:)  = interp1(u1(:,1),u1(:,4),0);
                results.GCth2(indice,:)     = interp1(u1(:,1),u1(:,2),0);
                results.GCth2dot(indice,:)  = interp1(u1(:,1),u1(:,5),0);
                results.GCr(indice,:)       = interp1(u1(:,1),u1(:,3),0);
                results.GCrdot(indice,:)    = interp1(u1(:,1),u1(:,6),0);
%                 results.rdotIncrease(indice,:) = Qpost(3) - u1(end,6);
                results.touchDownAngle(indice,:) = u1(end,1)+u1(end,2);
                
                %% Walking speed
%                 [~,~,~,~,X3start,~,~,~,Xhs] = locations(t2(1,:),u2(1,:),Param);

                [~,~,~,~,~,~,~,~,Xhs] = locations(t1(end,:),u1(end,:),Param,step,'SS');
                Speed(indice) = (Xhs)/(t1(end)-t2(1));
                
                %                                 fprintf('Walking speed is %gm/s\n',Speed)
            end
            
            %% MAXIMUM STEPS REACHED
            if step == maxsteps
                
                %% Final results
                results.t1 = t1;
                results.t2 = t2;
                results.u1 = u1;
                results.u2 = u2;
                results.v1 = u2v(u1);
                results.v2 = u2v(u2);
                
                fitness = step;
                
                [METest1,~,~] = Energy_test(t1(end),u1(end,:),Param,step);
                [METest2,~,~] = Energy_test(t2(end),u2(end,:),Param,step);
                
                fprintf('\nStable solution:\nq1 = %g;\nq2 = %g;\nq3 = %g;\nq4 = %g;\n',q)
                fprintf('q1dot = %g;\nq2dot = %g;\nq3dot = %g;\nq4dot = %g;\n',qdot)
                fprintf('ME is %g/%g\n',METest1,METest2)
                fprintf('Walking Speed is %.3g/%.3g m/s\n',Speed(indice-1),Speed(indice))
                
                
                %% output stable solution
%                 % This is repeated later. Consider cleaning up.
%                 [ Qpost ] = momentumSpringREDUX( u1(end,1),u1(end,2),u1(end,3),0,u1(end,4),u1(end,5),u1(end,6),0);
%                 qdot = [Qpost'];
%                 q = [u1(end,1),u1(end,2),u1(end,3),0];
                
%                 [METest1,~,~] = Energy_test(t1(end),u1(end,:),Param,step);
%                 [METest2,~,~] = Energy_test(t2(end),u2(end,:),Param,step);
%                 
%                 fprintf('\nStable solution:\nq1 = %g;\nq2 = %g;\nq3 = %g;\nq4 = %g;\n',q)
%                 fprintf('q1dot = %g;\nq2dot = %g;\nq3dot = %g;\nq4dot = %g;\n',qdot)
%                 fprintf('ME is %g/%g\n',METest1,METest2)
%                 fprintf('Walking Speed is %.3g/%.3g m/s\n',Speed(step-maxsteps+11),Speed(step-maxsteps+10))
%                 
                if animation == 1 && writeMovie == 1
                    writerObj = VideoWriter('DW.avi');
                    writerObj.FrameRate = 20;     % Must do this before 'open(writeObj)'
                    open(writerObj);
                    writeVideo(writerObj,total_M);
                    close(writerObj);
                end
                
                bb = 1;
                break
                %                                 return
            end % maxsteps reached
            
            %% initial conditions for the next step
            
            % Energy plots
            if energyPlots ==1
                [ME1,KE1,PE1] = Energy_test(t1,u1,Param,step);
                
                figure(2)
                hold on
                plot(t1,ME1,'k')
                plot(t1,PE1,'b')
                plot(t1,KE1,'r')
                
                [ME2,KE2,PE2] = Energy_test(t2,u2,Param,step);
                
                fprintf('ME at start \t = \t %g\nME at end \t\t = \t %g\n',ME2(1),ME2(end))
                
                figure(2)
                hold on
                plot(t2,ME2,'k')
                plot(t2,PE2,'b')
                plot(t2,KE2,'r')
                if abs(ME2(end) - ME_0) > 0.5
                    warning('Mechanical Energy is not conserved')
                end
            end

            
        case {'noTakeOff','noTouchDown','fallen','prematureLiftOff'}
            Speed = nan(11,1);
            fitness = step-1;
            GRFs = NaN;
            results = NaN;
            break
    end
    
end
% For Basin of Attraction?
% save(['Gain1p316_',num2str(bifurcation),'.mat'],'q1Vector','q3Vector','q1dotVector','q3dotVector','basin','Speed')

% save('basinParametersRoll12.mat','q1Vector','q1dotVector','q3dotVector','basin')

end

% function [S1th,S2th,S3th,S4th,xh,yh,xm,ym,xf,yf]= rollshape(s1,s2,s3,s4,rh,rm,rf,h)
% 
% S2th=atan(2/rm*(s2-h));
% S3th=atan(2/rm*(s3-h));
% 
% xm=h;
% ym=-1/rm*h^2;
% xh=tan(S2th)*(rm-rh)/2+xm;
% yh=rm/4*tan(S2th)^2+ym-rh/4*tan(S2th)^2;
% xf=tan(S3th)*(rm-rf)/2+xm;
% yf=rm/4*tan(S3th)^2+ym-rf/4*tan(S3th)^2;
% 
% S1th=atan(2/rh*(s1-xh));
% S4th=atan(2/rf*(s4-xf));
% end
