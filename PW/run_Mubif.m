function run_Mubif
% code used for unbalanced paper.
% gain values

clc
global ini

% To collect the t and u of each step. Takes a lot of data so not 
% recommended for large amounts of data
collectU = 0;

StableMax = 1000;
tic
%VALUES FOUND FROM GOSWAMI PG. 31
Beta=0.6;         % length ratio (b/a)        % default Beta=1,   Mew=2
Mew=3.6;          % mass ratio (mh/m)
M_total=80;
L=1;

%Unbalanced values
mass_pro_ratio = 1;      % (M_pro/M_phy)
length_pro_ratio = 1;    % (c/b) or (pro_upper/phy_upper) 0.73

alpha=-2*pi/180;     %slope angle
maxsteps=100;        % maximum number of steps taken


%% Physiological foot parameters
x_A=0;
y_A=L/(1+Beta);
m_A=M_total/(2+Mew);

%% Prosthetic foot parameters
x_B=0;          % x'2
y_B=L-((L-y_A)*length_pro_ratio);        % y'2
% m_B=m_A*mass_pro_ratio;

%% Other parameters
xc=0;
yc=L;
mc=m_A*Mew;


%%
s1 = -0.4712;      % hindfoot length
s2 = -0.01;       %
s3 = 0.01;
s4 = 0.4712;       % forefoot length

rh = 1.5;
rm = 1.5;
rf = 1.5;

length = 0.25;
ratio = 0.1;

% Don't know what this was used for???
s1_i = -(length*ratio)/(1+ratio);
s4_i = length/(1+ratio);

% Radius of curvature
r_A = 0.3*yc;
r_B = 0.3*yc;

%% PHYSIOLOGICAL FOOT PARAMETERS
s1_A=s1;      % hindfoot length
s2_A=s2;       %
s3_A=s3;
s4_A=s4;       % forefoot length
rh_A=r_A;      % hindfoot gain
rm_A=r_A;      % midfoot gain
rf_A=r_A;      % forefoot gain
h_A=0;            % horizontal distance

%% PROSTHETIC FOOT PARAMETERS
s1_B=s1;      % hindfoot length
s2_B=s2;       %
s3_B=s3;
s4_B=s4;       % forefoot length
rh_B=r_B;      % hindfoot gain
rm_B=r_B;      % midfoot gain
rf_B=r_B;      % forefoot gain
h_B=0;            % horizontal distance



%% initial condition range
%     IC = [  -0.1,-1.2;  % th1
%              0.1,3   ;  % th1dot
%             -12, 6   ]; % th2dot
% th1   % th1dot   % th2dot
% Lb=[-0.1      0.7          -10   ];
% Ub=[-0.5      3          5    ];
Lb=[-0.1      0.1          -2   ];     % -12
Ub=[-1.2      2.5          2    ];        % 6

% Change to x1:0.001:x2 for example, instead of linspace
MuBs = 1.9:0.025:15;  % change to 0.001
Bif = size(MuBs,2);

%% Maintain moment of inertia of the prosthetic leg
Inertia = m_A*(yc-y_A)^2;

%% This is where the bifurcation starts
% First initial condition is randomly selected
nest(1,:)=Lb+(Ub-Lb).*rand(size(Lb));
for bifurcation = 1:Bif
    
    x_axis = MuBs(bifurcation);
    
    m_B = mc/MuBs(bifurcation);
    
    
    % Make it so S1 and S4 are -90deg and 90deg
    s1_B = r_B*-pi/2;
    s4_B = r_B*+pi/2;
    
    % function dictate rollover shape values
    [S1th_A,S2th_A,S3th_A,S4th_A,xh_A,yh_A,xm_A,ym_A,xf_A,yf_A]=rollshape(s1_A,s2_A,s3_A,s4_A,rh_A,rm_A,rf_A,h_A);
    [S1th_B,S2th_B,S3th_B,S4th_B,xh_B,yh_B,xm_B,ym_B,xf_B,yf_B]=rollshape(s1_B,s2_B,s3_B,s4_B,rh_B,rm_B,rf_B,h_B);
    
    Stable = 0;
    
    for i = 1:StableMax
        %% INITIAL THETA VALUES
        th1_0	=	nest(i,1);
        % th2_0	=	golden_th2contact(th1_i);
        th1dot_0=	nest(i,2);
        th2dot_0=	nest(i,3);
        for step=1:maxsteps
            
            if mod(step,2)
                %       odd          disp('foot 1')
                ini=[x_A y_A m_A x_B y_B m_B xc yc mc alpha...
                    rh_A rm_A rf_A h_A xh_A yh_A xm_A ym_A xf_A yf_A S1th_A S2th_A S3th_A S4th_A...
                    rh_B rm_B rf_B h_B xh_B yh_B xm_B ym_B xf_B yf_B S1th_B S2th_B S3th_B S4th_B];
            else
                %      even           disp('foot 2')
                ini=[x_B y_B m_B x_A y_A m_A xc yc mc alpha...
                    rh_B rm_B rf_B h_B xh_B yh_B xm_B ym_B xf_B yf_B S1th_B S2th_B S3th_B S4th_B...
                    rh_A rm_A rf_A h_A xh_A yh_A xm_A ym_A xf_A yf_A S1th_A S2th_A S3th_A S4th_A];
            end
            
            if step == 1
                th2_0	=	golden_th2contact(th1_0);
            end
            
            [t,u,th1sw,th2sw,th1dotsw,th2dotsw,event,stp,Dist]=lagrangian_approach(th1_0,th2_0,th1dot_0,th2dot_0,step);
            
            switch event
                case {'none','fallforward','fallbackward'}
                    StepsReached(i) = step;
                    break
                case 'hlstr'

                    %                     Avg_Vel = Dist / stp;       % Average velocity
                    %                     fprintf('\nConditions for next heel strike, end of step %g:\nStep Period=%g',step,stp)
                    %                     fprintf('s\nth1_0\t=\t%g;\nth2_0\t=\t%g;\nth1dot_0=\t%g;\nth2dot_0=\t%g;\n',th1sw,th2sw,th1dotsw,th2dotsw)
                    % animation(t,u,step)
                    if step == 132
                        plots(t,u,step)
                    end
                    
                    th1_0=th1sw; th2_0=th2sw; th1dot_0=th1dotsw; th2dot_0=th2dotsw;     %   CHANGES AT EACH STEP
                    
                    %         if step == 24
                    %             plots(t,u,heel,step);     % plots graphs
                    %         %fprintf('\nConditions at heel strike:\nth1sw=%g\nth2sw=%g\nth1dotsw=%g\nth2dotsw=%g',th1_0,th2_0,th1dot_0,th2dot_0)
                    %         end
                    
                    %% Collect data for every step
                    % uAll is all t and u for every step. Useful for
                    % phase-plane diagrams but uses a lot of data when
                    % saving
                    if collectU
                        uAll(step,:) = {t,u};
                    end
                    DistAll(step) = Dist;
                    tf(step)=stp;       %time periods of each step
                    inter(step)=th2sw;
                    angvel(step)=th1dotsw;
                    

                    %% Mechanical Energy
                    [KE,PE(step)]=check_energy(t(end),u(end,:));
                    MEall(step) = KE+PE(step);
                    
                    % Calculates energy lost at collision at the end of
                    % step by comparing ME to the next step
                    if step > 1
                        [~,PEnew]=check_energy(t(1),u(1,:));
                        diffE(step-1) = MEall(step) - MEall(step-1) - (PEnew-PE(step-1));
                    end

                    
%                     IC1(step) = th1sw;
%                     IC2(step) = th2sw;
%                     IC3(step) = th1dotsw;
%                     IC4(step) = th2dotsw;
                    %     if step==maxsteps-1
                    %         stepperiod(ii,1)=stp;
                    %     end
                    if step==maxsteps

                        
                        fprintf('\nConditions for next heel strike, end of step %g:\nStep Period=%g',step,stp)
                        fprintf('s\nth1_0\t=\t%g;\nth2_0\t=\t%g;\nth1dot_0=\t%g;\nth2dot_0=\t%g;\n',th1sw,th2sw,th1dotsw,th2dotsw)
                        
                        %% Plot bifurcations, not saved
                        figure(101)
                        hold on
                        h(:,1) = plot(x_axis,tf(end-14:2:end),'r.'); % intact
                        h(:,2) = plot(x_axis,tf(end-15:2:end),'b.'); % prosthetic
                        
                        figure(102)
                        hold on
                        h(:,1) = plot(x_axis,inter(end-14:2:end),'r.'); % intact
                        h(:,2) = plot(x_axis,inter(end-15:2:end),'b.'); % prosthetic
                                     
                        figure(103)
                        hold on
                        h(:,1) = plot(x_axis,DistAll(end-14:2:end),'r.'); % intact
                        h(:,2) = plot(x_axis,DistAll(end-15:2:end),'b.'); % prosthetic
                    
%                         figure(104)
%                         hold on
%                         h(:,1) = plot(r_B,Dist(end-14:2:end)./tf(end-14:2:end),'b.'); % intact
%                         h(:,2) = plot(r_B,Dist(end-15:2:end)./tf(end-14:2:end),'r.'); % prosthetic
                       
                        %% Save all variables
                        diffE(step) = NaN;
                        AllVar(bifurcation,:) = {x_axis, DistAll, tf, inter, angvel, MEall, diffE};
                        
                        if collectU
                            uAllbif(bifurcation,:) = {x_axis,uAll};
                        end
                        % Stability has been reached
                        Stable = 1;
                    end
            end
        end
        if Stable == 1
            % If stable point is found, next bifurcation uses the previous
            % stable initial conditions
            nest(1,:)=nest(i,:);
            break
        end
        % If stability is not found, seect another random set of initial
        % conditions within the range.
        nest(i+1,:)=Lb+(Ub-Lb).*rand(size(Lb));
    end
    
    %% Stability not reached
    
    toc
end
figure(101)
xlabel('Prosthetic foot curvature')
ylabel('Stance time')
% legend(h(1,:),'intact','prosthetic')
hl(1) = plot(NaN,NaN,'b-');
hl(2) = plot(NaN,NaN,'r-');
legend(hl(1,:),'unaffected leg','modified leg')

figure(102)
xlabel('Prosthetic foot curvature')
ylabel('Inter-leg angle')
% legend(h(1,:),'intact','prosthetic')
hl(1) = plot(NaN,NaN,'b-'); % No need to repeat this
hl(2) = plot(NaN,NaN,'r-');
legend(hl(1,:),'unaffected leg','modified leg')

figure(103)
xlabel('Prosthetic foot curvature')
ylabel('Step length')
% legend(h(1,:),'intact','prosthetic')
hl(1) = plot(NaN,NaN,'b-'); % No need to repeat this
hl(2) = plot(NaN,NaN,'r-');
legend(hl(1,:),'unaffected leg','modified leg')

Xlabel = 'Mass ratio of leg B, \mu_B';
Sym = NaN;

save('MuBifSmallV2.mat')

end

function [S1th,S2th,S3th,S4th,xh,yh,xm,ym,xf,yf]= rollshape(s1,s2,s3,s4,rh,rm,rf,h)

S2th=s2/rm;
S3th=s3/rm;

xm=h;
ym=-1/rm*h^2;
xh=tan(S2th)*(rm-rh)/2+xm;
yh=rm/4*tan(S2th)^2+ym-rh/4*tan(S2th)^2;
xf=tan(S3th)*(rm-rf)/2+xm;
yf=rm/4*tan(S3th)^2+ym-rf/4*tan(S3th)^2;

S1th=s1/rm;
S4th=s4/rm;

end
