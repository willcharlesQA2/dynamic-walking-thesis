function [ TP1,FP1 ] = GRFrun( alpha,k,ME,Param )

% Runs the code again to obtain a single GRF from the input and plots a
% graph

[ t_all,GRF_all ] = GRFfromPS( alpha,k,ME,Param );

figure
hold on
tend = -t_all{2-1}.t2(1); % last t value
for i = 2:3; % number of steps
    plot(tend+t_all{i-1}.t2,GRF_all{i-1}.FY2L)
    tend = tend+t_all{i-1}.t2(end);
    plot(tend+t_all{i}.t1,GRF_all{i}.FY1)
    plot(tend+t_all{i}.t2,GRF_all{i}.FY2T)
end

xlabel('time')
ylabel('GRF')
title(['\alpha = ',num2str(alpha*-180/pi),'deg; k = ',num2str(k/1000),'kN/m; ME = ',num2str(ME),'J'])


%% Find GRF of one leg
tend = -t_all{i-1}.t2(1); % last t value
for i = 3; % number of steps
%     plot(tend+t_all{i-1}.t2,GRF_all{i-1}.FY2L)
%     tend = tend+t_all{i-1}.t2(end);
%     plot(tend+t_all{i}.t1,GRF_all{i}.FY1)
%     plot(tend+t_all{i}.t2,GRF_all{i}.FY2T)
tend = -t_all{i-1}.t2(1); % last t value
    TP1 = tend+t_all{i-1}.t2;
    tend = tend+t_all{i-1}.t2(end);
    TP1 = [TP1(1:end-1);tend+t_all{i}.t1(1:end-1);tend+t_all{i}.t2];
    FP1 = [GRF_all{i-1}.FY2L(1:end-1);GRF_all{i}.FY1(1:end-1);GRF_all{i}.FY2T];
end

end