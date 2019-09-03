%% Use the 'Data Cursor' to select a point from the BOA plot and run this 
 % file to get the GRF curve

datacursormode on
h = datacursormode;
s = getCursorInfo(h);

[ t_all,GRF_all ] = GRFfromPS( s.Position(1)*-pi/180,s.Position(2)*1000,s.Position(3),Param );

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
title(['\alpha = ',num2str(s.Position(1)),'; k = ',num2str(s.Position(2)),'; ME = ',num2str(s.Position(3))])


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
figure
hold on

xlabel('time')
ylabel('GRF')
% title(['\alpha = ',num2str(s.Position(1)),'; k = ',num2str(s.Position(2)),'; ME = ',num2str(s.Position(3))])
plot(TP1,FP1)