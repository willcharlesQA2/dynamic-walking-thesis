%% Use the 'Data Cursor' to select a point from the BOA plot and run this 
 % file to get the GRF curve

datacursormode on
h = datacursormode;
s = getCursorInfo(h);

[ t_all,GRF_all ] = GRFfromBOA( s.Position(1)*-pi/180,s.Position(2)*1000,s.Position(3),Param );

figure
hold on
plot(t_all{1},GRF_all{1})
plot(t_all{2}+t_all{1}(end),GRF_all{2})

xlabel('time')
ylabel('GRF')
title(['\alpha = ',num2str(s.Position(1)),'; k = ',num2str(s.Position(2)),'; ME = ',num2str(s.Position(3))])