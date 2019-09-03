function ds=f_arc(theta,~)
%% Used to numerically integrate to find s(theta_1) i.e. arc length
global Param
[~,~,dxth,dyth,~,~]=xth_yth(theta,1,1,Param);
ds=(dxth^2+dyth^2)^0.5;
end

