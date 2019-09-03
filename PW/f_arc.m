function ds=f_arc(theta,~)
%% Used to numerically integrate to find s(theta_1) i.e. arc length

[~,~,dxth,dyth,~,~]=xth_yth_arc(theta,1);
ds=(dxth^2+dyth^2)^0.5;
end