function [xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(th1,leg)
global ini
%% finds x(theta) and y(theta) for stance leg.
%% can also be used to find contact point of swing leg at heel strike

%,xh,yh,rh,xm,ym,rm,xf,yf,rf,S1th,S2th,S3th,S4th)
% global xh yh rh xm ym rm xf yf rf S1th S2th S3th S4th 
% Physiological foot parameters
% s1=-0.03;      % hindfoot length
% s2=-0.01;       %
% s3=0.011;
% s4=0.17;       % forefoot length
% rh=0.5;      % hindfoot gain
% rm=0.5;      % midfoot gain
% rf=0.5;      % forefoot gain
% h=0;            % horizontal distance
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
 %%
if leg==1       % STANCE LEG
    i=0;
elseif leg==2   % SWING LEG
    i=14;
end
    rh=ini(11+i);
    rm=ini(12+i);
    rf=ini(13+i);
%     h=ini(14+i);
    xh=ini(15+i);
    yh=ini(16+i);
    xm=ini(17+i);
    ym=ini(18+i);
    xf=ini(19+i);
    yf=ini(20+i);
    S1th=ini(21+i);
    S2th=ini(22+i);
    S3th=ini(23+i);
    S4th=ini(24+i);
    
    %% Ignore everything above, this is the curvature
    fr = rm;
     
if th1<=S1th          %start point theta(S1)
    theta=S1th;
    xft=xh; yft=yh; rft=rh;
elseif (S1th<th1) && (th1<S2th)
    theta=th1;
    xft=xh; yft=yh; rft=rh;
elseif (S2th<=th1) && (th1<=S3th)
    theta=th1;
    xft=xm; yft=ym; rft=rm;
elseif (S3th<th1) && (th1<S4th)
    theta=th1;
    xft=xf; yft=yf; rft=rf;
elseif th1>=S4th      %end point theta(S4)
    theta=S4th;
    xft=xf; yft=yf; rft=rf;
else 
    disp(th1)
    error('There is a problem with the rollover function')
end

xth =  fr*sin(theta);
yth = -fr*cos(theta) + fr;
% dxth=sec(theta)^2*(rft/2);
% dyth=dxth*tan(theta);       % dyth=sec(theta)^2*(xth-xm?); may be used instead
% ddxth=(rft)*sec(theta)^2*tan(theta);
% % ddyth=dxth*sec(theta)^2+2*sec(theta)^2*tan(theta)*(xth);
% ddyth=ddxth*tan(theta)+dxth*sec(theta)^2;

    dxth = fr*cos(theta);
    dyth = fr*sin(theta);

    ddxth = -fr*sin(theta);
    ddyth =  fr*cos(theta);




end