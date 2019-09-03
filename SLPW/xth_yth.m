function [xth,yth,dxth,dyth,ddxth,ddyth]=xth_yth(th1,leg,arc,Param)
%% finds x(theta) and y(theta) for stance leg. At the start and end of the 
% foot, rolling is LOCKED (dxth,dyth = 0), rather than xth_yth_arc.
% can also be used to find contact point of swing leg at heel strike


 %%

 if leg == 1
     fr = Param.fr1;
 else
     fr = Param.fr2;
 end
 
%     h=ini(14+i);

    S1th=Param.S1th;
    S2th=Param.S2th;
    S3th=Param.S3th;
    S4th=Param.S4th;
    
    
     
if th1<=S1th          %start point theta(S1)
    theta=S1th;
elseif (S1th<th1) && (th1<S2th)
    theta=th1;
elseif (S2th<=th1) && (th1<=S3th)
    theta=th1;
elseif (S3th<th1) && (th1<S4th)
    theta=th1;
elseif th1>=S4th      %end point theta(S4)
    theta=S4th;
else 
    disp(th1)
    error('There is a problem with the rollover function')
    % if error occurs check ddsth. i.e. if dxth=0, ddsth=inf??
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

if arc == 0
    if theta == S1th || theta== S4th
        dxth=0; dyth=0; ddxth=0;  ddyth=0;
    end
end

end