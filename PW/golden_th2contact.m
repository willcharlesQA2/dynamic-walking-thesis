% ------------------------GOLDEN SECTION METHOD----------------------------
% -------------------------------------------------------------------------
% Copyright (c) 2009, Katarzyna Zarnowiec, all rights reserved 
% mailto: katarzyna.zarnowiec@gmail.com
% -------------------------------------------------------------------------
function result = golden_th2contact(th1)

% figure; hold on;

a=0-th1;                            % start of interval (was 0.4)
b=pi/2-th1;                            % end of interval (was 2)
epsilon=1e-12;               % accuracy value
iter= 100;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations

% th1	=	-0.253217;

x1=a+(1-tau)*(b-a);             % computing x values
x2=a+tau*(b-a);

f_x1=th2contact(th1,x1);                     % computing values in x points
f_x2=th2contact(th1,x2);

% plot(x1,f_x1,'rx')              % plotting x
% plot(x2,f_x2,'rx')

while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        
        f_x1=th2contact(th1,x1);
        f_x2=th2contact(th1,x2);
        
%         plot(x1,f_x1,'rx');
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        
        f_x1=th2contact(th1,x1);
        f_x2=th2contact(th1,x2);
        
%         plot(x2,f_x2,'rx')
    end
    
    k=k+1;
end


% chooses minimum point
if(f_x1<f_x2)
%     sprintf('x_min=%f', x1)
%     sprintf('f(x_min)=%f ', f_x1)
%     plot(x1,f_x1,'ro')
    result = x1;
else
%     sprintf('x_min=%f', x2)
%     sprintf('f(x_min)=%f ', f_x2)
%     plot(x2,f_x2,'ro')
    result = x2;
end
end


function [Y_abs] = th2contact(th1,th2)    

%% Find the angle made by the swing leg to make contact with the floor with
% a given support leg angle.
global ini

xc = ini(7);
yc = ini(8);

[xth,yth,~,~,~,~]=xth_yth(th1,1); 
[xth2,yth2,~,~,~,~]=xth_yth(th1+th2,2); %   SWING LEG ROLLOVER SHAPE
xrp=(xth2-xc)*cos(th2)+(yth2-yc)*sin(th2)+xc;  % contact point in [x,y]
yrp=-(xth2-xc)*sin(th2)+(yth2-yc)*cos(th2)+yc; % reference frames.

%Xswp=(xrp-xth)*cos(th1)+(yrp-yth)*sin(th1)+sth; %X location of swing leg at contact point %CHANGE STH
Yswp=-(xrp-xth)*sin(th1)+(yrp-yth)*cos(th1);      %Y location of swing leg at contact point

Y_abs = abs(Yswp);

end