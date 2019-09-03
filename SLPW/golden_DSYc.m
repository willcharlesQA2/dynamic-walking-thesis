% ------------------------GOLDEN SECTION METHOD----------------------------
% -------------------------------------------------------------------------
% Copyright (c) 2009, Katarzyna Zarnowiec, all rights reserved 
% mailto: katarzyna.zarnowiec@gmail.com
% Adapted for use in double-support equations
% -------------------------------------------------------------------------
function result = golden_DSYc(th1,r1,Param)

% figure; hold on;

a=-th1;                            % start of interval
b=-th1-pi/2;                            % end of interval
epsilon=1e-12;               % accuracy value
iter= 100;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations

% th1	=	-0.253217;

x1=a+(1-tau)*(b-a);             % computing x values
x2=a+tau*(b-a);

f_x1=Simultaneous_DS(th1,r1,x1,Param);                     % computing values in x points
f_x2=Simultaneous_DS(th1,r1,x2,Param);

% plot(x1,f_x1,'rx')              % plotting x
% plot(x2,f_x2,'rx')

while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-tau)*(b-a);
        
        f_x1=Simultaneous_DS(th1,r1,x1,Param);
        f_x2=Simultaneous_DS(th1,r1,x2,Param);
        
%         plot(x1,f_x1,'rx');
    else
        a=x1;
        x1=x2;
        x2=a+tau*(b-a);
        
        f_x1=Simultaneous_DS(th1,r1,x1,Param);
        f_x2=Simultaneous_DS(th1,r1,x2,Param);
        
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


function [Zero] = Simultaneous_DS(th1,q3,th2,Param)    

%% Find the angle made by the swing leg to make contact with the floor with
% a given support leg angle.

LH1 = Param.LH1;
% Xst = Param.Xst;
% S1th = Param.S1th;
% S4th = Param.S4th;
% h = Param.h;

%need these for the roll-over shape
[xth,yth,~,~,~,~]=xth_yth(th1,1,0,Param); 
[xth2,yth2,dxth2,dyth2,ddxth2,ddyth2]=xth_yth(th1+th2,2,0,Param); %   SWING LEG ROLLOVER SHAPE
% xth = 0;
% yth = 0;
% xth2 = 0;
% yth2 = 0;



%% arc length
% [ sphi2, ~,~ ] = arcLength(th1+th2,dxth2,dyth2,ddxth2,ddyth2,S1th,S4th,h,1);
    

    yH1 = LH1 + q3;
% % From here
%     X3 = -xth*cos(th1) + (yH1 - yth)*sin(th1) + sth;
%     Y3 = xth*sin(th1) + (yH1 - yth)*cos(th1);
%     
%     phi2 = th1+th2;
%     
% % Zero = X3/sin(phi2) - Y3/cos(phi2) + (xth2 + xth2*tan(phi2)^2)/tan(phi2) - sphi2/sin(phi2);
% 
% Zero = abs(X3 - Y3*tan(phi2) + xth2*cos(phi2) + xth2*sin(phi2)*tan(phi2) - sphi2 - Xst);
% Zero = abs(zero);
%% New
yH2 = LH1;

xH2 = 0;

xc =  (xth2-xH2)*cos(th2) + (yth2-yH2)*sin(th2);

yc = -(xth2-xH2)*sin(th2) + (yth2-yH2)*cos(th2) + yH1;

% Yc
Yc = -(xc-xth)*sin(th1) + (yc-yth)*cos(th1);
Zero = abs(Yc);
end