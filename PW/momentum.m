function Phi=momentum(th1,th2,th1dot,th2dot,Xswp)

phi1=th1+th2;       % th1 post heel strike
phi2=-th2;          % th2 post heel strike

global x1 y1 xc yc m1 mc m2

[X1p,Y1p,Xcp,Ycp,X2p,Y2p,dX1dth1,dXcdth1,dX2dth1,dX2dth2,...
    dY1dth1,dYcdth1,dY2dth1,dY2dth2,~] = coordinates2(th1,th2,1);    % coordiantes pre-strike

[X1,Y1,Xc,Yc,X2,Y2,dX1dphi1,dXcdphi1,dX2dphi1,dX2dphi2,...
    dY1dphi1,dYcdphi1,dY2dphi1,dY2dphi2,sth] = coordinates2(phi1,phi2,2);  % coordiantes post-strike

XcX1pre=(x1-xc)*cos(th1)+(y1-yc)*sin(th1);     %vector from Xc to X1 pre-impact
YcY1pre=-(x1-xc)*sin(th1)+(y1-yc)*cos(th1);         %vector from Yc to Y1 pre-impact

XcX1post=(x1-xc)*cos(phi1)+(y1-yc)*sin(phi1);     %vector from Xc to X1 post-impact     %%CHECK THIS!!!!
YcY1post=-(x1-xc)*sin(phi1)+(y1-yc)*cos(phi1);    %vector from Yc to Y1 post-impact

XcX1test = X1p-Xcp;
YcY1test = Y1p-Ycp;

XcX1post=XcX1pre;
YcY1post=YcY1pre;
% zb=[Xc-sth,X2-sth,X1-sth]
% zz=[Xcp-Xswp,X1p-Xswp,X2p-Xswp]
% Qpre*thdot=Qpost*phidot
% i.e. Qpre*[th1dot,th2dot]'=Qpost*[phi1dot,phi2dot]'

Qpre=[mc*(dYcdth1*(Xcp-Xswp)-dXcdth1*Ycp)+m1*(dY1dth1*(X1p-Xswp)-dX1dth1*Y1p)+m2*(dY2dth1*(X2p-Xswp)-dX2dth1*Y2p),...    %th1dot
    m2*(dY2dth2*(X2p-Xswp)-dX2dth2*Y2p);       %th2dot
    m1*(XcX1pre*dY1dth1-YcY1pre*dX1dth1),   0];
%% +STH REMEMBER POINT OF SWING LEG!!!
Qpost=[mc*(dYcdphi1*(Xc-sth)-dXcdphi1*Yc)+m2*(dY1dphi1*(X1-sth)-dX1dphi1*Y1)+m1*(dY2dphi1*(X2-sth)-dX2dphi1*Y2),...    %phi1dot
    m1*(dY2dphi2*(X2-sth)-dX2dphi2*Y2);          %phi2dot
    m1*(XcX1post*dY2dphi1-YcY1post*dX2dphi1),   m1*(XcX1post*dY2dphi2-YcY1post*dX2dphi2)];     

% H=inv(Qpost)*Qpre;
H =     Qpost\Qpre;
Phi =   H*[th1dot;th2dot];  %Phi = [phi1dot,phi2dot]'


end