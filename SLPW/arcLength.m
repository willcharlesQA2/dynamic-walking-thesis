function [ sth, dsth, ddsth ] = arcLength(q1,calculateSth,leg,Param)
%arcLength calculates sth, dsth, ddsth
%   Detailed explanation goes here
if leg == 1
    fr = Param.fr1;
    S1th = Param.S1th;
    S4th = Param.S4th;
elseif leg == 2
    fr = Param.fr2;
    S1th = Param.S1th2;
    S4th = Param.S4th2;
end


%%CALCULATE ARC LENGTH
if q1<=S1th          %start point theta(S1)
    th_arc=S1th;
elseif q1>=S4th
    th_arc=S4th;      %end point theta(S4)
else
    th_arc=q1;
end

if calculateSth == 1
    sth = fr*th_arc;
else
    sth = NaN;
end

% sth = Arc length i.e. s_{\theta}
if th_arc == S1th || th_arc==S4th
    dsth=0; ddsth=0;
else
    dsth=fr;
    ddsth=0;
end
end

