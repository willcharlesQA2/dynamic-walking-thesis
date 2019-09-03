function out = f(th,X,Y,Param)
out = (Y-Param.fr)*tan(th) + Param.IC + Param.fr*th - X;
end

% fun = @f;
% x0 = [-1,2];
% z = fzero(fun,x0);