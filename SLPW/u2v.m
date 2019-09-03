function [v] = u2v(u)
% Changes the swing leg from th2 to phi2
% Single support
if size(u,2) == 6
    v = [u(:,1),u(:,1)+u(:,2),u(:,3),u(:,4),u(:,4)+u(:,5),u(:,6)];
end
if size(u,2) == 8
    v = [u(:,1),u(:,1)+u(:,2),u(:,3),u(:,4),u(:,5),u(:,5)+u(:,6),u(:,7),u(:,8)];
end