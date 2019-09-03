function animation_gif(total_M,filename)
% Creates a gif file using the frame obtaine by get_frames and filename
% e.g. 'testAnimated.gif'
DelayTime = 1/20;

ntot = size(total_M);
for n = 1:ntot(2)
    im = frame2im(total_M(n));
    [A,map] = rgb2ind(im,256);
    if n == 1;
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',DelayTime);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',DelayTime);
    end
end

end