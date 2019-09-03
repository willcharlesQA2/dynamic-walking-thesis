function [TP,FP,pks,locs,FPx] = GRFplot2(GRF_all,Bifs)


% F1 = Stance leg   SS
% F2 = Trailing leg DS
% F3 = Lead leg     DS

if nargin == 1
    Bifs = 5;
end
% pks = NaN(Bifs-1,2);
% locs = pks;
for loopF = 1:Bifs-1 % number of steps. Stars at 3 so that odd/even bifurcations are not messed up.
    tend = -GRF_all{loopF}.t2(1); % last t value
    TP1 = tend+GRF_all{loopF}.t2;
%     tend = GRF_all{loopF}.t2(end);
    TP{loopF} = [TP1(1:end-1);tend+GRF_all{loopF}.t1(1:end-1);tend+GRF_all{loopF+1}.t2];
    FP{loopF} = [GRF_all{loopF}.FY2L(1:end-1);GRF_all{loopF}.FY1(1:end-1);GRF_all{loopF+1}.FY2T];
    
    FPx{loopF} = [GRF_all{loopF}.FX2L(1:end-1);GRF_all{loopF}.FX1(1:end-1);GRF_all{loopF+1}.FX2T];
%                                 plot(TP1,FP1)
    [pksTemp,locs] = findpeaks(FP{loopF},TP{loopF});
    
    %% Test to see if pks is the same size
    if loopF > 1 
        if size(pksTemp) == size(pks{1})
            pks{loopF} = pksTemp;
        else
            warning('Different number of peaks for each step')
            pks{loopF} = NaN(size(pks{loopF-1}));
        end
    else
        pks{loopF} = pksTemp;
    end
%     pk1(loopF,:) = pks(1);
%     pk2(loopF,:) = pks(2);
%     loc1(loopF,:) = locs(1);
%     loc2(loopF,:) = locs(2);
%     GRFPeak1(loopF) = pks(1);
%     GRFPeak2(loopF) = pks(2);
end

% %% MUST BE SYMMETRIC WALKING
% LW = 3;
% figure
% hold on
% 
% %% Grey double-support phase
% DSt1 = t2(1)-t2(1);
% SSt1 = t2(end)-t2(1);
% DSt2 = t1(end)-t2(1);
% SSt2 = t1(end) + t2(end) - 2*t2(1);
% fill([DSt1,DSt1,SSt1,SSt1],[ylim fliplr(ylim)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3)
% fill([DSt2,DSt2,SSt2,SSt2],[ylim fliplr(ylim)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.3)
% 
% xlabel('Time (s)')
% ylabel('GRF %BW')

end