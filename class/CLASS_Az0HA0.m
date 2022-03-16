classdef CLASS_Az0HA0
    methods(Static)
%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
function [AxOut] = FabricOrientation_i(PAHV)
    [~,iB] = sort(PAHV,2);
    iB = iB(:,1:2);
    iB = sort(iB,2);
    iB = [iB iB+180];
    iFO = iB(1,2);
    for i = 1:size(iB,1)-1
        for j = 1:4
            stdT(j) = std([iFO ; iB(i+1,j)]);
        end
        [~,im] = min(stdT);
        iFO(i+1,1) = iB(i+1,im);
    end
    if mean(iFO) <= 90
        iFO = [iFO iFO+90];
    else
        iFO = [iFO-90 iFO];
    end
    iFO_Mean = round(mean(iFO),0);
    iFO(iFO>180) = iFO(iFO>180)-180;
    iFO(iFO<1) = iFO(iFO<1)+180;
    iFO_Mean(iFO_Mean>180) = iFO_Mean(iFO_Mean>180)-180;
    iFO_Mean(iFO_Mean<1) = iFO_Mean(iFO_Mean<1)+180;
    FO = iFO-1;
    FO_Mean = iFO_Mean-1;
    AxOut.iFO = iFO;
    AxOut.iFO_Mean = repmat(iFO_Mean,size(iFO,1),1);
    AxOut.FO = FO;
    AxOut.FO_Mean = repmat(FO_Mean,size(iFO,1),1);
end
%%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
function dE0 = HorizontalAnisotropy_i(Z,Zmdl,AX,dEhor)
    psbl_dE = nan(size(Z));
    for i = 1:length(Z)
        psbl_dE(i,1) = dEhor(i,AX(i,1));
    end
    dE1 = abs(psbl_dE);
    dE1(dE1>0.5) = nan; % vertical v3 assumption (v2-v1 can not be larger than 0.5)
    dE1 = fillmissing(dE1,'nearest');
    dE1 = smoothdata(dE1,'movmean');
    j1 = 1;
    for i = 1:size(Zmdl,1)
        [~,j2] = min(abs(Zmdl(i,1)-Z));
        dE0(i,1) = nanmean(dE1(j1:j2,1));
        j1 = j2+1;
    end
end
%%  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    end
end