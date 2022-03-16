function [Xmat,Ymat,Pmat] = FUNC_Column2MatrixProfile(CHL,CX,CY,PAR)
% CHL: Column Horizontal Length
% CX: Column X
% CY: Column Y
% PAR: Column Parameter
    x0 = sort([CX-(CHL/2) CX+(CHL/2)]);
    x1 = repmat(x0,size(CY,1),1);
    Xmat = nan(size(x1,1),1);
    Ymat = nan(size(x1,1),1);
    Pmat = nan(size(x1,1),1);
    NS = length(CX); % number of radar sites
    for j = 1:NS
        j1 = (j*2)-1;
        j2 = (j*2);
        Xmat = [Xmat x1(:,j1:j2) nan(size(x1,1),1)];
        Ymat = [Ymat repmat(CY(:,j),1,2) nan(size(x1,1),1)];
        Pmat = [Pmat repmat(PAR(:,j),1,2) nan(size(x1,1),1)];
    end
end