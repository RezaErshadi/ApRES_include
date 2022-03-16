function [CrcY] = FUNC_Depth2Elevation(Y,SE)
    NS = size(Y,2);
    CrcY = nan(size(Y,1),NS);
    for i =1:NS
        CrcY(1,i) = SE(1,i) - Y(1,i);
        a = nan(size(Y,1),2);
        a(:,1) = Y(:,i);
        dff = diff(Y(:,i));
        for j = 1:length(dff)
            CrcY(j+1,i) = CrcY(j,i)-dff(j);
        end
    end
end