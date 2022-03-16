function cm = FUNC_ZeroCenteredSeismic(mat)
    aaa = linspace(min(mat(:)),max(mat(:)),101);
    cmBlue = getPyPlot_cMap('Blues',length(aaa(aaa<0)));
    cmRed = getPyPlot_cMap('Reds',length(aaa(aaa>0)));
    cm(:,:,1) = [flip(cmBlue);[1 1 1];cmRed];

