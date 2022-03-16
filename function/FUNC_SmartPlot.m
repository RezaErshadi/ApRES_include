function [fg,ax] = FUNC_SmartPlot()
pltdim = [0.5000    0.0256    0.5003    0.9213];
CLASS_FixedPlot.SetFigureSize(pltdim(1),pltdim(2),pltdim(3),pltdim(4));    
fg = gcf;   
fg.Color = "white";
if figvis ~= 1
    set(fg,'Visible','off');
end