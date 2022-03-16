function [fg,ax,pl] = FUNC_SingleQuickPlot(x,y,c,cax)
if isempty(c)
    figure,
    pl = plot(x,y,'-');
    set(gca,'YDir','reverse')
else
    cm = getPyPlot_cMap('seismic',100);
    figure,
    pl = imagesc(x,y,c);
    shading 'interp';
    set(gca,'YDir','reverse');
    box(gca,'on')
    caxis(gca,cax);
    colormap(gca,cm);
end

fg = gcf;
ax = gca;
