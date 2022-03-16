% Fixed plot and relevant plotting functions
classdef CLASS_FixedPlot
methods(Static)
% Plot types:
% 1-Power Return (HH,VH,HV,VV) --> 2
% 2-Power Anomaly (HH,VH,HV,VV) --> 1
% 3-Phase Difference (HH,VH,HV,VV) --> 1
% 4-Coherence Phase (HHVV) --> 1
% 5-Coherence Magnitude (HHVV) --> 1
% 6-Coherence Real & Imaginary Part (HHVV) --> 1
% 7-Coherence Phase Derivative (HHVV) --> 1
% 8-Coherence Scaled Phase Derivative (HHVV) --> 1
%%
function[fg,ax,cb] = AdvancePlot(fg,Dta,ao,Z,w2plt,bd,pltdim)
    % An advance function to plot any kind of radar data
    if isempty(fg) == 1
        CLASS_FixedPlot.SetFigureSize(pltdim(1),pltdim(2),pltdim(3),pltdim(4));    fg = gcf;   fg.Color = "white";
    end
%      dtaParams = {1:HH0,   2-VH0,  3-HV0,  4-VV0,...
%                   5:PA_HH,  6:PA_VH,    7:PA_HV,    8:PA_VV,...
%                   9:PHI_HH, 10:PHI_VH,  11:PHI_HV, 12:PHI_VV,...
%                   13:absC, 14:argC, 15:reC, 16: imC,
%                   17:dphi_dz, 18:dEhor};
%------------------ Ax Title
    ttl = { 'PR(HH)','PR(VH)','PR(HV)','PR(VV)',...
            'P_{HH}','P_{VH}','P_{HV}','P_{VV}',...
            "d\phi_{HH}/d\theta'","d\phi_{VH}/d\theta'","d\_{HV}/d\theta'","d\_{VV}/d\theta'",...
            '|C_{HHVV}|','\phi(C_{HHVV})','re(C_{HHVV}|','im(C_{HHVV})',...
            'd\phi(C_{HHVV})/dZ','\Psi'};
%------------------ Colorbar Title
    cbttl = {   '[-]','[-]','[-]','[-]',... 
                '[dB]','[dB]','[dB]','[dB]',...
                '[rad]','[rad]','[rad]','[rad]',...
                '[-]','[rad]','[-]', '[-]',...
                '[rad/m]','[-]'};
%------------------ Colorbar Limit
    cax1 = [];
    cax2 = [-5 5];
    cax3 = 0.1.*[-pi pi];
    cax4 = [0 1];
    cax5 = [-pi pi];
    cax6 = 0.02.*[-pi pi];
    cax7 = [0 1];
    cax8 = [-1 1];
    cax =  {cax1,cax1,cax1,cax1,... 
            cax2,cax2,cax2,cax2,...
            cax3,cax3,cax3,cax3,...
            cax4,cax5,cax8,cax8,... 
            cax6,cax7};
%------------------         
    fntsz = [16,14]; % Title and Other font sizes
    xlbl = "\alpha [deg]"; % X lable
    ylbl = "Depth [m]"; % Y label
    tpm=(ao(1:end-1)+ao(2:end))/2; % X center
    Zm=(Z(1:end-1)+Z(2:end))/2; % Y center
%------------------ Number of Rows and Columns
    sbpltrow = [[1:18];
                [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,3]];
    nsf = length(w2plt); % Number of SubPlots
    rw = sbpltrow(nsf,2);
    col = ceil(nsf/rw);
%------------------
    for k = 1:length(w2plt)
        i = w2plt(k); % parameter flag
        ax{k} = subplot(rw,col,k);
        cDta = Dta{i}; % current parameter data
        cm = getPyPlot_cMap('seismic',100); % default colorbar (RedBlue Python Seismic)
        plttyp = 1; % Default Plot Type
        if i <= 4 
            plttyp = 2; % Plot type for Power Return
        end
        if i == 18
            ll = -0.01; % girdle strength lower limit
            cmlength = length(-1:abs(ll):1)-1;
            dE = cDta;
            dE(dE<ll) = nan;
            dE(dE>1) = nan;    
            cDta = dE;
            cm1 = getPyPlot_cMap('tab20b',(cmlength)/2);
            cm2 = getPyPlot_cMap('bone',(cmlength)/2);
            cm = [cm1;cm2];
        end
        if size(cDta,2) == length(ao);x = ao;else;x = tpm;end
        if size(cDta,1) == length(Z);y = Z;else;y = Zm;end
        [cb{k},ax{k}] = CLASS_FixedPlot.MakeAdvancePlot...
            (ax{k},ao,Z,cDta,ttl{i},xlbl,ylbl,cbttl{i},cax{i},cm,fntsz,plttyp,bd);
    end
    fg.InvertHardcopy = 'off';
end  
%%
function [fg,ax,cb] = StandardFigure(fg,Dta,ao,Z,bd,pltdim)
    if isempty(fg) == 1
        CLASS_FixedPlot.SetFigureSize(pltdim(1),pltdim(2),pltdim(3),pltdim(4));    fg = gcf;   fg.Color = "white";
    end 
    w2plt = [1,5,7,13,14,18];
    ax{1} = subplot(2,3,1);   ax{2} = subplot(2,3,2);   ax{3} = subplot(2,3,3);
    ax{4} = subplot(2,3,4);   ax{5} = subplot(2,3,5);   ax{6} = subplot(2,3,6);
%------------------ Ax Title
    ttl = { 'Power Return',...
            'Power Anomaly (HH)','Power Anomaly (VH)',...
            'Magnitude (C_{HHVV})','Phase (C_{HHVV})',...
            'Scaled Derivative (C_{HHVV})'};
%------------------ Colorbar Title
    cbttl = {   '[-]',... 
                'P_{HH} [dB]','P_{HV} [dB]',...
                'abs(C_{HHVV} [-]','arg(C_{HHVV}) [rad]',...
                '\Psi [-]'};
%------------------ Colorbar Limit
    cax1 = [];
    cax2 = [-5 5];
    cax3 = [0 1];
    cax4 = [-pi pi];
    cax5 = [0 1];
    cax =  {cax1,... 
            cax2,cax2,...
            cax3,cax4,... 
            cax5};
%------------------         
    fntsz = [16,14]; % Title and Other font sizes
    xlbl = "\alpha [deg]"; % X lable
    ylbl = "Depth [m]"; % Y label
    tpm=(ao(1:end-1)+ao(2:end))/2; % X center
    Zm=(Z(1:end-1)+Z(2:end))/2; % Y center
%------------------ Number of Rows and Columns
    rw = 2;
    col = 3;
%------------------
    PR = 20.*log10(abs([Dta{1}(:,1) Dta{2}(:,1) Dta{3}(:,1) Dta{4}(:,1)])); % Power Return
%     PR = PR - mean(mean(PR));
%     PR = [PR(:,1)-mean(PR(:,1)) PR(:,2)-mean(PR(:,2)) PR(:,3)-mean(PR(:,3)) PR(:,4)-mean(PR(:,4))];
    k = 1;
    ax{1} = subplot(rw,col,k);
    for i = 1:4
        plot(ax{1},PR(:,i),Z,'.')
        hold(ax{1},'on')
    end
    if ~isempty(bd)
        plot(ax{1},xlim(ax{1}),[bd bd],'-w','LineWidth',2)
        plot(ax{1},xlim(ax{1}),[bd bd],'--k','LineWidth',2)
    end
    set(ax{1},'YDir','reverse')
    grid(ax{1},'on')
    xlim(ax{1},[min(PR(:)) max(PR(:))])
    xlabel(ax{1},"Power Return [dB]");
    ylabel(ax{1},ylbl);
    ylim(ax{1},[Z(1) Z(end)]);
    set(ax{1},'FontSize',fntsz(2),'Layer','top');
    title(ax{1},ttl{k},'FontSize',fntsz(1))
    
    for k = 2:length(w2plt)
        i = w2plt(k); % parameter flag
        ax{k} = subplot(rw,col,k);
        cDta = Dta{i}; % current parameter data
        cm = getPyPlot_cMap('seismic',100); % default colorbar (RedBlue Python Seismic)
        plttyp = 1; % Default Plot Type
        if i <= 4 
            plttyp = 2; % Plot type for Power Return
        end
        if i == 13
            cmlength = length(0:0.01:1)-1;
            cmRed = 2*round(0.4*cmlength,0);
            cmBlue = 2*cmlength-cmRed;
            cmB = getPyPlot_cMap('seismic',cmBlue);
            cmR = getPyPlot_cMap('seismic',cmRed);
            cm1 = cmB(1:cmBlue/2,:);
            cm2 = cmR(cmRed/2+1:end,:);
            cm = flip([cm1;cm2]);
        end
        if i == 18
            ll = -0.01; % girdle strength lower limit
            cmlength = length(-1:abs(ll):1)-1;
            dE = cDta;
            dE(dE<ll) = nan;
            dE(dE>1) = nan;    
            cDta = dE;
            cm1 = flip(getPyPlot_cMap('tab20b',(cmlength)/2));
            cm2 = getPyPlot_cMap('Greys',cmlength);
            cm2 = cm2((cmlength)/2+1:end,:);
            cm = [cm1;cm2];
        end
        if size(cDta,2) == length(ao);x = ao;else;x = tpm;end
        if size(cDta,1) == length(Z);y = Z;else;y = Zm;end
        [cb{k},ax{k}] = CLASS_FixedPlot.MakeAdvancePlot...
            (ax{k},ao,Z,cDta,ttl{k},xlbl,ylbl,cbttl{k},cax{k},cm,fntsz,plttyp,bd);
    end
    for i = [2 3 5 6]
        yticklabels(ax{i},'')
        ylabel(ax{i},[])
    end
    fg.InvertHardcopy = 'off';
end
%%
function[cb,ax] = MakeAdvancePlot(ax,ao,Z,dta,ttl,xlbl,ylbl,cbttl,cax,cm,fntsz,plttyp,bd)
    imAlpha=ones(size(dta));
    imAlpha(isnan(dta))=0;
    switch plttyp
        case 1
            imagesc(ax,ao,Z,dta,'AlphaData',imAlpha)
            shading(ax,'interp');
            set(ax,'YDir','reverse');
            box(ax,'on')
            xtl = [0 45 90 135 180 225 270 315 360];
            xticks(ax,xtl);
            cb=colorbar;
            cb.Label.String = cbttl;
            caxis(ax,cax);
            colormap(ax,cm);
            xlabel(ax,xlbl);
            ylabel(ax,ylbl);
            title(ttl,'FontSize',fntsz(1))
            cb.FontSize  = fntsz(2);
            set(ax,'FontSize',fntsz(2),'Layer','top');
        case 2
            plot(ax,20*log10(abs(dta)),Z) 
            set(ax,'YDir','reverse');
            xlabel(ax,"Power Return [dB]");
            ylabel(ax,ylbl);
            ylim([Z(1) Z(end)]);
            set(ax,'FontSize',fntsz(2),'Layer','top');
            title(ttl,'FontSize',fntsz(1))
            cb = [];
    end
    grid(ax,'on');
    hold(ax,'on')
    if isempty(bd)~=1
        plot(ax,xlim(ax),[bd bd],'-w','LineWidth',2)
        plot(ax,xlim(ax),[bd bd],'--k','LineWidth',2)
%     else
%         plot(ax,xlim(ax),[Z(end) Z(end)],'--k','LineWidth',2)
    end
    hold(ax,'off')
end  
%%
function f = SetFigureSize(ss1,ss2,w,h)
    f = figure;
    set(f,'Color',[1 1 1]);
    set(f, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]); % full screen figure
    set(f, 'Units', 'centimeters');
    scrn = get(f, 'OuterPosition'); % get the size of the screen in CM

    wdt = scrn(3) * w;
    hgt = scrn(4) * h;

    s1 = scrn(3) * ss1;
    s2 = scrn(4) * ss2;

    set(f, 'OuterPosition', [s1, s2, wdt, hgt]); % change the figure size to the new size
    set(f, 'Units', 'Normalized');
end
%%
function c = redblue(m)
    if nargin < 1, m = size(get(gcf,'colormap'),1); end
    if (mod(m,2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        m1 = m*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        r = [r; ones(m1,1)];
        g = [g; flipud(g)];
        b = flipud(r);
    else
        % From [0 0 1] to [1 1 1] to [1 0 0];
        m1 = floor(m*0.5);
        r = (0:m1-1)'/max(m1,1);
        g = r;
        r = [r; ones(m1+1,1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end
    c = [r g b];
end
%% 
function [grad,im]=colorGradient(c1,c2,depth)
    error(nargchk(2, 3, nargin));
    %If c1 or c2 is not a valid RGB vector return an error.
    if numel(c1)~=3
        error('color c1 is not a valir RGB vector');
    end
    if numel(c2)~=3
        error('color c2 is not a valir RGB vector');
    end
    if max(c1)>1&&max(c1)<=255
        %warn if RGB values are given instead of Intensity values. Convert and
        %keep procesing.
        warning('color c1 is not given as intensity values. Trying to convert');
        c1=c1./255;
    elseif max(c1)>255||min(c1)<0
        error('C1 RGB values are not valid.')
    end
    if max(c2)>1&&max(c2)<=255
        %warn if RGB values are given instead of Intensity values. Convert and
        %keep procesing.
        warning('color c2 is not given as intensity values. Trying to convert');
        c2=c2./255;
    elseif max(c2)>255||min(c2)<0
        error('C2 RGB values are not valid.')
    end
    %default depth is 64 colors. Just in case we did not define that argument.
    if nargin < 3
        depth=64;
    end
    %determine increment step for each color channel.
    dr=(c2(1)-c1(1))/(depth-1);
    dg=(c2(2)-c1(2))/(depth-1);
    db=(c2(3)-c1(3))/(depth-1);
    %initialize gradient matrix.
    grad=zeros(depth,3);
    %initialize matrix for each color. Needed for the image. Size 20*depth.
    r=zeros(20,depth);
    g=zeros(20,depth);
    b=zeros(20,depth);
    %for each color step, increase/reduce the value of Intensity data.
    for j=1:depth
        grad(j,1)=c1(1)+dr*(j-1);
        grad(j,2)=c1(2)+dg*(j-1);
        grad(j,3)=c1(3)+db*(j-1);
        r(:,j)=grad(j,1);
        g(:,j)=grad(j,2);
        b(:,j)=grad(j,3);
    end
    %merge R G B matrix and obtain our image.
    im=cat(3,r,g,b);
end
%% 
    end
end