function bn = fmcw_findbed(r,a,bsrm,bedmethod,ampthresh)

% Find radar basal reflector depth bed from amplitude profile 
%
% args:
% r = range (m)
% a = amplitude
% bsrm = bed search range (m), 2 element vector
% bedmethod = search medthod
% ampthresh = amplitude threshold (only valid for method ampthresh)

doplot = 0;
ispeak = zeros(size(a));
ispeak(2:end-1) = diff(a(1:end-1))>0 & diff(a(2:end))<0; % all peaks
br1 = r>bsrm(1) & r<bsrm(2);
switch bedmethod
    case 'maxAmp' % choose strongest reflector in range
        [~,bn] = max(abs(a).*br1);
        
    case 'ampThresh' % first peak where amp exceeds ampthresh
        bn = find(db(a)>=ampthresh & r>=bsrm(1) & r<= bsrm(2) & ispeak,1,'first');
        
    case 'xcor' % Match bed return with nominal bed return
        % Make nominal bed return
        dr = diff(r(1:2));
        beddecayscale = 200;
        n = 2*round(beddecayscale/dr);
        x = [-n*dr:dr:n*dr];
        %y = ones(size(x));
        y=exp(-x/beddecayscale);
        y(x<0)=1e-1;
        y = y/1e4;
        
        % xcor to match
        [c,laggs] = xcorr(abs(a(br1)),y);
        [~,mciub] = max(c);
%         cmax=max(c);
%         mciub=find(c>cmax*0.95,1,'first');
        ublag = laggs(mciub);  
        bn1le = find(br1,1,'first') + ublag + (length(x)-1)/2-1; % leading edge of bed
        bn = find(ispeak & r>=r(bn1le),1,'first'); %first peak (or trough) after leading edge
        
        % Plot it
        if doplot
            figure
            plot(r,log10(abs(a)),'r') % actual return
            hold on
            plot(r(find(ispeak)),log10(abs(a(find(ispeak)))),'or')
            plot(r(bn1le)+x,log10(y),'g') % synthetic bed used for match
            plot([r(bn1le) r(bn1le)],log10([1e-10 1e-5]),'col',[0.8 0.8 0.8]) % leading edge of bed
            plot([r(bn) r(bn)],log10([1e-10 1e-5]),'col',[0.6 0.6 0.6]) % bed
            title('a')
        end
end
if isempty(bn)
    error('bed not found')
end