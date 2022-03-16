function signal = FUNC_SortFullAzimuth(Signal,AntOr,AntAng,or)
sig_i = AntOr == or;
ORsignal = Signal(:,sig_i);
sig_or = AntAng(sig_i);
[~,sort_i] = sort(sig_or);
signal = ORsignal(:,sort_i);
signal(:,end) = [];