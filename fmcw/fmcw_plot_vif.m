function fmcw_plot_vif(prof)

if isstring(prof) || ischar(prof)
    prof = fmcw_load(prof);
end

mv = mean(prof.vif, 1);
N = size(prof.vif, 1);
rmse = zeros(1, N);
stde = rmse;

Nr = ceil(N/5);

for k = 1:size(prof.vif, 1)
subplot(Nr+1,5,k)
hold off
plot(prof.t, prof.vif(k,:), 'b')    
hold on
plot(prof.t, mv, 'k')
ylim([0 2.5])
plot(prof.t, abs(mv - prof.vif(k,:)), 'r')
rmse(k) = rms(mv-prof.vif(k,:));
stde(k) = std(prof.vif(k,:));
title(sprintf("CHIRP #%d - RMSE: %f", k, rmse(k)))
end

subplot(Nr+1,1,Nr+1)
hold off
stem(rmse, 'rx')
title('RME Error')

