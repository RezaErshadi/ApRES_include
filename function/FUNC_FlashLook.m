function [] = FUNC_FlashLook(filepath,p,att,material,maxRange,tdc,BedRange)
%%
DtaOrg = fmcw_load(filepath,1);
[data,f,z,twt,~] = FUNC_readApRESfile(filepath,p,att,material,tdc);
t = DtaOrg.t;
vlt = DtaOrg.vif;
amp = abs(data{end});
if ~isempty(BedRange)
    iBed = fmcw_findbed(z,amp,BedRange,'xcor');
    IceThickness = z(iBed);
    ThicknessTime = t(iBed);
    fprintf("Ice Thickness = %f at %f seconds\n",IceThickness,ThicknessTime)
end
[~,imr] = min(abs(maxRange-z));
z(imr+1:end) = [];
amp(imr+1:end) = [];
pwr = 20.*log10(amp);
% ---
figure,
subplot(1,2,1)
plot(vlt,t)
hold on
plot([0 0],[min(t) max(t)],'-r','linewidth',2)
plot([2.5 2.5],[min(t) max(t)],'-r','linewidth',2)
plot([0.5 0.5],[min(t) max(t)],'-g','linewidth',2)
plot([2 2],[min(t) max(t)],'-g','linewidth',2)
xlim([-0.5 3])
xlabel("Voltage")
ylabel("Time [s]")
set(gca,'FontSize',16)
%---
subplot(1,2,2)
yyaxis left
plot(pwr,z,'.-')
if ~isempty(BedRange)
    hold on
    plot([min(pwr)-10 max(pwr)+10],[IceThickness IceThickness],'-k','linewidth',2)
end
xlim([min(pwr)-10 max(pwr)+10])
ylim([0 maxRange])
set(gca,'YDIR','reverse')
xlabel("Power [dB]")
ylabel("Depth [m]")
yyaxis right
ylabel("Time [s]")
% plot(pwr,t(1:imr),'.-')
% xlim([min(pwr)-10 max(pwr)+10])
ylim([0 t(imr)])
set(gca,'YDIR','reverse')
grid on
set(gca,'FontSize',16)
end