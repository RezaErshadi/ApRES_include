function DtaMean = FUNC_SimpleRead(filePath,material)
% ps = filesep;
% filePath = strcat(DataList.folder,ps,DataList.name);
%%
p=1; % pad factor (i.e. level of interpolation to use during fft)
winFun=@blackman; % window function handle
frange=[2e8 , 4e8]; %Normally [2e8,4e8]
%% Density From field observations
rhos=407.0613; % surface density (kg/m^3)
rhoi=907.7165; % ice density (kg/m^3)
Lrho=39.5512;  % half depth decay (m) (assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho))
nI=1.68;
n=3;
%%
% vif = [];
% chirpNum = [];
% chirpAtt = [];
% chirpTime = [];
% for i = 1:length(filePath)
%     TempDta(i) = fmcw_load(filePath(i),1);
%     TempDta(i) = fmcw_cull_freq(TempDta(i),frange);
%     t(i,1) = datetime(TempDta(i).TimeStamp,'ConvertFrom','datenum');
%     DtaLoad.vif = [vif ; TempDta.vif];
%     DtaLoad.chirpNum = [chirpNum ; TempDta.chirpNum];
%     DtaLoad.chirpAtt = [chirpAtt ; TempDta.chirpAtt];
%     DtaLoad.chirpTime = [chirpTime ; TempDta.chirpTime];
%     DtaLoad.ChirpsInBurst = size(DtaLoad.vif,1); 
%     
%     vif = DtaLoad.vif;
%     chirpNum = DtaLoad.chirpNum;
%     chirpAtt = DtaLoad.chirpAtt;
%     chirpTime = DtaLoad.chirpTime;
% end

DtaLoad = fmcw_load(filePath,1);
DtaLoad = fmcw_burst_split_by_att(DtaLoad);
DtaLoad = DtaLoad(1);
DtaLoad = fmcw_cull_freq(DtaLoad,frange);
DtaMean = fmcw_burst_mean(DtaLoad);
DtaMean.TIME = datetime(DtaLoad.TimeStamp,'ConvertFrom','datenum');
if material == "ice"
    maxRange = 4000;
    [Range,Rfine,SpecCor,spec] = fmcw_range(DtaMean,p,maxRange,winFun);
%     ShouldBeZero=@(d,dI,L,RhoSp) -dI+d+L*(nI-1)/nI*(1-RhoSp)*(exp(-d/L)-1);
%     TrueDepthfun=@(RhoSp,L,dI) FUNC_bisection(@(d) ShouldBeZero(d,dI,L,RhoSp),0,max(dI)+20);
%     TrueDepth = TrueDepthfun(rhos/rhoi,Lrho,Range);
    TrueDepth = Range;
elseif material == "soil"
    maxRange = 50;
    [Range,~,SpecCor,~] = fmcw_range(DtaMean,p,maxRange*n/sqrt(DtaMean.er)+100,@blackman);
    TrueDepth = Range*sqrt(DtaMean.er)/n;
    DtaMean.WaveLength = DtaMean.lambdac*sqrt(DtaMean.er)/n;
end
TrueDepth(1) = 1e-20;
DtaMean.Z = TrueDepth';
DtaMean.Signal = SpecCor.';