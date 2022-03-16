function [OutputParameters,SiteNumber,fg,ax,cb] = FUNC_QuickStandardPlot(InputParameters,SiteName)
%% Unpack the input parameters
fields = fieldnames(InputParameters);
for i = 1:length(fields)
    eval(strcat(fields{i},' = InputParameters.',fields{i},';'))
end
%% Read The Data
DtaDir = strcat(PP,ps,ProjectName,ps,SiteName);
[hh,vh,hv,vv,Z,dZ,maz,TimeStamp,f,Bed] = FUNC_ReadApRES(DtaDir,maxRange,BedRange);
disp(strcat('loading time: ',string(round(toc,2)),' [s]'))
%%
SiteInfo = FUNC_ReadSiteInfo(InfoDir,ps,ProjectName,SiteName,TimeStamp,Bed);
SiteNumber = SiteInfo.SiteNumber;
AntRot = SiteInfo.AntennaCorrection;
sf = [SiteInfo.HH SiteInfo.VH SiteInfo.HV SiteInfo.VV];
%% Corrections
% Keeping the original signals
ORGhh = hh;
ORGvh = vh;
ORGhv = hv;
ORGvv = vv;
% --- applying Synthesizing Factors
hh = sf(1).*ORGhh;
vh = sf(2).*ORGvh;
hv = sf(3).*ORGhv;
vv = sf(4).*ORGvv;
% ---
%% Signal to Parameters
[HH,VH,HV,VV] = CLASS_S2P.AzimuthSynthesizer(hh,vh,hv,vv,ao,AntRot);
ObsDta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
PAHH = ObsDta{5};
PAHV = ObsDta{7};
CP = ObsDta{14};
Psi = ObsDta{18};
%% Quick Extraction: Fabric Orientation & Horizontal Anisotropy
Zmdl = (ResZmdl:ResZmdl:Z(end))'; Zmdl(end) = Z(end);
AxOut = CLASS_Az0HA0.FabricOrientation_i(PAHV);
HAOut.FO = CLASS_Az0HA0.HorizontalAnisotropy_i(Z,Zmdl,AxOut.iFO,Psi); % Initial Horizontal Anisotropy
HAOut.MeanFO = CLASS_Az0HA0.HorizontalAnisotropy_i(Z,Zmdl,AxOut.iFO_Mean,Psi); % Initial Horizontal Anisotropy
%% Standard Plot
fg = [];
ax = [];
cb = [];
if FigVis ~= 0
    [fg,ax,cb] = CLASS_FixedPlot.StandardFigure([],ObsDta,ao,Z,mean(Bed),pltdim);
    fg.InvertHardcopy = 'off';
end
%% Save the Inversion Results
OutputParameters.HH = HH;
OutputParameters.VH = VH;
OutputParameters.HV = HV;
OutputParameters.VV = VV;
OutputParameters.ObsDta = ObsDta;
OutputParameters.SiteInfo= SiteInfo;
OutputParameters.f = f;
OutputParameters.ao = ao;
OutputParameters.Z = Z;
OutputParameters.dZ = dZ;
OutputParameters.Zmx = max(Z);
OutputParameters.Bed = Bed;
OutputParameters.ObsDta = ObsDta;
OutputParameters.Zmdl = Zmdl;
OutputParameters.AxOut = AxOut;
OutputParameters.HAOut = HAOut;
end