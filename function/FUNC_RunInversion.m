function [InversionOutput,SiteNumber,fg,ax] = FUNC_RunInversion(InputParameters,SiteName)
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
% Synthesizing data
[HH,VH,HV,VV] = CLASS_S2P.AzimuthSynthesizer(hh,vh,hv,vv,ao,AntRot);
% Using HH, VH, HV and VV to calculate Polarimetric Prameters:
% Power Return (1:HH,2:VH,3:HV,4:VV)
% Power Anomaly (5:HH,6:VH,7:HV,8:VV)
% Phase Difference (9:HH,10:VH,11:HV,12:VV)
% Coherence Magnitude (HHVV)
% Coherence Phase (HHVV)
% Coherence Real Part & Imaginary Part (HHVV)
% Coherence Phase Derivative (HHVV)
% Coherence Scaled Phase Derivative (HHVV)
ObsDta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
PAHH = ObsDta{5};
PAHV = ObsDta{7};
CP = ObsDta{14};
Psi = ObsDta{18};

EV = [];
if ~isempty(EigenValuesFile)
    EvDir = strcat(PP,ps,ProjectName,ps,'_Info',ps,EigenValuesFile);
    EV = FUNC_ReadEigenValues(EvDir);
end
% Plot observed radar dara
[fg,ax] = CLASS_InvPlot.InversionPlot(PAHH,PAHV,CP,Psi,Z,ao,mean(Bed),FigVis);
% plot ice core observations (Eigenvalues)
if ~isempty(EV)
    %------ plot ice core Eigenvalues
    plot(ax{9},EV(:,2),EV(:,1),'xk');
    plot(ax{9},EV(:,3),EV(:,1),'dk');
    plot(ax{9},EV(:,4),EV(:,1),'*k');
    %------ plot ice core lambda2-lambda1
    plot(ax{10},EV(:,5),EV(:,1),'xk');
end
drawnow
%% Inverse Approach - First Step
%------ define model resolution and depth vector
Zmdl = (ResZmdl:ResZmdl:Z(end))'; Zmdl(end) = Z(end);
%------ define inversion resolution and depth vector
Zinv = (resInv:resInv:Z(end))'; Zinv(end) = Z(end);
%------ get the initial fabric orientation
AxOut = CLASS_Az0HA0.FabricOrientation_i(PAHV);
plot(ax{3},AxOut.FO(:,1),Z,'.g','MarkerFaceColor','g')
plot(ax{3},AxOut.FO(:,2),Z,'.g','MarkerFaceColor','g')
plot(ax{4},AxOut.FO(:,1),Z,'.g','MarkerFaceColor','g')
plot(ax{4},AxOut.FO(:,2),Z,'.g','MarkerFaceColor','g')
%------ get the initial horizontal anisotropy
HA0 = CLASS_Az0HA0.HorizontalAnisotropy_i(Z,Zmdl,AxOut.iFO,Psi); % Initial Horizontal Anisotropy
plot(ax{10},HA0,Zmdl,'.g','MarkerFaceColor','g');
drawnow 
%------ set the initial reflection ratio
r0 = zeros(size(Zmdl));
%------ get the initial Eigenvectors
switch FabricOrientationApproach
case "Manual"
    MHB(end,1) = Z(end);
    j1 = 1;
    for i = 1:size(MHB,1)
        [~,j2] = min(abs(MHB(i,1)-Zmdl));
        [~,ii] = min(abs(MHB(i,2)-AxOut.FO_Mean(1,:)));
        v1_0(j1:j2,1) =  AxOut.FO_Mean(1,ii);
        j1 = j2+1;
    end
case "SemiAuto"
    MHB(end,1) = Z(end);
%     [v1_0,ax] = FUNC_SpecificLayers_Ver2(MHB,HA0,r0,Zmdl,Z,ObsDta,dZ,ax);
%     [v1_0,ax] = FUNC_SpecificLayers_Ver3(AxOut.FO,HA0,r0,Zmdl,Z,ObsDta,dZ,ax);
%     [v1_0,ax] = FUNC_SpecificLayers_Ver4(AxOut,HA0,r0,Zmdl,Z,ObsDta,dZ,ax);
    [v1_0,ax] = FUNC_SpecificLayers_Ver5(MHB,HA0,r0,Zmdl,ObsDta,dZ,ax);
case "Auto"
    CPcpy = CP;
    cpTemp1 = CPcpy(1:end-1,:).*CPcpy(2:end,:);
    cpTemp1 = [CPcpy(1,:) ; cpTemp1];
    cpTemp2 = cpTemp1<0;
    cpTemp3 = [Z sum(cpTemp2,2)];
    cpTemp3(cpTemp3(:,1)<100 | cpTemp3(:,1)>500,:) = [];
    [~,imhb] = max(cpTemp3(:,2));
    Ztop = cpTemp3(imhb,1);
    MHB = [Ztop ; Z(end)];
    [v1_0,ax] = FUNC_SpecificLayers(MHB,HA0,r0,Zmdl,ObsDta,dZ,ax,[14]);
end
v2_0 = v1_0+90; v2_0(v2_0>180) = v2_0(v2_0>180)-180;
plot(ax{11},v1_0,Zmdl,'.g','MarkerFaceColor','g');
plot(ax{11},v2_0,Zmdl,'.g','MarkerFaceColor','g'); 
drawnow
%------ Run a forward model with initial values
ax = FUNC_PlotOptimizedModel(ax,HA0,r0,v1_0,Zmdl,Z,dZ,ao);
%--------------------------------
in_HA = HA0;
in_r = r0;
in_v1 = v1_0;
%% Inverse Approach (Invert fabric orientation and reflection ratio)
switch InvOrder
case "v1&r"
    [out_v1,out_v2,ax] = funcInv_v1...
        (v1AM,v1ST,v1NLC,in_HA,in_r,in_v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,v1CF,options1,ax);
    in_v1 = out_v1;
    [out_r,ax] = funcInv_r...
        (rAM,rST,rNLC,in_HA,in_r,in_v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,rCF,options1,ax);
case "r&v1"
    [out_r,ax] = funcInv_r...
        (rAM,rST,rNLC,in_HA,in_r,in_v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,rCF,options1,ax);
    in_r = out_r;
    [out_v1,out_v2,ax] = funcInv_v1...
        (v1AM,v1ST,v1NLC,in_HA,in_r,in_v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,v1CF,options1,ax);
end
%% EigenValues
in_HA = HA0;
in_r = out_r;
[EigVal,out_HA,out_r] = FUNC_EstimateEigenvalues(ax,in_HA,in_r,Zmdl,ResZmdl);
plot(ax{9},EigVal(:,1),Zinv,'.-r','MarkerFaceColor','r','LineWidth',2,'MarkerSize',15);
plot(ax{9},EigVal(:,2),Zinv,'.-b','MarkerFaceColor','b','LineWidth',2,'MarkerSize',15);
plot(ax{9},EigVal(:,3),Zinv,'.-g','MarkerFaceColor','g','LineWidth',2,'MarkerSize',15);
plot(ax{10},out_HA,Zinv,'.-r','LineWidth',2,'MarkerSize',15);
plot(ax{12},out_r,Zinv,'.-r','LineWidth',2,'MarkerSize',14)
%% -------------------------------- OPTIMIZED MODEL
HA = out_HA;
r = out_r;
v1 = out_v1;
v2 = out_v2;
for i=5:8
    cla(ax{i});
end
OptimizedModel = CLASS_FM.BeginForwardModel(Zmdl,HA,r,v1,dZ,0);
EstPar = OptimizedModel.Dta;
estPA_HH = EstPar{5};
estPA_HV = EstPar{7};
estCP = EstPar{14};
estdPsi = EstPar{18};
ax = CLASS_InvPlot.UpdateInversionPlot(ax,estPA_HH,estPA_HV,estCP,estdPsi,Z,[],ao,[],[]);
fg.InvertHardcopy = 'off';
%% Save the Inversion Results
InversionOutput.SiteInfo= SiteInfo;

InversionOutput.TimeStamp = TimeStamp;
InversionOutput.f = f;
InversionOutput.ao = ao;
InversionOutput.Z = Z;
InversionOutput.dZ = dZ;
InversionOutput.Zmx = max(Z);
InversionOutput.Bed = Bed;
InversionOutput.ObsDta = ObsDta;

InversionOutput.Zmdl = Zmdl;
InversionOutput.Zinv = Zinv;

InversionOutput.coreData = EV;

InversionOutput.AxOut = AxOut;

InversionOutput.v1_0 = v1_0;
InversionOutput.v2_0 = v2_0;
InversionOutput.HA0 = HA0;
InversionOutput.r0 = r0;

InversionOutput.v1 = v1;
InversionOutput.v2 = v2;
InversionOutput.HA = HA;
InversionOutput.r = r;

InversionOutput.EigVal = EigVal;
end
%%
% Invert fabric orientation
function [v1,v2,ax] = funcInv_v1...
    (AM,ST,NLC,HA,r,v1In,Zmdl,Zinv,Z,dZ,ao,ObsDta,CF,options1,ax)
    %----------                
    options1.StepTolerance = ST;
    switch AM
    case "Pointwise"
        [v1,v2,ax] = CLASS_Inversion_Pointwise.Invert_v1...
            (HA,r,v1In,Zmdl,Zinv,dZ,ObsDta,CF,options1,ax);
    case "Legendre"
        [v1,v2,ax] = CLASS_Inversion_Legendre.Invert_v1...
            (NLC,HA,r,v1In,Zmdl,Zinv,dZ,ObsDta,CF,options1,ax);
    end
    ax = FUNC_PlotOptimizedModel(ax,HA,r,v1,Zmdl,Z,dZ,ao);
    plot(ax{11},v1,Zmdl,'.-b','LineWidth',2,'MarkerSize',15);
    plot(ax{11},v2,Zmdl,'.-r','LineWidth',2,'MarkerSize',15);
    drawnow
end
% Invert reflection ratio
function [r,ax] = funcInv_r...
    (AM,ST,NLC,HA,rIn,v1,Zmdl,Zinv,Z,dZ,ao,ObsDta,CF,options1,ax)
    %----------
    options1.StepTolerance = ST;
    switch AM
    case "Pointwise"
        [r,ax] = CLASS_Inversion_Pointwise.Invert_r...
            (HA,rIn,v1,Zmdl,Zinv,dZ,ObsDta,CF,options1,ax);
    case "Legendre"
        [r,ax] = CLASS_Inversion_Legendre.Invert_r...
            (NLC,HA,rIn,v1,Zmdl,Zinv,dZ,ObsDta,CF,options1,ax);
    end
    ax = FUNC_PlotOptimizedModel(ax,HA,r,v1,Zmdl,Z,dZ,ao);
    plot(ax{12},r,Zmdl,'.g','MarkerFaceColor','g'); 
    drawnow
end
