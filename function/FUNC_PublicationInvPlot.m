function [fg,ax] = FUNC_PublicationInvPlot(FilePath,FD)
load(FilePath);
fields = fieldnames(InversionOutput);
for i = 1:length(fields)
    eval(strcat(fields{i},' = InversionOutput.',fields{i},';'))
end
OptimizedModel = CLASS_FM.BeginForwardModel(Zmdl,HA,r,v1,dZ,0);

OBS_pahh = ObsDta{5};
OBS_pahv = ObsDta{7};
OBS_cp = ObsDta{14};
OBS_psi = ObsDta{18};

OPT_pahh = OptimizedModel.Dta{5};
OPT_pahv = OptimizedModel.Dta{6};
OPT_cp = OptimizedModel.Dta{14};
OPT_psi = OptimizedModel.Dta{18};

[fg,ax] = CLASS_InvPlot.InversionPlot(OBS_pahh,OBS_pahv,OBS_cp,OBS_psi,Z,ao,[],1);

ax = CLASS_InvPlot.UpdateInversionPlot(ax,OPT_pahh,OPT_pahv,OPT_cp,OPT_psi,Z,[],ao,[],[]);

p1 = plot(ax{3},AxOut.FO(:,1),Z,'.','color',[0.309, 1, 0.039],'MarkerFaceColor',[0.309, 1, 0.039]);
p2 = plot(ax{3},AxOut.FO(:,2),Z,'.','color',[0.309, 1, 0.039],'MarkerFaceColor',[0.309, 1, 0.039]);
p3 = plot(ax{4},AxOut.FO(:,1),Z,'.','color',[0.309, 1, 0.039],'MarkerFaceColor',[0.309, 1, 0.039]);
p4 = plot(ax{4},AxOut.FO(:,2),Z,'.','color',[0.309, 1, 0.039],'MarkerFaceColor',[0.309, 1, 0.039]);

p5=[]; p6=[]; p7=[]; p11=[];
if ~isempty(coreData)
    p5 = plot(ax{9},coreData(:,2),coreData(:,1),'xk');
    p6 = plot(ax{9},coreData(:,3),coreData(:,1),'dk');
    p7 = plot(ax{9},coreData(:,4),coreData(:,1),'*k');
    p11 = plot(ax{10},coreData(:,5),coreData(:,1),'xk');
end

p8=[]; p9=[]; p10=[];
if ~isempty(EigVal)
    p8 = plot(ax{9},EigVal(:,1),Zmdl,'.-r','MarkerFaceColor','r','LineWidth',2,'MarkerSize',15);
    p9 = plot(ax{9},EigVal(:,2),Zmdl,'.-b','MarkerFaceColor','b','LineWidth',2,'MarkerSize',15);
    p10 = plot(ax{9},EigVal(:,3),Zmdl,'.-g','MarkerFaceColor','g','LineWidth',2,'MarkerSize',15);
end

p12 = plot(ax{10},HA,Zmdl,'.-r','LineWidth',2,'MarkerSize',15);

p13=[];
if ~isempty(FD)
    p13 = rectangle(ax{11},'Position',[FD(1) Z(1) FD(2) Z(end)],'FaceColor',[0.984, 0.6, 0.054 0.75],'EdgeColor','k');
end

p13legend = line(NaN,NaN,'LineWidth',5,'Color',[0.984, 0.6, 0.054 1]);

p14 = plot(ax{11},v1,Zmdl,'.-b','LineWidth',2,'MarkerSize',15);
p15 = plot(ax{11},v2,Zmdl,'.-r','LineWidth',2,'MarkerSize',15);

p16 = plot(ax{12},r,Zmdl,'.-r','LineWidth',2,'MarkerSize',15);

lgnbck = [0.9 0.9 0.9];
legend(ax{9},[p5 p6 p7 p8 p9 p10],...
    {'Meas. \lambda1','Meas. \lambda2','Meas. \lambda3','Est. \lambda1','Est. \lambda2','Est. \lambda3'},...
    'FontSize',11,'FontWeight','bold',...
    'Position',[0.239726038097457 0.235385996390996 0.0736161736811196 0.0902972423267899],'Color',lgnbck)
legend(ax{10},[p11 p12],{'Measured','Estimated'},'FontSize',11,'FontWeight','bold',...
    'Position',[0.411543714183217 0.294499723655584 0.0788217959351394 0.0309672939456428],'Color',lgnbck)
legend(ax{11},[p14 p15 p13legend],{'Est. v1','Est. v2' 'Flow dir.'},'FontSize',11,'FontWeight','bold',...
    'Position',[0.666566093801804 0.280434297820817 0.0720544870049136 0.0452331260263281],'Color',lgnbck)

fg.InvertHardcopy = 'off';
end