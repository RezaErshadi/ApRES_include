function [v1_0,ax] = FUNC_SpecificLayers_Ver4(AxOut,HA0,r0,Zmdl,Z,ObsDta,dZ,ax)
    dm = 50;
    MHB = (dm:dm:Z(end))'; MHB(end) = Z(end);
    
    FO = AxOut.FO;
    FO_Mean = AxOut.FO_Mean;
    
%     FO = [FO FO+180];
%     i1 = 1;
%     for i = 1:size(Zmdl,1)
%         [~,i2] = min(abs(Zmdl(i)-Z));
%         nFO(i,:) = mean(FO(i1:i2,:));
%         i1 = i2+1;
%     end
%     nFO = nFO(:,3);
%     if mean(nFO) <= 90
%         nFO = [nFO nFO+90];
%     else
%         nFO = [nFO-90 nFO];
%     end
%     nFO(nFO>180) = nFO(nFO>180)-180;
%     nFO(nFO<1) = nFO(nFO<1)+180;
%     nFO = round(nFO,0);
    
    Ax0 = repmat(FO_Mean,length(Zmdl),1);

    tic
    for j = 1:size(Ax0,2)
        OP0 = CLASS_FM.BeginForwardModel(Zmdl,HA0,r0,Ax0(:,j),dZ,0);
        EstPar = OP0.Dta;
        SignCheck_CP(:,j) = func_SignCheck(ObsDta{14},EstPar{14},MHB,Z);
    end
    toc
    [~,im] = max(mean(SignCheck_CP));
    Ax0 = [Ax0(:,im)-30 Ax0(:,im) Ax0(:,im)+30];
    txt = repmat('Ax0(1,:),',1,size(MHB,1));
    txt(end) = [];
    eval(strcat('b = combvec(',txt,');')) 
    
    Ax0Comb = nan(length(Zmdl),length(b));
    i1 = 1;
    for i = 1:size(MHB,1)
        [~,i2] = min(abs(MHB(i)-Zmdl));
        Ax0Comb(i1:i2,:) = repmat(b(i,:),(i2-i1+1),1);
        i1 = i2+1;
    end
    
    SignCheck_CP = nan(size(MHB,1),length(b));
    SignCheck_HV = nan(size(MHB,1),length(b));
    SignCheck_HH = nan(size(MHB,1),length(b));
    tic
    parfor j = 1:size(Ax0Comb,2)
        OP0 = CLASS_FM.BeginForwardModel(Zmdl,HA0,r0,Ax0Comb(:,j),dZ,0);
        EstPar = OP0.Dta;
        SignCheck_CP(:,j) = func_SignCheck(ObsDta{14},EstPar{14},MHB,Z);
        SignCheck_HV(:,j) = func_SignCheck(ObsDta{7},EstPar{7},MHB,Z);
        mf_HV = FUNC_GetTheMisfit(ObsDta,EstPar,7,1); % hv misfit
        misfit_HV(1,j) = norm(mf_HV);
    end
    toc
    mfHV = round(rescale(misfit_HV,0,100),0);
    scCP = round(rescale(mean(SignCheck_CP),0,100),0);
    scHV = round(rescale(mean(SignCheck_HV),0,100),0);
    meanAll = round(mean([mfHV;scHV;scCP]),0);
    temp = [b;mfHV;scHV;scCP;meanAll;round(rescale(meanAll,0,100),0)];
    
    [~,II] = sort(temp(end,:),'descend');
    C = temp(:,II);
    
%     mfHV = misfit_HV;
%     [mfHV,I] = sort(mfHV);
%     scCP = SignCheck_CP(:,I);
%     lowCP = mean(scCP)<90;
%     scHV = SignCheck_HV(:,I);
%     C = [b(:,I) ; mean(scCP) ; mean(scHV) ; mfHV ; mean([scCP;scHV;scHH])];
%     C(:,lowCP) = [];
%     [~,II] = sort(C(end,:),'descend');
%     CC = C(:,II);
    
    i1 = 1;
    for i = 1:size(MHB,1)
        [~,i2] = min(abs(MHB(i)-Zmdl));
        v1_0(i1:i2,:) = C(i,1);
        i1 = i2+1;
    end
end

function A = func_SignCheck(obs,mdl,MHB,Zmdl)
    obs(obs<0) = -1;
    obs(obs>=0) = 1;
    mdl(mdl<0) = -1;
    mdl(mdl>=0) = 1;  
    a = sum(obs == mdl,2);
    [n,m] = size(obs);
    acc_prc = round((a).*(100)./(m),1);   
    
    i1 = 1;
    for i = 1:size(MHB,1)
        [~,i2] = min(abs(MHB(i,1)-Zmdl));
        A(i,1) = mean(acc_prc(i1:i2));
        i1 = i2+1;
    end 
end
















