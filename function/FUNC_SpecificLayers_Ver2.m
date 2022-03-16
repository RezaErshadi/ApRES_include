function [v1_0,ax] = FUNC_SpecificLayers_Ver2(MHB,HA0,r0,Zmdl,Z,ObsDta,dZ,ax)
    ao1 = MHB(1,2)-45:45:MHB(1,2)+45; % new azimuth vector (Coarse)
    ao1(ao1<0) = ao1(ao1<0)+180;
    ao1(ao1>180) = ao1(ao1>180)-180;
    txt = repmat('ao1,',1,size(MHB,1));
    txt(end) = [];
    eval(strcat('b = combvec(',txt,');'))
    i1 = 1;
    Ax0Comb = nan(length(Zmdl),length(b));
    for i = 1:size(MHB,1)
        [~,i2] = min(abs(MHB(i)-Zmdl));
        Ax0Comb(i1:i2,:) = repmat(b(i,:),(i2-i1+1),1);
        i1 = i2+1;
    end
    misfit_HV = nan(1,size(Ax0Comb,2));
    SignCheck = nan(size(MHB,1)+1,size(Ax0Comb,2));
    tic
    parfor j = 1:size(Ax0Comb,2)
        OP0 = CLASS_FM.BeginForwardModel(Zmdl,HA0,r0,Ax0Comb(:,j),dZ,0);
        EstPar = OP0.Dta;
        SignCheck(:,j) = func_SignCheck(ObsDta{14},EstPar{14},MHB,Z);
        mf_HV = FUNC_GetTheMisfit(ObsDta,EstPar,7,1); % hv misfit
        misfit_HV(1,j) = norm(mf_HV);
    end
    toc
    ErrorValue = [misfit_HV ; SignCheck];  

    a = [b ; ErrorValue];
    aa = a;
    aa(:,aa(end,:)~=0) = [];
    
    [~,I] = sort(aa(size(MHB,1)+1,:));
    C = aa(:,I);
    
    i1 = 1;
    for i = 1:size(MHB,1)
        [~,i2] = min(abs(MHB(i)-Zmdl));
        v1_0(i1:i2,:) = C(i,1);
        i1 = i2+1;
    end
end

function v1_0 = func_CheckPhaseSign(v1_0,Zmdl,HA0,r0,dZ,ObsC,old_i2)
    old_v1_0 = v1_0(1:old_i2(end-1));
    New_v1_0 = v1_0(old_i2(end-1)+1:old_i2(end));
    
    p = New_v1_0;
    p90 = p+90; p90(p90>180) = p90(p90>180)-180;
    % find the correct sign from the phase
    OP0 = CLASS_FM.BeginForwardModel(Zmdl,HA0,r0,v1_0,dZ,0);
    mf1 = FUNC_GetTheMisfit(ObsC,OP0.Dta,14,1);
    
    OP0 = CLASS_FM.BeginForwardModel(Zmdl,HA0,r0,[old_v1_0;p90],dZ,0);
    mf2 = FUNC_GetTheMisfit(ObsC,OP0.Dta,14,1);
    
    misfit_p = norm(mf1);
    misfit_p90 = norm(mf2);
    
    if misfit_p90<misfit_p
        pp = p90;
    else
        pp = p;
    end
    % the flip effect
    p90f = (pp - (2*(pp-90)))+90; p90f(p90f>180) = p90f(p90f>180)-180;
    tic
    OP0 = CLASS_FM.BeginForwardModel(Zmdl,HA0,r0,[old_v1_0;pp],dZ,0);
    mf1 = FUNC_GetTheMisfit(ObsC,OP0.Dta,7,1);
    toc
    
    OP0 = CLASS_FM.BeginForwardModel(Zmdl,HA0,r0,[old_v1_0;p90f],dZ,0);
    mf2 = FUNC_GetTheMisfit(ObsC,OP0.Dta,7,1);
    
    misfit_pp = norm(mf1);
    misfit_p90f = norm(mf2);
    
    if abs(misfit_p90f-misfit_pp)<1e-1
        v1_0 = [old_v1_0;pp];
    else
        if misfit_p90f<misfit_pp
            v1_0 = [old_v1_0;p90f];
        else
            v1_0 = [old_v1_0;pp];
        end
    end 
end
function A = func_SignCheck(CP_obs,CP_mdl,MHB,Zmdl)
    CP_obs(CP_obs<0) = -1;
    CP_obs(CP_obs>=0) = 1;
    CP_obs(CP_mdl<0) = -1;
    CP_mdl(CP_mdl>=0) = 1;  
    a = sum(CP_obs == CP_mdl,2);
    [n,m] = size(CP_obs);
    tr = m/4;
    b  = a-tr;
    a(b>=0) = 1;
    a(b<0) = -1;
    
    i1 = 1;
    for i = 1:size(MHB,1)
        [~,i2] = min(abs(MHB(i,1)-Zmdl));
        A(i,1) = mean(a(i1:i2,1));
        i1 = i2+1;
    end 
    
    B = A;
    B(B<0) = nan;
    B = sum(isnan(B));
    
    A = [A;B];
    
end
















