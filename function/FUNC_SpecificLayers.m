function [v1_0,ax] = FUNC_SpecificLayers(MHB,HA0,r0,Zmdl,ObsDta,dZ,ax,CostFunc)
    v1_0 = [];
    old_i2 = 0;
    ao1 = 0:5:179; % new azimuth vector (Coarse)
    ErrorValue = nan(length(MHB),length(ao1)); % nan matrix for error values
    for i = 1:length(MHB) % for the range specified in MHB:
        [~,i2] = min(abs(MHB(i)-Zmdl));
        cz = Zmdl(1:i2); % current range
        cHA0 = HA0(1:i2); % current HA0
        cr0 = r0(1:i2); % current r0
        AxComb = ao1; 
        Ax0Comb = repmat(AxComb,length(cz),1);
        if i ~= 1
            Ax0Comb(1:length(v1_0),:) = repmat(v1_0,1,length(ao1));
        end 
        ll = length(0:dZ:cz(end));
        for jj = 1:length(ObsDta)
            ObsC{jj} = ObsDta{jj}(1:ll,:);
        end 
        pause(1);
        % run forward models with fixed HA0 and r0 and different v1_0
        parfor j = 1:size(AxComb,2)
            OP0 = CLASS_FM.BeginForwardModel(cz,cHA0,cr0,Ax0Comb(:,j),dZ,0);
            EstPar = OP0.Dta;
            misfit = [];
            for k = 1:length(CostFunc)
                mf = FUNC_GetTheMisfit(ObsC,EstPar,CostFunc(k),1); % HH misfit
                misfit(k) = norm(mf);
            end
            ErrorValue(i,j) = sum(misfit);
        end
        [~,im1] = min(ErrorValue(i,:));
        v1_0 = Ax0Comb(:,im1);
        old_i2 = [old_i2 i2];
        v1_0 = func_CheckPhaseSign(v1_0,cz,cHA0,cr0,dZ,ObsC,old_i2);
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