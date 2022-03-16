function [v1_0,ax] = FUNC_SpecificLayers_Ver5(MHB,HA,r,Zmdl,ObsDta,dZ,ax)
    ao1 = 0:5:179;
    ErrorValue = nan(length(MHB),length(ao1)); % nan matrix for error values
    for i = 1:length(MHB) % for the range specified in MHB:
        disp(MHB(i))
        [~,i2] = min(abs(MHB(i)-Zmdl));
        cZmdl = Zmdl(1:i2); % current range
        cHA = HA(1:i2); % current HA0
        cr = r(1:i2); % current r0
        AxComb = ao1; 
        Ax0Comb = repmat(AxComb,length(cZmdl),1);
        if i ~= 1
            Ax0Comb(1:length(v1_0),:) = repmat(v1_0,1,length(ao1));
        end 
        ll = length(0:dZ:cZmdl(end));
        for jj = 1:length(ObsDta)
            ObsC{jj} = ObsDta{jj}(1:ll,:);
        end 
        % run forward models with fixed HA0 and r0 and different v1_0
        parfor j = 1:size(AxComb,2)
            v1_test = Ax0Comb(:,j);
            OP0 = CLASS_FM.BeginForwardModel(cZmdl,cHA,cr,v1_test,dZ,0);
            EstPar = OP0.Dta;
            SignCheck(i,j) = FUNC_CoherencePhaseSignFit(ObsC{14},EstPar{14});
        end
        [~,im1] = max(SignCheck(i,:));
        v1_0 = Ax0Comb(:,im1);
    end
end
















