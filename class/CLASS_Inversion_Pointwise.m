classdef CLASS_Inversion_Pointwise
    methods(Static)
%% Invert Fabric Orientation
function [v1,v2,ax] = Invert_v1...
    (HA,r,v1In,Zmdl,Zinv,dZ,ObsDta,CF,options1,ax)
    %----------
    mp0 = FUNC_Resample2Coarse(Zmdl,Zinv,v1In);
    bndr(:,1) = 0.*ones(size(mp0));
    bndr(:,2) = 180.*ones(size(mp0));
    tic
    fun = @(intls)CLASS_Inversion_Pointwise.fmincon_Pointwise_v1...
        (intls,HA,r,Zmdl,Zinv,dZ,ObsDta,CF);
    mp = fmincon(fun,mp0,[],[],[],[],bndr(:,1),bndr(:,2),[],options1);
    toc
    v1 = FUNC_Resample2Fine(Zmdl,Zinv,mp);
    v1(v1>179) = v1(v1>179)-180; v1(v1<0) = v1(v1<0)+180;
    v2 = v1+90; v2(v2>179) = v2(v2>179)-179;
end
% -------------------------------------------------------------------------
function ErrorValue = fmincon_Pointwise_v1...
    (intls,HA,r,Zmdl,Zinv,dZ,ObsDta,CF)
    %----------
    v1Est = FUNC_Resample2Fine(Zmdl,Zinv,intls);
    OP = CLASS_FM.BeginForwardModel(Zmdl,HA,r,v1Est,dZ,0);
    EstPar = OP.Dta;
    for i = 1:length(CF)
        msfTyp = 1;
        mf = FUNC_GetTheMisfit(ObsDta,EstPar,CF(i),msfTyp);
        misfit(i) = norm(mf);
    end
    ErrorValue = sum(misfit);
    fprintf("****************************************************** \n");
    for i = 1:length(Zinv)
        [~,i2] = min(abs(Zinv(i)-Zmdl));
        fprintf("Depth: %.0f ---> Estimated v1: %.2f \n",Zinv(i),v1Est(i2));
    end
    fprintf("****************************************************** \n");
end
%% Invert Reflection Ratio
function [r,ax] = Invert_r...
    (HA,rIn,v1,Zmdl,Zinv,dZ,ObsDta,CF,options1,ax)
    %----------
    mp0 = FUNC_Resample2Coarse(Zmdl,Zinv,rIn);
    bndr(:,1) = -30.*ones(size(Zinv));
    bndr(:,2) = 30.*ones(size(Zinv));
    tic
    fun = @(intls)CLASS_Inversion_Pointwise.fmincon_Pointwise_r...
        (intls,HA,v1,Zmdl,Zinv,dZ,ObsDta,CF);
    mp = fmincon(fun,mp0,[],[],[],[],bndr(:,1),bndr(:,2),[],options1);
    toc
    r = FUNC_Resample2Fine(Zmdl,Zinv,mp);
end
% -------------------------------------------------------------------------
function ErrorValue = fmincon_Pointwise_r...
    (intls,HA,v1,Zmdl,Zinv,dZ,ObsDta,CF)
    %----------
    rEst = FUNC_Resample2Fine(Zmdl,Zinv,intls);
    OP = CLASS_FM.BeginForwardModel(Zmdl,HA,rEst,v1,dZ,0);
    EstPar = OP.Dta;
    for i = 1:length(CF)
        msfTyp = 1;
        mf = FUNC_GetTheMisfit(ObsDta,EstPar,CF(i),msfTyp);
        misfit(i) = norm(mf);
    end
    ErrorValue = sum(misfit);
    fprintf("****************************************************** \n");
    for i = 1:length(Zinv)
        [~,i2] = min(abs(Zinv(i)-Zmdl));
        fprintf("Depth: %.0f -------> Estimated r: %.2f \n",Zinv(i),rEst(i2));
    end
    fprintf("****************************************************** \n");
end
%%
    end
end