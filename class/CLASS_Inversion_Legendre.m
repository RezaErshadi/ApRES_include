classdef CLASS_Inversion_Legendre
    methods(Static)
%% Invert Fabric Orientation
function [v1,v2,ax] = Invert_v1...
    (N,HA,r,v1In,Zmdl,Zinv,dZ,ObsDta,CF,options1,ax)
    %----------
    % Convert the intial fabric orientation to legendre coefficients with N as number
    % of legendre coefficients
    mp0 = FUNC_legendrefit(v1In, N, 'inv');
    % Inversion (Legendre Coefficients)
    tic
    fun = @(intls)CLASS_Inversion_Legendre.fmincon_Legendre_v1...
        (intls,HA,r,Zmdl,Zinv,dZ,ObsDta,CF);
    mp = fmincon(fun,mp0,[],[],[],[],[],[],[],options1);
    toc
    % convert the inverted legendre coefficients to fabric orientation
    v1 = FUNC_LegendreCoeff2Data(Zmdl,mp);
    v1(v1>179) = v1(v1>179)-180; v1(v1<0) = v1(v1<0)+180;
    v2 = v1+90; v2(v2>179) = v2(v2>179)-179;
end
% -------------------------------------------------------------------------
function ErrorValue = fmincon_Legendre_v1...
    (intls,HA,r,Zmdl,Zinv,dZ,ObsDta,CF)
    %----------
    v1Est = FUNC_LegendreCoeff2Data(Zmdl,intls);
    OP = CLASS_FM.BeginForwardModel(Zmdl,HA,r,v1Est,dZ,0);
    EstPar = OP.Dta;
    for i = 1:length(CF)
        msfTyp = 1;
        mf = FUNC_GetTheMisfit(ObsDta,EstPar,CF(i),msfTyp);
        misfit(i) = norm(mf);
    end
    ErrorValue = sum(misfit);
    if sum(v1Est>360 | v1Est<-180)
        ErrorValue = 1e20;
    end
    fprintf("****************************************************** \n");
    for i = 1:length(Zinv)
        [~,i2] = min(abs(Zinv(i)-Zmdl));
        fprintf("Depth: %.0f ---> Estimated v1: %.2f \n",Zinv(i),v1Est(i2));
    end
    fprintf("****************************************************** \n");
end
%% Invert Reflection Ratio
function [r,ax] = Invert_r...
    (N,HA,rIn,v1,Zmdl,Zinv,dZ,ObsDta,CF,options1,ax)
    %----------
    % Convert the intial reflection ratio to legendre coefficients with N as number
    % of legendre coefficients
    mp0 = FUNC_legendrefit(rIn, N, 'inv');
    % Inversion (Legendre Coefficients)
    tic
    fun = @(intls)CLASS_Inversion_Legendre.funfminconLegendre_r...
        (intls,HA,v1,Zmdl,Zinv,dZ,ObsDta,CF);
    OPO = fmincon(fun,mp0,[],[],[],[],[],[],[],options1);
    toc
    % convert the inverted legendre coefficients to reflection ratio
    r = FUNC_LegendreCoeff2Data(Zmdl,OPO);
end
% -------------------------------------------------------------------------
function ErrorValue = funfminconLegendre_r...
    (intls,HA,v1,Zmdl,Zinv,dZ,ObsDta,CF)
    %----------
    rEst = FUNC_LegendreCoeff2Data(Zmdl,intls);
    OP = CLASS_FM.BeginForwardModel(Zmdl,HA,rEst,v1,dZ,0);
    EstPar = OP.Dta;
    for i = 1:length(CF)
        msfTyp = 1;
        mf = FUNC_GetTheMisfit(ObsDta,EstPar,CF(i),msfTyp);
        misfit(i) = norm(mf);
    end
    ErrorValue = sum(misfit);
    if sum(rEst>30 | rEst<-30)
        ErrorValue = 1e20;
    end
    fprintf("****************************************************** \n");
    for i = 1:length(Zinv)
        [~,i2] = min(abs(Zinv(i)-Zmdl));
        fprintf("Depth: %.0f ---> Estimated r: %.2f \n",Zinv(i),rEst(i2));
    end
    fprintf("****************************************************** \n");
end
%% -------------------------------------------------------------------------
    end
end