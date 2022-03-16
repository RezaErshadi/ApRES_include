function TrueDepth = FUNC_TrueDepthCorrection(Range,material,tdc)
if tdc == true
    if material == "ice"
        % surface density (kg/m^3)
        rhos=407.0613; 
        % ice density (kg/m^3)
        rhoi=907.7165;
        % half depth decay (m) (assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho))
        Lrho=39.5512;  
        nI=1.68;
        n=3;
        ShouldBeZero=@(d,dI,L,RhoSp) -dI+d+L*(nI-1)/nI*(1-RhoSp)*(exp(-d/L)-1);
        TrueDepthfun=@(RhoSp,L,dI) FUNC_bisection(@(d) ShouldBeZero(d,dI,L,RhoSp),0,max(dI)+20);
        TrueDepth = TrueDepthfun(rhos/rhoi,Lrho,Range);
    elseif material == "soil"


    end
else
    TrueDepth = Range;
end