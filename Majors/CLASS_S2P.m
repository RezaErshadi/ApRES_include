classdef CLASS_S2P
    methods(Static)
%% Synthesizing data for any azimuthal angle
    function [HH,VH,HV,VV] = AzimuthSynthesizer(Shh,Svh,Shv,Svv,theta_prime,Rotation)
        %---------------------------------------------------------
        hh = Shh(:,1);
        vh = Svh(:,1);
        hv = Shv(:,1);
        vv = Svv(:,1);
        %---------------------------------------------------------
        HH = zeros(size(Shh));
        VV = zeros(size(Svv));
        HV = zeros(size(Shv));
        VH = zeros(size(Svh));
        %---------------------------------------------------------
%         a = (real(hv) + real(vh))./2;
%         b = (imag(hv) + imag(vh))./2;
%         avg_HV_VH = complex(a,b);
%         hv = avg_HV_VH;
%         vh = avg_HV_VH;
        vh = hv;
        %---------------------------------------------------------
        theta_prime = theta_prime+Rotation;
        for kk=1:length(theta_prime)
            t= theta_prime(kk)*pi/180;
            HH(:,kk) = (hh*cos(t)^2) + (vv*sin(t)^2) - (sin(t)*cos(t)*(hv+vh));
            VH(:,kk) = (cos(t)^2*vh) - (sin(t)^2*hv) + (sin(t)*cos(t)*(hh-vv));
            HV(:,kk) = (cos(t)^2*hv) - (sin(t)^2*vh) + (sin(t)*cos(t)*(hh-vv));
            VV(:,kk) = (vv*cos(t)^2) + (hh*sin(t)^2) + (sin(t)*cos(t)*(hv+vh));
        end
    end  
%% Extracting parameters from the synthesized signal
    function Dta = Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,dtatyp)                            
        dZ = mean(diff(Z));
        %--------------------------------------------------------- Power anomaly & lateral phase difference
        [PA_HH,PHI_HH] = CLASS_S2P.AmpPhs(HH,ao);
        [PA_VH,PHI_VH] = CLASS_S2P.AmpPhs(VH,ao);
        [PA_HV,PHI_HV] = CLASS_S2P.AmpPhs(HV,ao);
        [PA_VV,PHI_VV] = CLASS_S2P.AmpPhs(VV,ao);
        %--------------------------------------------------------- Polarimetric Coherence
        C_HHVV = CLASS_S2P.Chhvv(HH,VV,C_DepthWin,dZ,dtatyp);
        absC = abs(C_HHVV);
        argC = angle(C_HHVV);
        reC = real(C_HHVV);
        imC = imag(C_HHVV);
        if isempty(C_ConvWin)
            C_ConvWin = C_DepthWin;
        end
        %--------------------------------------------------------- Polarimetric Coherence depth gradient
        gradC = CLASS_S2P.Chhvv_Z_derivative(C_HHVV,C_ConvWin,dZ);
        Psi = (2*299792458*sqrt(3.15)*gradC) ./ (4*pi*f*0.034);
        %--------------------------------------------------------- Output structure
        % PA,PD : HH VH HV VV
        % Coh: abs, arg, re, im , grad, Psi
        Sig.HH = HH; Sig.VH = VH; Sig.HV = HV; Sig.VV = VV;
        PA.HH = PA_HH; PA.VH = PA_VH; PA.HV = PA_HV; PA.VV = PA_VV;
        PD.HH = PHI_HH; PD.VH = PHI_VH; PD.HV = PHI_HV; PD.VV = PHI_VV;
        C.absC = absC; C.argC = argC; C.reC = reC; C.imC = imC;
        C.gradC = gradC; C.Psi = Psi;
        %--------------------------------------------------------- Denoising
        [PA,PD] = CLASS_Denoising.RunDenoising(DenoisingFlag,PA,PD,dZ);
        PA_HH = PA.HH; PA_VH = PA.VH; PA_HV = PA.HV; PA_VV = PA.VV;
        PHI_HH = PD.HH; PHI_VH = PD.VH; PHI_HV = PD.HV; PHI_VV = PD.VV;
        %--------------------------------------------------------- Store
        Dta = {  HH,VH,HV,VV,... % 1 2 3 4
                 PA_HH,PA_VH,PA_HV,PA_VV,... % 5 6 7 8
                 PHI_HH,PHI_VH,PHI_HV,PHI_VV,... % 9 10 11 12
                 absC,argC,reC,imC,gradC,Psi}; % 13 14 15 16 17 18
    end   
%% Power Anomaly & Lateral Phase Difference
    function [delta_P,dphi_dtheta] = AmpPhs(s,ao)
        delta_P = 20.*log10(abs(s)./nanmean(abs([s]),2)); % Power anomaly
        delta_P(isinf(delta_P)) = nan;
        s = s.';
        PhaseDiff = nan(size(s,1)-1,size(s,2));
        for i=1:length(ao)-1
             PhaseDiff(i,:)=angle(s(i+1,:).*conj(s(i,:)));
        end
        xdiffphase=cos(PhaseDiff');
        ydiffphase=sin(PhaseDiff'); 
        dphi_dtheta=atan2(ydiffphase,xdiffphase);
    end
%%
    function Chhvv = Chhvv(Shh,Svv,DptAvgC,dZ,dtatype)
        stwnd = DptAvgC/dZ;
        Chhvv1 = movsum(Shh.*conj(Svv),[stwnd],1,'omitnan');
        Chhvv2 = sqrt(movsum(abs(Shh).^2,[stwnd],1,'omitnan'));
        Chhvv3 = sqrt(movsum(abs(Svv).^2,[stwnd],1,'omitnan'));
        Chhvv = Chhvv1 ./ (Chhvv2 .* Chhvv3);
        if dtatype == "radar"
            Chhvv = conj(Chhvv); % assumption from Jordan (converting the synthetic received seignal to de-ramped signal)
        end
    end 
%%      
    function dphi_dz = Chhvv_Z_derivative(Chhvv,ConvWin,dZ)
        ConvWin =round(ConvWin/dZ,0);
        win=gausswin(ConvWin);
        win=win/sum(win);
        R = real(Chhvv);
        I = imag(Chhvv); 
%         R = CLASS_Denoising.DenoisePCA(real(Chhvv),1);
%         I = CLASS_Denoising.DenoisePCA(imag(Chhvv),1);
        Rz=[];
        Iz=[];
        Rlow=[];
        Ilow=[];
        for j=1:size(Chhvv,2)
            Rz(:,j)=diff(conv(R(:,j),win,'same'))/dZ;
            Iz(:,j)=diff(conv(I(:,j),win,'same'))/dZ;
            Rlow(:,j)=conv(R(:,j),win,'same');
            Ilow(:,j)=conv(I(:,j),win,'same');
        end
        Iz(size(I,1),:)=NaN;
        Rz(size(R,1),:)=NaN;
        dphi_dz =((R.*Iz-I.*Rz)./(Rlow.^2+Ilow.^2));
    end
    end
end