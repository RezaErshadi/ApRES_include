function [HH,VH,HV,VV] = QuadpoleSynthesizer(Shh,Svh,Shv,Svv,ao,ac)
% Shh,Svh,Shv,Svv --> complex signal
% ao --> Azimuthal Orientation
% ac --> Azimuthal Correction (set the a0=0 to True North)
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
% if hv and VH are different, then average them
if sum(vh == hv) ~= length(vh)
    a = (real(hv) + real(vh))./2;
    b = (imag(hv) + imag(vh))./2;
    avg_HV_VH = complex(a,b);
    hv = avg_HV_VH;
    vh = avg_HV_VH;
end
%---------------------------------------------------------
ao = ao+ac;
for kk=1:length(ao)
    t= ao(kk)*pi/180;
    HH(:,kk) = (hh*cos(t)^2) + (vv*sin(t)^2) - (sin(t)*cos(t)*(hv+vh));
    VH(:,kk) = (cos(t)^2*vh) - (sin(t)^2*hv) + (sin(t)*cos(t)*(hh-vv));
    HV(:,kk) = (cos(t)^2*hv) - (sin(t)^2*vh) + (sin(t)*cos(t)*(hh-vv));
    VV(:,kk) = (vv*cos(t)^2) + (hh*sin(t)^2) + (sin(t)*cos(t)*(hv+vh));
end