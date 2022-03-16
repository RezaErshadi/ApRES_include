function AntInfo = FUNC_ExtractAntennaInfo(FileName,ix)
% This function only analyse the name of the ApRES data file
% nameing order: T & R
% the orientations are defined from the perspective of a person who is
% standing between the two antenna while T is on the left side and R in on
% the right side of that person
HHpos = ["HH" "hh" "rr" "ll" "rl" "lr"];
VHpos = ["VH" "vh" "ur" "ul" "dr" "dl"];
HVpos = ["HV" "hv" "ru" "rd" "lu" "ld"];
VVpos = ["VV" "vv" "uu" "dd" "ud" "du"];
allpos = [HHpos VHpos HVpos VVpos];
%%
disp(FileName);
[~,~,ext] = fileparts(FileName);
FileName = erase(FileName,ext);
splt = split(FileName,"_");

Efield = string(splt{contains(splt,allpos)});
for i = 1:length(allpos)
    k = contains(Efield,allpos(i));
    if k == 1
        k = i;
        break
    end
end
Efield = erase(Efield,replace(Efield,allpos(k),""));
if Efield == ""
    Efield = "N/D";
end

deg = string(splt{contains(splt,"deg")});
deg = str2double(erase(deg,"deg"));
if isnan(deg)
    deg = 0;
end

if contains(Efield,HHpos) == 1
    AntInfo = "HH";
elseif contains(Efield,VHpos)
    AntInfo = "VH";    
elseif contains(Efield,HVpos)
    AntInfo = "HV";  
elseif contains(Efield,VVpos)
    AntInfo = "VV";
else
    AntInfo = "N/D";
end

AntInfo = {ix FileName char(AntInfo) char(Efield) deg};