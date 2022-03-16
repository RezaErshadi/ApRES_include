function fnInfo = FUNC_ScanFileName(FileName)
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

Efield1 = string(splt(logical(sum(splt == allpos,2))));
% Efield1 = string(splt{contains(splt,allpos)});
if Efield1 == ""
    Efield1 = "ND";
end

deg = string(splt{contains(splt,"deg")});
deg = str2double(erase(deg,"deg"));
if isnan(deg)
    deg = 0;
end

if contains(Efield1,HHpos) == 1
    Efield2 = "HH";
elseif contains(Efield1,VHpos)
    Efield2 = "VH";    
elseif contains(Efield1,HVpos)
    Efield2 = "HV";  
elseif contains(Efield1,VVpos)
    Efield2 = "VV";
else
    Efield2 = "ND";
end

Rpos = string(splt{contains(splt,"Rpos")});
Tpos = string(splt{contains(splt,"Tpos")});
if Rpos == "" || Tpos == ""
    TRdist = "ND";
else
    Rpos = str2double(erase(Rpos,'Rpos'));
    Tpos = str2double(erase(Tpos,'Tpos'));
    TRdist = abs(Tpos-Rpos); 
end

fnInfo = {FileName char(Efield2) char(Efield1) deg char(TRdist)};