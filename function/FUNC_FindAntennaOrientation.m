function AntInfo = FUNC_FindAntennaOrientation(FileName)
HHpos = ["HH" "hh" "h1h1" "h1h2" "h2h1" "h2h2" "><" ">>" "<<" "<>"];
VHpos = ["VH" "vh" "v1h1" "v1h2" "v2h1" "v2h2" "v>" "v<" "^>" "^<"];
HVpos = ["HV" "hv" "h1v1" "h1v2" "h2v1" "h2v2" "<v" "<^" ">v" ">^"];
VVpos = ["VV" "vv" "v1v1" "v1v2" "v2v1" "v2v2" "^^" "^v" "vv" "v^"];
%%
disp(FileName);
[~,~,ext] = fileparts(FileName);
FileName = erase(FileName,ext);
splt = split(FileName,"_");
Efield = string(splt{contains(splt,"Efield")});
Efield = erase(Efield,"Efield");
deg = string(splt{contains(splt,"deg")});
deg = erase(deg,"deg");

if contains(FileName,HHpos) == 1
    AntInfo = "HH";
elseif contains(FileName,VHpos)
    AntInfo = "VH";    
elseif contains(FileName,HVpos)
    AntInfo = "HV";  
elseif contains(FileName,VVpos)
    AntInfo = "VV";  
    
AntInfo = [AntInfo ; Efield ; deg];
end