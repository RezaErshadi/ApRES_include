function AntAng = FUNC_FindAntennaAngle(FileName)
pause(1)
splt = split(FileName,'_');
AntAng = str2double(splt(4));
disp(FileName)
disp(AntAng)