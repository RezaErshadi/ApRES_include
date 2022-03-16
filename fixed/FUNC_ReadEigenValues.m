function EV = FUNC_ReadEigenValues(FP)
    fid = fopen(FP,'r');
    Dta = textscan(fid, '%f %f %f %f %s', 'headerlines', 1, 'Delimiter',',');
    fclose(fid);
    ZEV = Dta{:,1}; %depth
    E1 = Dta{:,2}; %E1
    E2 = Dta{:,3}; %E2
    E3 = Dta{:,4}; %E2
    GirdleStrength = E2-E1;
    PoleStrength = E3-E2;
    EV = [ZEV,E1,E2,E3,GirdleStrength,PoleStrength];
end