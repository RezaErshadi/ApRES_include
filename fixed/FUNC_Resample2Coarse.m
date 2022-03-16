function Otp = FUNC_Resample2Coarse(zFine,zCoarse,InpFine)
    Otp = nan(size(zCoarse,1),size(InpFine,2));
    i1 = 1;
    for i = 1:size(zCoarse,1)
        [~,i2] = min(abs(zCoarse(i,1)-zFine));
        Otp(i,1) = mean(InpFine(i1:i2));
        i1 = i2+1;
    end      
end