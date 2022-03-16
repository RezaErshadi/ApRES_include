function Otp = FUNC_Resample2Fine(zFine,zCoarse,InpCoarse)
    Otp = nan(size(zFine,1),size(InpCoarse,2));
    i1 = 1;
    for i = 1:length(zCoarse)
        [~,i2] = min(abs(zCoarse(i)-zFine));
        Otp(i1:i2,1) = InpCoarse(i,:);
        i1 = i2+1;
    end
end