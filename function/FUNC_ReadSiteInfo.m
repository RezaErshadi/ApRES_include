function SiteInfo = FUNC_ReadSiteInfo(InfoDir,ps,ProjectName,SiteName,TimeStamp,Bed)
SiteInfo = readtable(strcat(InfoDir,ps,ProjectName,'.csv'));
ii = SiteInfo.SiteName == SiteName;
SiteInfo = SiteInfo(ii,:);
% Update the table
if isnat(SiteInfo.TimeStamp) && ~isempty(TimeStamp)
    SiteInfo.TimeStamp = mean(table2array(TimeStamp(:,2)));
    writetable(SiteInfo,strcat(InfoDir,ps,ProjectName,'.csv'),'Delimiter',',')
end
if ~isempty(Bed)
    if isnan(SiteInfo.BedDepthRadar)
        SiteInfo.BedDepthRadar = round(mean(Bed),2);
        if ~isnan(SiteInfo.SurfaceElevationREMA)
            SiteInfo.BedElevationRadar = round(SiteInfo.SurfaceElevationREMA - SiteInfo.BedDepthRadar,2);
        end
        writetable(SiteInfo,strcat(InfoDir,ps,ProjectName,'.csv'),'Delimiter',',')
    end
end
