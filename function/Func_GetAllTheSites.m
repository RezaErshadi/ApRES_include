function AllSites = Func_GetAllTheSites(PP,ps,ProjectName)
    dd = dir(strcat(PP,ps,ProjectName));
    for i = 1:length(dd)
        AllSites(i,1) = string(dd(i).name);
    end
    AllSites(AllSites(:,1)==".",:) = [];
    AllSites(AllSites(:,1)=="..",:) = [];
    AllSites(contains(AllSites(:,1),"_"),:) = [];