function FoldSumm = FUNC_FolderSummary(AntInfo)

VV = func_sortdat(AntInfo,"VV");
HH = func_sortdat(AntInfo,"HH");
VH = func_sortdat(AntInfo,"VH");
HV = func_sortdat(AntInfo,"HV");


FoldSumm = [HH VH HV VV];

function Vout = func_sortdat(AntInfo,or)
i = AntInfo(2,:) == or; % find the place of all the data with same orientation (e.g. HH)
ii = find(AntInfo(2,:) == or); % find also their true original index for later
temp1 = AntInfo(:,i); % new variable with only the selected orientations data
temp1 = [temp1 ; ii]; % record their original index for later
% sort them based on their azimuthal degree
[~,is1] = sort(str2double(temp1(4,:))); 
temp1 = temp1(:,is1);

Vout = [];
unqOr = unique(temp1(3,:));
for j = 1:length(unqOr)
    jj = temp1(3,:) == unqOr(j);
    Vout = [Vout temp1(:,jj)];
end









