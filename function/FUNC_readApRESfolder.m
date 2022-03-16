function [DATA,f,z,twt] = FUNC_readApRESfolder(fldrPath,p,att,tdc,material,lff,sff)
FuncVersion = "11102021";
ps = filesep;
ws_SaveLoad = strcat(fldrPath,ps,"FastFileVER"+FuncVersion+".mat");
%% Check for FastFile
IsFastFile = [];
if lff == true
    IsFastFile = dir(string(fldrPath)+ps+"*FastFileVER"+FuncVersion+".mat");
end
%% Read/Load Data
if ~isempty(IsFastFile) % Load FastFile
    disp(IsFastFile.name);
    wb = waitbar(0,'Loading FastFile');
    load(strcat(IsFastFile.folder,ps,IsFastFile.name));
    %---- waitbar
    waitbar(1,wb);
    close(wb);
    %---- waitbar
else % Read data for the first time
    %---- waitbar
    wb = waitbar(0,'Locating Files');
    %---- waitbar
    DataList = dir(string(fldrPath)+ps+"*.dat");
    if isempty(DataList)
        DataList = dir(string(fldrPath)+ps+"*.DAT");
    elseif isempty(DataList)
        DataList = dir(string(fldrPath)+ps+"*.mat");
    end
    
    nFiles = size(DataList,1);
    z = [];
    twt = [];
    %---- waitbar
    waitbar(0, wb,'Reading Files');
    %---- waitbar
    for i = 1:nFiles
        %---- waitbar
        tic
        %---- waitbar
        tempFilePath = strcat(fldrPath,ps,DataList(i).name);
        [DATA(:,i),f,Z,TWT,~] = FUNC_readApRESfile(tempFilePath,p,att,material,tdc);
        z = mean([z Z],2);
        twt = mean([twt TWT],2);
        fnInfo(i,:) = [i DATA(1:6,i)'];
        %---- waitbar
        t2(i) = toc;
        avgTime = mean(t2);
        remainedFiles = nFiles-i;
        remainedTime = remainedFiles*avgTime;
        waitbar(i/(nFiles-1.0),wb,...
            sprintf("%i files & %i seconds left (%.2f s/f)",...
            remainedFiles,round(remainedTime,0),avgTime));
        %---- waitbar
    end
    clear("tempInfo","tempFilePath","Z","TWT");
    fnInfo = sortrows(fnInfo,[3,4,5,6]);
    reindexing = cell2mat(fnInfo(:,1));
    DATA = DATA(:,reindexing); 
    %---- waitbar
    close(wb);
    %---- waitbar
%% Save the FastFile
    if sff == true
        wb = waitbar(0,'Saving FastFile');
        save(ws_SaveLoad,'DATA','z','twt');
        waitbar(1,wb);
        close(wb);
    end
end
