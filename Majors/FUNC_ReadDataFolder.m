function [Data,z,f] = FUNC_ReadDataFolder(DataDir,FF,material)
FuncVersion = "15032022";
ps = filesep;
ws_SaveLoad = fullfile(DataDir,"_FastFileVER"+FuncVersion+".mat");
%% Check for FastFile
IsFastFile = [];
if FF ~= 0
    IsFastFile = dir(string(DataDir)+ps+"*FastFileVER"+FuncVersion+".mat");
end
%% Read/Load Data
if ~isempty(IsFastFile) % Load FastFile
    wb = waitbar(0,'Loading FastFile');
    load(fullfile(IsFastFile.folder,IsFastFile.name));
    waitbar(1,wb);
    close(wb);
else % Read data for the first time
    DataList = dir(string(DataDir)+ps+"*.dat");
    if isempty(DataList)
        DataList = dir(string(DataDir)+ps+"*.DAT");
    elseif isempty(DataList)
        DataList = dir(string(DataDir)+ps+"*.mat");
    end

    for i = 1:size(DataList,1)
        Data(i,:) = FUNC_ExtractAntennaInfo(DataList(i).name,i);
    end
    Data = sortrows(Data,[3,4,5]);
    nFiles = size(Data,1);
    z = [];
    f = [];
    wb = waitbar(0,'Reading Files');
    for i = 1:nFiles
        tic
        ix = Data{i,1};
        filePath = fullfile(DataList(ix).folder,DataList(ix).name);
        DtaMean = FUNC_SimpleRead(filePath,material);
        Data{i,6} = char(DtaMean.TIME);
        Data{i,7} = DtaMean.Signal;
        z = mean([z DtaMean.Z],2);
        f = mean([f DtaMean.fc]);
        t2(i) = toc;
        avgTime = mean(t2);
        remainedFiles = nFiles-i;
        remainedTime = remainedFiles*avgTime;
        waitbar(i/(nFiles-1.0),wb,...
            sprintf("%i files & %i seconds left (%.2f s/f)",...
            remainedFiles,round(remainedTime,0),avgTime));
    end
    close(wb);
    Data = sortrows(Data,[3,4,5,6]);
    Data(:,1) = [];
    Data = Data';
%% Save the FastFile
    wb = waitbar(0,'Saving FastFile');
    save(ws_SaveLoad,'Data','z','f');
    waitbar(1,wb);
    close(wb);
end
