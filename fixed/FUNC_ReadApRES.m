function [hh,vh,hv,vv,Z,dZ,m_az,TimeStamp,f,Bed] = FUNC_ReadApRES(DtaDir,maxRange,BedRange)
%% READ ME
% This giant function makes the life much easier to raed ApRES files.
% The ApRES files for one site (including polarimetric direction and
% azimuthal orientation) must be named properly.
% After reading the files it will also save them as one file to load
% them much faster for the next run (with the same inputs)
% If the input parameters such as the maximum range were changed, it will
% create another fast file with the new inputs. This saves time for reading
% long time series files.

% Input:
% DtaDir -> Path of the folder with *DAT or *MAT files 
% Name of the files inside DtaDir must follow these rules:
% it should contain azimuthal angle (e.g. _0deg_ )
% it should contain polarisation orientation (e.g. _HH_ )

% maxRange: Maximum depth [m] to load (can't be more than 5000 m)
% BedRange: Optional, if you know the bed range otherwise set it as empty

% Output:
%%
FuncVersion = '_FastFileV26022021';
%%
ps = filesep; % File separator for current platform
%%
p=1; % pad factor (i.e. level of interpolation to use during fft)
winFun=@blackman; % window function handle
frange=[2e8 , 4e8]; %Normally [2e8,4e8]
% Density From field observations
rhos=407.0613; % surface density (kg/m^3)
rhoi=907.7165; % ice density (kg/m^3)
Lrho=39.5512;  % half depth decay (m) (assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho))
nI=1.68;
% bulkAlignRange=[80 100];
% dzPowerNorm=10;
% dzdiffphase=10;
% RelPowerCutOff=999999; 
% diffphaseCutOff=99999;
% DepthCutoff=-99999;
%%
% You can create a MetaData.txt for your ApRES data. 
% It should have same structure as the Template
% MetaDataPath = string(DtaDir)+ps+"MetaData.txt";
% if isfile(MetaDataPath)
%     fid = fopen(MetaDataPath,'r');
%     MD1 = fscanf(fid,'%s');
%     fclose(fid);
%     MD2 = split(MD1,';');
%     i = 1;
%     while MD2{i} ~= "end"
%         a1 = split(MD2{i},':');
%         if contains(a1{1},'(')
%             a2 = split(a1{1},'(');
%             a2 = a2{1};
%         else
%             a2 = a1{1};
%         end
%         MetaData(i,1) = string(a2);
%         MetaData(i,2) = string(a1{2});
%         i = i+1;
%     end
% else
%     MetaData = [];
% end
%% find the ApRES files in the selected folder
DataList = [dir(string(DtaDir)+ps+"*.dat") ; dir(string(DtaDir)+ps+"*.DAT") ; dir(string(DtaDir)+ps+"*.mat")];
for i = 1:length(DataList)
    DataInfo(i,1) = string(DataList(i).name);
    DataInfo(i,2) = string(DataList(i).folder);
end
%% Read/Load data
% Set the load type to "read data"
loadtype = "NormalFile";
% check if there is a fast data file
isfast = find(contains(DataInfo(:,1),strcat('p',string(p),FuncVersion,'.mat')), 1);
if ~isempty(isfast)
    % Set the load type to "fast data"
    loadtype = "FastFile";
end
switch loadtype
    case "FastFile"
        disp('Loading FastFile...')
        fastpath = strcat(DataInfo(isfast,2),ps,DataInfo(isfast,1));
        load(fastpath);
    case "NormalFile"
        disp('Loading NormalFile...')
    %--------- find the azimuth and antenna orientation from the name of the file
    % Ignore files if they have no azimuthal orientation in the file name
    DataInfo(~contains(DataInfo(:,1),'deg'),:) = [];
    % Ignore files if they have no polarization plane in the file name
    temp = [contains(DataInfo(:,1),'HH') contains(DataInfo(:,1),'VH') ... 
            contains(DataInfo(:,1),'HV') contains(DataInfo(:,1),'VV')];
    DataInfo(sum(temp,2)==0,:) = [];
    % Split the name parts based on the position of "_"
    SplitDataName = split(DataInfo(:,1),'_');
    % Sort the data based on their recorded azimuth
    az = SplitDataName(contains(SplitDataName,'deg'));
    DataInfo(:,3)  = replace(az,'deg','');
    az = str2double(DataInfo(:,3));
    [az,iaz] = sort(az);
    m_az = unique(az'); % measured azimuth
    DataInfo = DataInfo(iaz,:);
    % Sort the data with same recorded azimuth based on their antenna polarization
    DataInfo(contains(DataInfo(:,1),'HH'),4)  = "HH";
    DataInfo(contains(DataInfo(:,1),'VH'),4)  = "VH";
    DataInfo(contains(DataInfo(:,1),'HV'),4)  = "HV";
    DataInfo(contains(DataInfo(:,1),'VV'),4)  = "VV";
    m_or = unique(DataInfo(:,4)); % measured orientation
    %
    hh = []; vh = []; hv = []; vv = [];
    Z = [];
    f = nan;
    indNoDta = [];
    sz = [4 length(m_az)+1];
    varTypes = ["string",repmat(["datetime"],1,length(m_az))];
    varNames = ["Polarization",string(m_az)];
    TimeStamp = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    TimeStamp(:,1) = {'HH';'VH';'HV';'VV'};
    % Load the data (all the recorded azimuths and polarizations)
    for i = 1:length(m_az) % read existing azimuths [°]
        for j = 1:length(m_or) % read existing orientations (HH,VH,HV,VV)
            % find all the files in the selected folder with same orientation and azimuth (different time)
            iload = find(DataInfo(:,3) == string(m_az(i)) & DataInfo(:,4) == m_or(j));
            % concatenate all the recorded data with same azimuth and
            % polarization at different times
            Dta = func_LoadThemAll(iload,DataInfo,p,winFun,ps,frange);
            if isempty(Z)
                Z = Dta.range;
            end
            if isnan(f) || Dta.f~=f
                f = nanmean([f,Dta.f]);
            end
            if isnan(Dta.specCor)
                indNoDta(end+1,1) = i;
            else
                switch m_or(j)
                    case "HH"
                        hh(:,i) = Dta.specCor.';
                        TimeStamp(1,i+1) = {Dta.TimeStamp};
                    case "VH"
                        vh(:,i) = Dta.specCor.';
                        TimeStamp(2,i+1) = {Dta.TimeStamp};
                    case "HV"
                        hv(:,i) = Dta.specCor.';
                        TimeStamp(3,i+1) = {Dta.TimeStamp};
                    case "VV"
                        vv(:,i) = Dta.specCor.';
                        TimeStamp(4,i+1) = {Dta.TimeStamp};
                end
            end
        end
    end
    if length(Z) == size(hh,1)
        hh(:,indNoDta) = nan;
    end
    if length(Z) == size(vh,1)
        vh(:,indNoDta) = nan;
    end
    if length(Z) == size(hv,1)
        hv(:,indNoDta) = nan;
    end
    if length(Z) == size(vv,1)
        vv(:,indNoDta) = nan;
    end
    % Save 
    ws_save = split(DtaDir,ps);
    ws_save = strcat(DtaDir,ps,ws_save{end},'_p',string(p),'_FastFileV26022021.mat');
    disp('Saving FastFile...')
    save(ws_save,'hh','vh','hv','vv','TimeStamp','f','Z','m_az');
end
% Update: MaxDepth , Bed, Depth Correction
    Z = func_TrueDepthCorrection(nI,rhos,rhoi,Lrho,Z);
    [~,iMD] = min(abs(Z-maxRange)); % maximum depth index
    Z = Z(1:iMD);
    Z(1) = 1e-20; Z = Z'; % depth vector [m]
    dZ = mean(diff(Z)); % depth resolution [m]
    hh = hh(1:iMD,:);
    vh = vh(1:iMD,:);
    hv = hv(1:iMD,:);
    vv = vv(1:iMD,:);
    % find the bed
    Bed = [];
    if ~isempty(BedRange)
        Bed = nan(4,length(m_az));
        for i = 1:length(m_az)
            Bed(1,i) = func_Bed(Z,hh(:,i),BedRange);
            Bed(2,i) = func_Bed(Z,vh(:,i),BedRange);
            Bed(3,i) = func_Bed(Z,hv(:,i),BedRange);
            Bed(4,i) = func_Bed(Z,vv(:,i),BedRange);
        end
    end
end

function Dta = func_LoadThemAll(iload,DataInfo,p,winFun,ps,frange)
    maxRange = 5000;
    range = [];
    specCor = [];
    f = [];
    TimeStamp = NaT;
    for k=1:length(iload) % average in case of several measurements with same setup (time laps)
        filePath = strcat(DataInfo(iload(k),2),ps,DataInfo(iload(k),1));
        try
            DtaLoad = fmcw_load(filePath,1); % load the ApRES file
            DtaLoad=fmcw_cull_freq(DtaLoad,frange);
            if size(DtaLoad.vif,1) > 1
                % Takes mean of fmcw_burst
                DtaMean = fmcw_burst_mean(DtaLoad); 
            else
                DtaMean = DtaLoad;
            end
            [TempRange,~,TempSpecCor,~] = fmcw_range(DtaMean,p,maxRange,winFun);
            range(k,:) = TempRange;
            specCor(k,:) = TempSpecCor; 
            f(k) = mean(DtaMean.f); % center frequency
            TimeStamp(k) = datetime(DtaLoad.TimeStamp,'ConvertFrom','datenum');
        catch ME
            fprintf("Error in loading file: " + DataInfo(iload(k),1) +"\n")
            rethrow(ME)
        end
    end
    TimeLaps = [];
    if ~isempty(range)
        Dta.range = mean(range,1);
        Dta.specCor = mean(specCor,1);
        Dta.f = nanmean(f);
        Dta.TimeStamp = nanmean(TimeStamp);
    else
        Dta.range = nan;
        Dta.specCor = nan;
        Dta.f = nan;
        Dta.TimeStamp = nan;
    end
end

function [Bed] = func_Bed(range,specCor,BedRange)
    Bed = nan;
    if ~isempty(BedRange)
        bn = fmcw_findbed(range,specCor,BedRange,'xcor');
        Bed = range(bn);
    end
end

function [TrueDepth] = func_TrueDepthCorrection(nI,rhos,rhoi,Lrho,inp)
    TrueDepth = [];
    if ~isempty(inp)
        ShouldBeZero=@(d,dI,L,RhoSp) -dI+d+L*(nI-1)/nI*(1-RhoSp)*(exp(-d/L)-1);
        TrueDepthfun=@(RhoSp,L,dI) FUNC_bisection(@(d) ShouldBeZero(d,dI,L,RhoSp),0,max(dI)+20);
        TrueDepth = TrueDepthfun(rhos/rhoi,Lrho,inp);
    end
end

function [specCor,ax] = func_Signal2Noise(specCor,range,maxRange,BedDepth,trshld,plt)
% masking were the signal to noise ratio drops
    SignaldB = abs(specCor);
    NoiseFloor=max(SignaldB(range>maxRange-trshld));
    NoiseFloordb=20*log10(NoiseFloor);
    ax = [];
    if plt == 1
        figure; hold all;
        ax = gca;
        plot(20*log10(abs(specCor)),range,'.r');
        set(ax,'Ydir','reverse')
        set(ax,'Color',[0.6 0.6 0.6]);
        if ~isempty(BedDepth)
            plot([NoiseFloordb,NoiseFloordb],[0,BedDepth],'-b','linewidth',2);
            plot([NoiseFloordb,NoiseFloordb],[BedDepth maxRange],'--b','linewidth',2);
            plot([min(20*log10(SignaldB)),max(20*log10(SignaldB))],[BedDepth,BedDepth],'--k','linewidth',2);
        else
            plot([NoiseFloordb,NoiseFloordb],[0,maxRange],'--k','linewidth',2);
        end
        plot(20*log10(abs(specCor(abs(specCor)>NoiseFloor))),range(abs(specCor)>NoiseFloor),'.g');
        title('Noise | Signal')
        xlabel('Power [dB]')
        ylabel('Depth [m]')
        grid on
        xlim([min(20*log10(SignaldB)),max(20*log10(SignaldB))])
        ylim([range(1) range(end)])
    end
    specCor(SignaldB<NoiseFloor)=nan;
end







