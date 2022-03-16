function [DATA,f,z,twt,DATAstruct] = FUNC_readApRESfile(filePath,p,att,material,tdc)

splt = split(filePath,filesep);
fnInfo = FUNC_ScanFileName(splt{end});
splt = split(splt{end},'.');
fn = splt{1};
%%
winFun=@blackman; % window function handle
frange=[2e8 , 4e8]; %Normally [2e8,4e8]
%% Read the file
% load the raw data
DtaLoad = fmcw_load(filePath);

% spplit the data by the number of attenuator
DtaLoad = fmcw_burst_split_by_att(DtaLoad);

% selecet oen attenuator
DtaLoad = DtaLoad(att);

% make sure the frequency is in range
% DtaLoad = fmcw_cull_freq(DtaLoad,frange);

% check the sub-bursts quality and remove the bad ones
% (needs its own function)

% average all the subbursts
DtaMean = fmcw_burst_mean(DtaLoad);

% add a new variable to the data structure: when this data was recorded
DtaMean.TIME = datetime(DtaLoad.TimeStamp,'ConvertFrom','datenum');
%% Depth and Time
if material == "ice"
    maxRange = 20000;
    [Range,~,SpecCor,~] = fmcw_range(DtaMean,p,maxRange,winFun);
elseif material == "soil"
    % add soil stuff (ask Susanne)
end

%------- add a two-way-traveltime function
DtaMean.TWT = nan(length(SpecCor),1);

% true Depth Correction
TrueDepth = FUNC_TrueDepthCorrection(Range,material,tdc);
TrueDepth(1) = 1e-20;

% add a new variable to the data structure: depth
DtaMean.Z = TrueDepth';

% add a new variable to the data structure: complex signal
DtaMean.Signal = SpecCor.';
%% Extract the important variables
z = DtaMean.Z;
twt = DtaMean.TWT;
f = DtaMean.fc;
%--- in Cell format
DATA = fnInfo';
DATA{end+1,1} = char(DtaMean.TIME);
DATA{end+1,1} = DtaMean.processing;
DATA{end+1,1} = [DtaMean.Attenuator_1 DtaMean.Attenuator_2];
DATA{end+1,1} = DtaMean.Signal;

%--- in structure format
DATAstruct.filename = fn;
DATAstruct.time = char(DtaMean.TIME);
DATAstruct.fullPath = DtaMean.filename;

DATAstruct.processing = DtaMean.processing;
DATAstruct.chirpAtt = [DtaMean.Attenuator_1 DtaMean.Attenuator_2];
DATAstruct.temerature = [DtaMean.Temperature_1 DtaMean.Temperature_1];
DATAstruct.f = DtaMean.fc;

DATAstruct.constants.K = DtaMean.K;
DATAstruct.constants.T = DtaMean.T;
DATAstruct.constants.B = DtaMean.B;
DATAstruct.constants.dt = DtaMean.dt;
DATAstruct.constants.er = DtaMean.er;
DATAstruct.constants.ci = DtaMean.ci;
DATAstruct.constants.lambdac = DtaMean.lambdac;

DATAstruct.signal = DtaMean.Signal;
DATAstruct.vif = DtaMean.vif;
DATAstruct.t = DtaMean.t;



