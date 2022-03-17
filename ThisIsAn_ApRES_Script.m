function P = ThisIsAn_ApRES_Script(fnfp)
tic
[filepath,scriptName,~] = fileparts(fnfp);
cd(filepath)
SplitPath = split(filepath,filesep); % separet the folders in the current path (splitting)
NumFolder = length(SplitPath); % number of folders in the current path
%% find the parent ApRES folder
for i = NumFolder:-1:1
    if SplitPath{i} == "ApRES"
        break
    else
        cd ..
    end
end
%% set the HomeDire as ApRES folder
HomeDir = pwd;
%% Clear the ApRES from the path
pth = path;
if ~contains(pth,"_ApRES_include")
    IncludeDir = fullfile(HomeDir,'_ApRES_include');
    addpath(genpath(IncludeDir));
end
%% Clear path
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(genpath(fullfile(HomeDir,'projects')));
%% find the project's name
prjN = SplitPath{find(SplitPath == "projects")+1};
%% Add the current project to the path
ProjectsPath = fullfile(HomeDir,'projects',prjN);
addpath(genpath(ProjectsPath));
%%
cd(HomeDir);
P.include = HomeDir;
P.Project = ProjectsPath;
P.data = fullfile(HomeDir,'projects',prjN,'data');
P.script = scriptName;
pause(0.1);
toc