clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
dataDir = scriptPath;

%% setup

load(fullfile(dataDir,'meta.mat'),'meta')
channelLabel = meta.channelLabel;
nucChannel = meta.nucChannel;
npos = meta.nPositions;
ntime = meta.nTime;
nchannels = meta.nChannels;

% define imageType as micropatterned 'MP' or 'disordered'
imageType = 'disordered';
% radii of micropatterns in micron for each condition (ignored if not 'MP')
radii = 700*ones(1,meta.nWells)/2; 

segDir = dataDir;

%% iterate over positions, time points, and channels, and read out intensities

if strcmp(imageType,'MP')
    positions(meta.nPositions) = Colony();
elseif strcmp(imageType,'disordered')
    positions(meta.nPositions) = Position();
end

tic
for pi = 1:meta.nPositions

    condi = find(pi >= meta.conditionStartPos,1,'last'); % condition index
          
    fprintf('position %d of %d\n',pi,npos)

    if strcmp(imageType,'MP')
        positions(pi) = Colony(meta, pi);
        positions(pi).setRadius(radii(condi), meta.xres);
        positions(pi).well = condi;
    elseif strcmp(imageType,'disordered')
        positions(pi) = Position(meta, pi);
    end
    
    if ~isempty(meta.fileNames)
        [barefname, id] = parseFilename(meta.fileNames{pi},segDir);
        positions(pi).filenameFormat = meta.fileNames{pi};
        if isempty(id)
            bare = [barefname,'_w%.4d_t%.4d'];
        else
            bare = [barefname,sprintf('_p%.4d',id),'_w%.4d_t%.4d'];
        end
    else
        [~,bare,~] = fileparts(meta.filenameFormat);
        id = pi - 1;
        bare = strrep(bare,'_p%.4d',sprintf('_p%.4d',id));
        positions(pi).filenameFormat = meta.filenameFormat;
    end
    positions(pi).setID(id+1);
    
    % shared properties between micropattern and disordered
    positions(pi).dataChannels = 1:positions(pi).nChannels;
    clear CData
    Cdata(ntime) = struct; %#ok<SAGROW>
    for ti = 1:ntime
        prefix = sprintf(bare,nucChannel,ti-1);
        segname = [prefix,'_masks.mat'];
        load(fullfile(segDir,segname),'masks','cellData','bgmask');
        fields = fieldnames(cellData);
        for fi = 1:length(fields)
            CData(ti).(fields{fi}) = cellData.(fields{fi}); %#ok<SAGROW>
        end
        ncells = length(masks);
        positions(pi).ncells(ti) = ncells;
        CData(ti).nucLevel = NaN(ncells, nchannels); 
        CData(ti).cytLevel = NaN(ncells, nchannels);
        CData(ti).NCratio = NaN(ncells, nchannels); 
        CData(ti).background = NaN(1,nchannels);
        for ci = 1:nchannels
            img = positions(pi).loadImage(dataDir,ci-1,ti);
            
            [nL,cL,ncR,bg] = readIntensityValues(img,masks,bgmask);
            CData(ti).nucLevel(:,ci) = nL; 
            CData(ti).cytLevel(:,ci) = cL; 
            CData(ti).NCratio(:,ci) = ncR;
            CData(ti).background(ci) = bg;
        end
        fprintf('.')
        if mod(ti,50) == 0
            fprintf('\n')
        end
    end
    fprintf('\n')
    positions(pi).cellData = CData;

    if strcmp(imageType,'MP')
        % margin around nominal radius to exclude cells from, in micron
        % 96h colonies extend beyond 350 um
        margin = 50; 
        positions(pi).setCenter(margin);
    end
end
toc

save(fullfile(dataDir,'positions.mat'),'positions')

%save csv files for each position to import to trackmate
for pi = 1:npos
    [dataTable, writename] = cellDataToTable(positions(pi),meta);
    writetable(dataTable,fullfile(dataDir,writename));
end

%% load data if processing has already been done
load(fullfile(dataDir,'positions.mat'),'positions')


%% load trackmate results, parse into a digraph, and add to the Position object

verified = cell(1,npos);
for pi = 1:npos
    trackmatename = fullfile(dataDir,sprintf('trackmate_p%.4d_edges.csv',pi-1));
    if exist(trackmatename,'file')
        fprintf('position %d\n',pi)
        T = readtable(trackmatename);
        G = parseTrackmateEdges(positions(pi),T);
        positions(pi).G = G;
        verified{pi} = false(positions(pi).ncells(ntime),1);
    end
end
save(fullfile(dataDir,'positions.mat'),'positions')

save(fullfile(dataDir,'positions.mat'),'positions')
save(fullfile(dataDir,'validatedTracking.mat'),'verified')

%% or do single-cell tracking in MATLAB & save results

%user parameters for single-cell tracking
max_linking_distance = 22; %maximum movement distance between adjacent frames in pixels
max_gap_distance = 35; %maximum movement distance for gap-closing over more than 1 frame in pixels
max_gap_time = 5; %maximum number of frames over which to do gap closing

trackopts = struct(...
    'max_linking_distance', max_linking_distance,...
    'max_gap_distance',     max_gap_distance,...
    'max_gap_time',         max_gap_time,...
    'validate_links',       false,...
    'dataDir',              dataDir,...
    'dejitter',             false,...
    'imsize',               [meta.ySize, meta.xSize]);

%do tracking, store links between cells in a directed graph (digraph)
verified = cell(1,npos);
for pi = 1:npos
    fprintf('position %d of %d\n',pi,npos)
    G = create_tracks(positions(pi), trackopts);
    positions(pi).G = G;
    verified{pi} = false(positions(pi).ncells(ntime),1);
end

save(fullfile(dataDir,'positions.mat'),'positions')
save(fullfile(dataDir,'validatedTracking.mat'),'verified')












