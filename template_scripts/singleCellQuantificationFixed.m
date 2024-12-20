clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
baseDir = scriptPath;

%% setup

dirs = {''};

dataDirs = fullfile(baseDir,dirs);
nrounds = length(dataDirs);
r1 = 1; %first round; should be 1 in general

channelLabel = cell(1,nrounds);
nucChannels = NaN(1,nrounds);

%load metadata from each round
metas = cell(nrounds,1);
for ri = 1:nrounds
    meta = load(fullfile(dataDirs{ri},'meta.mat'),'meta');
    metas{ri} = meta.meta;
    channelLabel{ri} = metas{ri}.channelLabel;
    nucChannels(ri) = metas{ri}.nucChannel;
end
npos = metas{r1}.nPositions;
nucChannel = nucChannels(r1);

% define imageType as micropatterned 'MP' or 'disordered'
imageType = 'disordered';
% radii of micropatterns in micron for each condition (ignored if not 'MP')
radii = 700*ones(1,metas{r1}.nWells)/2; 

segDir = dataDirs{r1}; %segment images in the first directory

disp('loading nuclear and cytoplasmic masks')
tic
masks = cell(1,npos); cellData = cell(1,npos); bgmasks = cell(1,npos);
for pi = 1:npos
    if isempty(metas{r1}.fileNames)
        [~,bare,~] = fileparts(metas{r1}.filenameFormat);
        id = pi - 1;
        prefix = sprintf(bare,id,nucChannel,0);
    else
        [barefname, id] = parseFilename(metas{r1}.fileNames{pi},segDir);
        if isempty(id)
            bare = [barefname,'_w%.4d_t%.4d'];
            prefix = sprintf(bare,nucChannel,0);
        else
            bare = [barefname,'_p%.4d_w%.4d_t%.4d'];
            prefix = sprintf(bare,id,nucChannel,0);
        end
    end
    segname = [prefix,'_masks.mat'];
    mask = load(fullfile(segDir,segname));
    masks{pi} = mask.masks; 
    cellData{pi} = mask.cellData; 
    bgmasks{pi} = mask.bgmask;
end
toc

%load info for alignment between rounds of iterative imaging
fname = fullfile(baseDir,'alignmentInfo.mat');
if exist(fname,'file')
    load(fname,'shiftyxs','zshifts','zscales')
else
    shiftyxs = []; zshifts = []; zscales = [];
end

%% iterate over colonies, folders, and channels, and read out intensities

if strcmp(imageType,'MP')
    positions(metas{r1}.nPositions) = Colony();
elseif strcmp(imageType,'disordered')
    positions(metas{r1}.nPositions) = Position();
end

tic
%parpool(4)
parfor pi = 1:metas{r1}.nPositions

    condi = find(pi >= metas{r1}.conditionStartPos,1,'last'); %#ok<PFBNS> % condition index
          
    fprintf('colony %d of %d\n',pi,npos)

    if strcmp(imageType,'MP')
        positions(pi) = Colony(metas{r1}, pi);
        positions(pi).setRadius(radii(condi), metas{r1}.xres); %#ok<PFBNS>
        positions(pi).well = condi;
    elseif strcmp(imageType,'disordered')
        positions(pi) = Position(metas{r1}, pi);
    end

    % shared properties between micropattern and disordered
    positions(pi).dataChannels = 1:positions(pi).nChannels;
    positions(pi).cellData = cellData{pi};
    positions(pi).ncells = length(masks{pi});
    
    if ~isempty(metas{r1}.fileNames)
        [barefname, id] = parseFilename(metas{r1}.fileNames{pi},segDir);
        positions(pi).filenameFormat = metas{r1}.fileNames{pi};
        positions(pi).setID(id+1);
    else
        positions(pi).filenameFormat = metas{r1}.filenameFormat;
        positions(pi).setID(pi);
    end

    for ri = 1:nrounds
        fprintf('round %d of %d\n',ri,nrounds)
        nchannels = metas{ri}.nChannels;
        ncells = positions(pi).ncells;
        nucLevel = NaN(ncells, nchannels); 
        cytLevel = NaN(ncells, nchannels);
        NCratio = NaN(ncells, nchannels); 
        BG = NaN(1,nchannels);
        for ci = 1:nchannels
            img = positions(pi).loadImage(dataDirs{ri},ci-1,1); %#ok<PFBNS>
            
            if ri > 1 %for rounds 2 and later, align the images to the first round before reading out intensities
                %apply rigid shift to the image stack
                shiftyx = shiftyxs(pi,:,ri-1); %#ok<PFBNS>
                img = xyalignImageStack(img,shiftyx);
                %apply z shift + scaling
                shiftz = zshifts(pi,ri-1); %#ok<PFBNS>
                scalez = zscales(pi,ri-1); %#ok<PFBNS>
                img = zShiftScale(img,shiftz,scalez,metas{r1}.nZslices);
            end
            
            [nL,cL,ncR,bg] = readIntensityValues(img,masks{pi},bgmasks{pi});
            nucLevel(:,ci) = nL; 
            cytLevel(:,ci) = cL; 
            NCratio(:,ci) = ncR;
            BG(ci) = bg;
        end
        if ri==1
            nucLevels = nucLevel;
            cytLevels = cytLevel;
            NCratios = NCratio;
            bgs = BG;
        else
            nucLevels = cat(2, nucLevels, nucLevel); 
            cytLevels = cat(2, cytLevels, cytLevel);
            NCratios = cat(2, NCratios, NCratio);
            bgs = cat(2, bgs, BG);
        end
    end

    positions(pi).cellData.nucLevel = nucLevels;
    positions(pi).cellData.cytLevel = cytLevels;
    positions(pi).cellData.NCratio = NCratios;
    positions(pi).cellData.background = bgs; 

    if strcmp(imageType,'MP')
        % margin around nominal radius to exclude cells from, in micron
        % 96h colonies extend beyond 350 um
        margin = 50; 
        positions(pi).setCenter(margin);
    end
end
toc

meta = metas{r1};
meta.channelLabel = cat(2,channelLabel{:});

save(fullfile(baseDir,'positions.mat'),'positions')
save(fullfile(baseDir,'meta_combined.mat'),'meta')

%% load data if processing has already been done

load(fullfile(baseDir,'positions.mat'),'positions')
load(fullfile(baseDir,'meta_combined.mat'),'meta')


%% Load live data, make LineageTrace object

%specify the path to the live data folder - here I assume that the fixed
%data and live data are in the same parent directory and the name of the
%folder containing live data is 'live'
[parentDir,~,~] = fileparts(baseDir);
liveDir = fullfile(parentDir,'live');
fixedDir = baseDir;

liveMeta = load(fullfile(liveDir,'meta.mat'));
liveMeta = liveMeta.meta;
livePos = load(fullfile(liveDir,'positions.mat'));
livePos = livePos.positions;

fixedMeta = load(fullfile(fixedDir,'meta_combined.mat'));
fixedMeta = fixedMeta.meta;
fixedPos = load(fullfile(fixedDir,'positions.mat'));
fixedPos = fixedPos.positions;

channelLabels = fixedMeta.channelLabel;

npos = length(livePos);


%% initialize lineage trace & map live to fixed
lt = LineageTrace(livePos, fixedPos, liveMeta, fixedMeta, liveDir, fixedDir);

close all
positionIdx = 1:npos; %select the position indices for which to link live to fixed data
maxDist = 18; %maximum distance to link a nucleus in the live data to one in the fixed data in pixels
mapPoints(lt, maxDist, struct('alignmentMethod','automatic',...
    'positionIdx',positionIdx));

%for each position for which linking was done, check the links by eye
for ii = positionIdx
    lt.checkAlignment(ii);
    cleanSubplot(18)
    title(sprintf('Position #%d',ii))
    pause
    clf
end
close all

% load validated tracking results & assign fate
load(fullfile(liveDir,'validatedTracking.mat'),'verified')

%fields of cellData for which to build histories of live cell data
fields = {'XY','NCratio','nucLevel','cytLevel','nucArea','nucZ',...
    'nucMajorAxis','nucMinorAxis'};

hists = cell(1,npos);
for pidx = 1:npos
    %build single-cell histories from the tracking digraph for each cell in
    %the last frame of the live data
    histories = graphSignalingHistories(livePos(pidx),livePos(pidx).G,fields,1:livePos(pidx).ncells(end));
    
    %%%% OPTIONAL - only keep histories for cells tracked all the way back to the first time point %%%%
    start_times = cellfun(@(x) x(1), {histories.Time});
    hists{pidx} = histories(start_times==1);
    
    XYfinal = lt.live_position(pidx).cellData(end).XY;
    mapped = lt.mapped_idxs{pidx};
    nchan = size(lt.fixed_position(pidx).cellData.nucLevel,2);
    
    %only keep histories mapped to fixed data
    %this only works if the positions were not cropped for tracking
    cellidxs = cellfun(@(x) x(end), {hists{pidx}.CellIdxs});
    hists{pidx} = hists{pidx}(ismember(cellidxs,mapped));
    
    for ii = 1:length(hists{pidx})
        %mark each cell as verified or not
        if ~isempty(verified{pidx})
            hists{pidx}(ii).verified = verified{pidx}(hists{pidx}(ii).CellIdxs(end));
        else
            hists{pidx}(ii).verified = false;
        end
        %determine index of final cell in the history
        xy = hists{pidx}(ii).XY(end,:);
        d = sum((XYfinal - xy).^2,2);
        [~,I] = min(d);
        if I ~= hists{pidx}(ii).CellIdxs(end)
            error('indexing is messed up somewhere')
        end
        fixedIdx = find(mapped == I);
        hists{pidx}(ii).fateMarkers = lt.fixed_position(pidx).cellData.nucLevel(fixedIdx,:);
    end
    lt.histories{pidx} = hists{pidx};
end

save(fullfile(baseDir,'lt.mat'),'lt')

%% Or load existing LineageTrace object
load(fullfile(baseDir,'lt.mat'))

%% make a matrices of signaling histories and fixed-cell data

%options for constructing single-cell histories
%if interpolatemissing is set to false, the other options are ignored
opts = struct('interpolatemissing',false,... %option to interpolate values at missing times (e.g., because of gap-closing the in tracking)
    'peakremoval','none',... %option to detect and remove peaks in nuclear intensity or nuclear-cytoplasmic ratio because of cell division
    'domedfilt',false); %option to median filter signaling histories

Mats = lt.histories2mats(1:npos,opts);

%matrix of single-cell histories of nuclear-cytoplasmic SMAD4 readout
SMAD4 = Mats.NCratio(:,:,2);
%vector of pSMAD1 readout values
pSMAD1 = Mats.fateMarkers(3,:);
%boolean vector indicating whether each signaling history is manually
%validated/corrected
verified = Mats.verified;













