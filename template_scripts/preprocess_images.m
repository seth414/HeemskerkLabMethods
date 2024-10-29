clear; close all;

%% setup
%preprocess_images is a script to set the preprocessing options for
%microscope image data. This is largely to format data to load the
%appropriate images for segmentation and previews and to stitch image
%montages
%ways of preprocessing:
%   -make MIPs
%   -write Z stacks -> just DAPI by default; can specify additional
%   channels if desired e.g. if foreground segmentation is needed
%   -stitch image grids/montages
%
%we need a boolean for each of these things to determine whether to
%do them; do any of these need to be exclusize with any others?

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
% raw data location, modify if not the same as the location of this script
dataDir = scriptPath; 

makeStitchedImages = true;
makeMIPs = false;
writeZslices = false;


% -----------metadata------------------------
%nuclear marker channel (e.g., DAPI, HOECHST; indexed from 0)
nucChannel = 0;
%description of conditions in order (should be alphabetical order that the
%system finds image files but usually that is the same as the order of
%imaging)
conditions = {'LDN0','LDN100'}; %cells are treated with 100 ng/mL BMP4 and variable LDN; LDN dose given here in nM
%names of channels, e.g. channelLabel = {'DAPI','SMAD23','pAKT','SOX17'};
channelLabel = {'H2B','SMAD4'}; %H2B::RFP nuclear marker in the first channel; GFP::SMAD4 in the second channel 
%if the number of positions per condition varies between conditions, set it
%manually; e.g., posPerCondition = [4 4 4 4 5]. If left empty ([]) then
%posPerCondition divides the number of positions by the number of
%conditions and assumes it is the same for each condition
posPerCondition = [];

manualMeta = struct();
manualMeta.nucChannel = nucChannel;
manualMeta.channelLabel = channelLabel;
manualMeta.conditions = conditions;
manualMeta.posPerCondition = posPerCondition;

meta = Metadata(dataDir, manualMeta);
save(fullfile(dataDir,'meta.mat'),'meta');

% -----------stitch options------------------------
%number of colonies/stitched grids; set this manually or leave empty to 
%determine automatically; must be manually set if stitching is done at the
%same time as image acquisition
ngrids = [];

gridsize = [2 2]; %size of the stitched image grid (height by width)
montageOverlap = 18; %percent overlap between adjacent images
montageFusionGrid = false; %set to true if the matlab script multiMontageFusion was used to generate positions for the grids

%stitching method; options: weightedsample, distweight, intensityweight, randomsample
%default and generally best option is weightedsample; distweight is
%similarly good and faster but stitching grid lines may be more visible
%with this option
stitchmode = 'weightedsample';
%file type of MIPs -> if segmenting or quantifying MIPs use .tif, otherwise
%.jpg results in much smaller files for visualization & quality control
mipext = '.jpg';

stitchOptions = struct(... 
    'gridsize', gridsize,...
    'montageOverlap', montageOverlap,...
    'stitchmode', stitchmode,...
    'ext', '',... %image filename extension (.nd2, .ims, .lif, .vsi, .czi, .tif) or leave empty to automatically look for files
    'montageFusionGrid', montageFusionGrid,... %true if multiMontageFusion script was used to generate positions
    'nucChannel', nucChannel,...
    'ngrids', ngrids,...
    'tmax',[],... %set maximum number of time points to stitch, or leave empty
    'channelLabel',{channelLabel},...
    'conditions',{conditions},...
    'mipext', mipext,...
    'snakeorder','leftright');%,... %leftright or bottomtop; sometimes Fusion-generated grids snake left to right, sometimes bottom to top
    % 'orderFromFilename',false); %use if panels are broken up from a Nikon-stitched large image and have suffixes describing position on the grid



% ------------MIP options--------------------------

% channels: channels to process, counting start from 0 here
% saveidx: for channels containing nuclear marker if multiple z-slices 
% tmax: cutoff on time points to process (e.g. when it went out of focus)
MIPoptions = struct(        'channels', 0:length(channelLabel)-1,...    
                            'saveidx',  false(1,length(channelLabel)),... 
                            'tmax',     []  );
                        
% ------------z stack options--------------------------

% channels: channels to process, counting start from 0 here
% saveidx: for channels containing nuclear marker if multiple z-slices 
% tmax: cutoff on time points to process (e.g. when it went out of focus)
zSliceOptions = struct('filepattern','',... %look for files automatically if empty, set manually if preferred, e.g. if only running on *FusionStitcher.ims files
    'nucChannel', nucChannel,...
    'writeChannels',  [],... %
    'tmax',     []  );


%% do the processing
%possibly also make this block a call to a function instead
%this would make it easier to, for instance, generate and save metadata in
%a consistent way

runPreprocessing(dataDir,stitchOptions,MIPoptions,zSliceOptions,...
    makeStitchedImages,makeMIPs,writeZslices)






