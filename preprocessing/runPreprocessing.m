function runPreprocessing(dataDir,stitchOptions,MIPoptions,zSliceOptions,...
    makeStitchedImages,makeMIPs,writeZslices)
%function to call separate preprocessing scripts given user inputs and
%metadata

% load(fullfile(dataDir,'meta.mat'),'meta');

if makeStitchedImages && (writeZslices || makeMIPs)
    error('if doing stitching, set writeZslices and makeMIPs to false')
end

if makeMIPs
    inputdir = dataDir;
    outputdir = fullfile(dataDir, 'MIP');
    if ~exist(outputdir,'dir'), mkdir(outputdir); end
    
    batchMIP(inputdir, outputdir, MIPoptions);
end

if makeStitchedImages
    stitchImageMontages(dataDir, stitchOptions);
end

if writeZslices
    writeZstacks(dataDir,zSliceOptions)
end

%handle positions per condition based on meta.nPositions set in the above
%blocks and meta.nWells determined from the condition list
load(fullfile(dataDir,'meta.mat'),'meta');
if isempty(meta.posPerCondition)
    %automatically determine the number of positions per condition
    %assuming it is the same for each condition
    meta.posPerCondition = repmat(meta.nPositions/meta.nWells,1,meta.nWells);
elseif length(meta.posPerCondition) == 1
    meta.posPerCondition = repmat(meta.posPerCondition,1,meta.nWells);
else
    %if manually set as an array of length > 1, don't modify
end

if any(floor(meta.posPerCondition) ~= meta.posPerCondition) || (sum(meta.posPerCondition) ~= meta.nPositions) || (length(meta.posPerCondition) ~= meta.nWells)
    error('posPerCondition is not set correctly')
end
save(fullfile(dataDir,'meta.mat'),'meta');


end