function writeZstacks(dataDir,opts)
%function to write tiff stacks of z slices for segmentation
%%% TODO - handle all positions being saved in a single multi-series nd2
%%% file that is not used for stitching %%%

if isfield(opts,'writeDir')
    %can specify a separate folder to which to write the z slices if
    %desired
    writeDir = opts.writeDir;
    if ~exist(writeDir,'dir'), mkdir(writeDir); end
else
    %by default write the z stacks in the data directory
    writeDir = dataDir;
end

if isfield(opts,'filepattern') && ~isempty(opts.filepattern)
    %manually specify the image file extension
    filepattern = opts.filepattern;
else
    %look for files from a list of microscope image file formats
    filepatterns = {'*FusionStitcher*.ims','*.ims','*.nd2','*.lif','*.vsi'};
    idx = 1;
    criterion = true;
    while criterion && idx <= length(filepatterns)
        filepattern = filepatterns{idx};
        listing = dir(fullfile(dataDir,filepattern));
        if ~isempty(listing)
            criterion = false;
        end
        idx = idx + 1;
    end
    disp(strcat("writing stacks for files matching pattern ",filepattern))
end

if isfield(opts,'nucChannel')
    nucChannel = opts.nucChannel;
else
    nucChannel = 0;
end

if ~isfield(opts,'tmax') || isempty(opts.tmax)
    tmax = Inf;
else
    tmax = opts.tmax;
end

%input additional channels for which to write z slices, indexed from 0
if ~isfield(opts,'writeChannels')
    writeChannels = [];
else
    writeChannels = opts.writeChannels;
end

%always write z slices for the nuclear channel
writeChannels = union(nucChannel,writeChannels);
nChannels = numel(writeChannels);

filelist = dir(fullfile(dataDir,filepattern));
%do we also want to find image files in subfolders like for Jo's paper?
%that would be implemented as follows:
%filelist = dir(fullfile(dataDir,['**/',filepattern]));
nss = ones(1,numel(filelist));
for fi = 1:numel(filelist)
    
    %name of the z slice file should follow 
    fname = filelist(fi).name;
    [~,~,ext] = fileparts(fname);
    [barefname, pidx] = parseFilename(fname,dataDir);

    disp(['writing z slices for ' fname]);
    tic
    r = bfGetReader(fullfile(dataDir, fname));
    r.setSeries(0);
    
    nZslices = r.getSizeZ();
    nt = r.getSizeT();
    
    if strcmp(ext,'.nd2') || strcmp(ext,'.lif')
        %nd2 and lif files may have multiple fields of view saved as
        %different series of one file
        ns = r.getSeriesCount;
    else
        %ims files may have multiple series for the same field of view that
        %just consist of downsampled versions of the same data
        ns = 1;
    end
    nss(fi) = ns;
    
    %account for nd2 and lif files with multiple series here
    for si = 1:ns
        if isempty(pidx) && ns > 1
            id = si;
        else
            id = pidx;
        end
        
        r.setSeries(si-1);
        
        if isempty(id)
            outfnamePattern = [barefname, '_w%.4d_t%.4d.tif'];
        else
            outfnamePattern = [barefname, sprintf('_p%.4d',id), '_w%.4d_t%.4d.tif'];
        end
        
        for ti = 1:min(tmax,nt)
            for cii = 1:nChannels
                ci = writeChannels(cii);
                for zi = 0:nZslices-1
                    if zi == 0
                        mode = 'overwrite';
                    else
                        mode = 'append';
                    end
                    im = bfGetPlane(r, r.getIndex(zi,ci,ti-1)+1);
                    %should we check for blank/mostly blank z slices and skip?
                    writename = fullfile(writeDir,sprintf(outfnamePattern,ci,ti-1));
                    imwrite(im,writename,'WriteMode',mode)
                end
            end
        end
    end
    toc
end

load(fullfile(dataDir,'meta.mat'),'meta')
meta.fileNames = {filelist.name};
%sum of the number of series in each file should give total positions
meta.nPositions = sum(nss); 
save(fullfile(dataDir,'meta.mat'),'meta')


end