clear; close all; clc

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
%all data folders are in subfolders of the directory containing the script
baseDir = fullfile(scriptPath);

%% setup

%list folders containing image data
dirs = {...
    'round1',...
    'round2',...
    };
dataDirs = fullfile(baseDir,dirs);
nrounds = length(dataDirs);

%load metadata from each round
channelLabel = cell(1,nrounds);
nucChannels = NaN(1,nrounds);
metas = cell(nrounds,1);
for ri = 1:nrounds
    meta = load(fullfile(dataDirs{ri},'meta.mat'),'meta');
    metas{ri} = meta.meta;
    channelLabel{ri} = metas{ri}.channelLabel;
    nucChannels(ri) = metas{ri}.nucChannel;
end
npos = metas{1}.nPositions;

%% do alignment
%directories for which to do alignment
%each new folder can be aligned as it is acquired instead of all at once
rounds = 2:nrounds;

shiftyxs = NaN(npos,2,nrounds-1);
zshifts = NaN(npos,nrounds-1);
zscales = NaN(npos,nrounds-1);

%iterate over positions
for pidx = 1:npos
    fprintf('position %d of %d\n',pidx,npos)
    %load DAPI image for the first channel
    if ~isempty(metas{1}.fileNames)
        [barefname, id] = parseFilename(metas{1}.fileNames{pidx},dataDirs{1});
        filenameFormat = metas{1}.fileNames{pidx};
    else
        filenameFormat = metas{1}.filenameFormat;
    end
    DAPI1 = loadImages(filenameFormat, dataDirs{1}, id, nucChannels(1), 1);
    
    nz1 = size(DAPI1,3);
    mip1 = max(DAPI1,[],3);
    lim1 = stitchedlim(mip1);
    
    for ri = rounds
        fprintf('round %d of %d\n',ri,nrounds)
        QCdir = fullfile(dataDirs{ri},'overlays');
        if ~exist(QCdir,'dir'), mkdir(QCdir); end
        
        if ~isempty(metas{ri}.fileNames)
            [barefname, id] = parseFilename(metas{ri}.fileNames{pidx},dataDirs{ri});
            filenameFormat = metas{ri}.fileNames{pidx};
        else
            filenameFormat = metas{ri}.filenameFormat;
        end
        DAPI2 = loadImages(filenameFormat, dataDirs{ri}, id, nucChannels(ri), 1);
        
        mip2 = max(DAPI2,[],3);
        lim2 = stitchedlim(mip2);
        
        prefix = sprintf([barefname,'_p%.4d'],pidx-1);

        %xy alignment
        shiftyx = findImageShift(mip1,mip2,'automatic');
        shiftyxs(pidx,:,ri-1) = shiftyx;
        mipa = alignImage(mip1,mip2,shiftyx);

        %save mip overlay
        overlay = makeMIPOverlay(mip1,mipa,lim1,lim2,0.075);
        savename = fullfile(QCdir,[prefix,'_MIPoverlay.png']);
        imwrite(overlay,savename);
        
        %align the entire DAPI stack in xy
        dapia = xyalignImageStack(DAPI2,shiftyx);

        %z alignment
        [shiftz, scalez] = zAlignImageStacks(DAPI1,dapia);
        zshifts(pidx,ri-1) = shiftz;
        zscales(pidx,ri-1) = scalez;
        %apply the z shift and scaling
        dapia = zShiftScale(dapia,shiftz,scalez,nz1);
        
        %save cross section overlay
        index = round(size(dapia,1)/2);
        crossSectionOverlay(DAPI1,dapia,metas{1},index,lim1,lim2,0.3)
        savename = fullfile(QCdir,[prefix,sprintf('_crossSectionX%d_overlay.png',index)]);
        saveas(gcf,savename)
%         exportgraphics(gcf, savename,'Resolution',300)
        
    end
    close all
end

save(fullfile(baseDir,'alignmentInfo.mat'),'shiftyxs','zshifts','zscales')

%% local functions

function IMa = xyalignImageStack(img,shiftyx)
m = size(img,1); n = size(img,2); nz = size(img,3); imclass = class(img);

IMa = zeros(m,n,nz,imclass);
for zi = 1:nz
    im = img(:,:,zi);
    ima = alignImage(im,im,shiftyx);    
    IMa(:,:,zi) = ima;
end

end

function imz = zShiftScale(img,shift,scale,nz1)
m = size(img,1); n = size(img,2); nz = size(img,3);
nznew = round(nz*scale);

[XX,YY,ZZ] = meshgrid(1:n,1:m,1:nz);
[Xnew,Ynew,Znew] = meshgrid(1:n,1:m,linspace(1 - shift,nz - shift,nznew));

imz = uint16(interp3(XX,YY,ZZ,single(img),Xnew,Ynew,Znew));

nznew = size(imz,3);
if nznew > nz1
    imz = imz(:,:,1:nz1);
elseif nznew < nz1
    imz = cat(3,imz,zeros(m,n,nz1-nznew,'uint16'));
end

end

function RGB = makeMIPOverlay(mip1,mip2,lim1,lim2,cfs)
m = size(mip1,1); n = size(mip1,2);

A = imadjust(mip1,lim1); B = imadjust(mip2,lim2);
overlay = cat(3,A,B,A);

%use arial font, but different matlab versions have different names for
%available fonts
fonts = listTrueTypeFonts;
font = 'Arial';
if sum(strcmp(font,fonts)) == 0
    font = 'Arial Unicode';
    if sum(strcmp(font,fonts)) == 0
        I = find(contains(fonts,'Arial'),1);
        font = fonts{I};
    end
end

%add label text
text = {'round 1', 'round 2'};
ypos = round([m*(1 - 0.75*cfs); m*(1 - 0.75*cfs)]);
xpos = round([0.005*n; 0.15*n]);
RGB = insertText(overlay,[xpos,ypos],text,'BoxColor','black',...
    'TextColor',{'m','g'},'BoxOpacity',0.3,'FontSize',96,...
    'Font',font);

end

function crossSectionOverlay(img1,img2,meta1,index,lim1,lim2,cfs)

xres = meta1.xres;
zres = meta1.zres;

cross1 = transpose(squeeze(img1(index,:,end:-1:1)));
zsize1 = round(size(cross1,1)*zres/xres);
cross1 = imresize(cross1,[zsize1,size(cross1,2)]);

cross2 = transpose(squeeze(img2(index,:,end:-1:1)));
zsize2 = round(size(cross2,1)*zres/xres);
cross2 = imresize(cross2,[zsize2,size(cross2,2)]);

if zsize1 > zsize2
    cross2 = [zeros(zsize1-zsize2,size(cross2,2),'uint16'); cross2];
elseif zsize2 > zsize1
    cross1 = [zeros(zsize2-zsize1,size(cross1,2),'uint16'); cross1];
end
zsize = size(cross1,1);
n = size(cross1,2);

C1 = imadjust(cross1,lim1); C2 = imadjust(cross2,lim2);
o1 = cat(3,C1,C2,C1);

titles = {'round 1','round 2',...
    "\color{magenta}round 1 \color{green} round 2"};
crosses = {C1,C2,o1};

figpos = figurePosition([n,2*zsize*length(crosses)]);
figure('Position',figpos)
for ii = 1:length(crosses)
    subplot_tight(length(crosses),1,ii)
    imshow(crosses{ii})
    cleanSubplot
    
    ypos = zsize*0.5*cfs; xpos = 0.01*n*cfs;
    text(xpos, ypos, titles{ii},...
        'Color','w','FontUnits','normalized','FontSize',cfs,...
        'FontWeight','bold')
end


end











