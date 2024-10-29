clear; close all; clc

%% load positions
scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
% imdir = fullfile(scriptPath,'imgs');
% imdir = 'D:\221119_Smad4GFP_BMPdoses_tracking';

dataDir = scriptPath;
load(fullfile(dataDir,'positions.mat'),'positions')
load(fullfile(dataDir,'meta.mat'),'meta')

ntime = positions(1).nTime;
nucChannel = meta.nucChannel;
npos = length(positions);
imsize = [meta.ySize, meta.xSize];

load(fullfile(dataDir,'validatedTracking.mat'),'verified')

%% choose a position index for which to start/continue validation
pidx = 1;

if isempty(verified{pidx})
    verified{pidx} = false(positions(pidx).ncells(ntime),1);
end

allpoints = cell2mat({positions(pidx).cellData.XY}');
G = positions(pidx).G;
%get X,Y information for tracks for visualization
disp('Building tracks')
tracks = graphtracks(G, ntime);
[X,Y] = trackXY(tracks, allpoints);

%% load time lapse of adjusted MIPs
%TODO: load MIPs directly instead of loading the stack and calculating the MIP

tic
%find contrast limits based on the image in the last time point
img = max(positions(pidx).loadImage(dataDir,nucChannel,ntime),[],3);
lim = stitchedlim(img);

%load nuclear image in each time point, adjust contrast, and downsample to
%uint8 to preserve memory
%make sure contrast adjustment is done before reducing the bit depth
imgs = zeros(imsize(1), imsize(2), ntime, 'uint8');
for ti = 1:ntime
    img = max(positions(pidx).loadImage(dataDir,nucChannel,ti),[],3);
    imgs(:,:,ti) = im2uint8(imadjust(img,lim));
    fprintf('.')
    if mod(ti,50) == 0
        fprintf('\n')
    end
end
fprintf('\n')
toc

%% Manual correction
%backtrack and correct any errors, marking cell done at the end and save
%tracking digraph for each position
tail_length = 12; %adjust the length of tracking tails
fs = 12; %change the font size for text overlaid on the image

ntracks = length(tracks);
%used for converting from [time, cell] indices to graph node indices
cumtimes = [0, cumsum(positions(pidx).ncells)];

%define colors for each track
colors = make_colors(ntracks);
% labelColors = distinguishable_colors(5,{'w','k'});
labelColor = [1 0 0];

close all
f = figure('WindowState','maximized');
%initialize a plot object for image display
p1 = imshow(imgs(:,:,ntime),'InitialMagnification','fit');
hold on
%initialize a scatterplot of cell centroid positions
xys = positions(pidx).cellData(ntime).XY;
f1 = scatter(xys(:,1),xys(:,2),25,'g','filled');
%initialize an array of Line objects to display tracking tails
p = initialize_lines(tracks, colors);

drawnow
tidx = ntime;
cleanSubplot
breakvar = false;
updateflag = false;
while ~breakvar
    tidxprev = tidx;
    %update image for new frame
    set(p1,'CData',imgs(:,:,tidx));
    %update scatter plot of cell centroid positions on the image
    xys = positions(pidx).cellData(tidx).XY;
    set(f1,'XData',xys(:,1),'YData',xys(:,2));
    
    idx = 1;
    %draw tracking tails
    for ii = 1:ntracks
        for k = 1:size(tracks{ii},2)
            sidx = max(1, tidx - tail_length);
            x = X{ii}(sidx:tidx,k);
            x = x(~isnan(x));
            y = Y{ii}(sidx:tidx,k);
            y = y(~isnan(y));
            if ~isempty(x)
                set(p(idx),'XData',x,'YData',y)
            else
                set(p(idx),'XData',1,'Ydata',1)
            end
            idx = idx + 1;
        end
    end
    
    if tidx == ntime
        %show indices of *all* cells colored based on whether they have
        %been verified multiply the RGB color triplet by 0.75 so it appears
        %slightly grayed out
        tt = cell(1,2);
        inds = find(verified{pidx}); ninds = find(~verified{pidx});
        tt{1} = text(xys(inds,1)+5,xys(inds,2)+5,...
            cellstr(num2str(inds(:))),'FontSize',fs,'FontWeight',...
            'bold','Color',0.75*labelColor);
        tt{2} = text(xys(ninds,1)+5,xys(ninds,2)+5,...
            cellstr(num2str(ninds(:))),'FontSize',fs,'FontWeight',...
            'bold','Color',labelColor);
    end
    
    title(sprintf('%d of %d', tidx, ntime))
    
    waitforbuttonpress;
    key = f.CurrentCharacter;
    
    if strcmp(key,'a')
        %go back one time point
        tidx = tidx - 1;
    elseif strcmp(key,'d')
        %advance one time point
        tidx = tidx + 1;
    elseif strcmp(key,'w')
        %advance ten time points
        tidx = tidx + 10;
    elseif strcmp(key,'s')
        %go back ten time points
        tidx = tidx - 10;
    elseif strcmp(key,'+')
        %zoom in 2X
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        set(gca,'xlim', mean(xl) + 0.25*(xl(2) - xl(1))*[-1 1])
        set(gca,'ylim', mean(yl) + 0.25*(yl(2) - yl(1))*[-1 1])
    elseif strcmp(key,'-')
        %zoom out 2X
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        set(gca,'xlim', mean(xl) + (xl(2) - xl(1))*[-1 1])
        set(gca,'ylim', mean(yl) + (yl(2) - yl(1))*[-1 1])
    elseif strcmp(key,'2')
        %pan down
        yl = get(gca,'ylim');
        set(gca,'ylim', yl + 0.25*(yl(2) - yl(1)))
    elseif strcmp(key,'8')
        %pan up
        yl = get(gca,'ylim');
        set(gca,'ylim', yl - 0.25*(yl(2) - yl(1)))
    elseif strcmp(key,'4')
        %pan left
        xl = get(gca,'xlim');
        set(gca,'xlim', xl - 0.25*(xl(2) - xl(1)))
    elseif strcmp(key,'6')
        %pan right
        xl = get(gca,'xlim');
        set(gca,'xlim', xl + 0.25*(xl(2) - xl(1)))
    elseif strcmp(key,'5')
        %reset field of view to default (full image)
        set(gca,'xlim', [0,size(imgs,1)] + [-0.5 0.5])
        set(gca,'ylim', [0,size(imgs,2)] + [-0.5 0.5])
    elseif strcmp(key,'e')
        %exit the while loop
        breakvar = true;
    elseif strcmp(key,'m')
        clear s t
        %make a change to the tracking results; add or remove a link
        [lindx, ~] = listdlg('ListString',...
            {'Add link','Remove link before','Remove link after','Cancel'});
        if lindx == 1 %add link
            %keep track of source node to add 
            s = pickpoint(xys,tidx,cumtimes,'source'); %source node
        elseif lindx == 2 %remove link before (into) selected node
            t = pickpoint(xys,tidx,cumtimes,'target'); %target node
            if indegree(G,t) == 1
                s = predecessors(G,t); %source node
                G = rmedge(G,s,t); %remove edge between these nodes
                updateflag = true;
            else
                disp('There are no incoming edges to this node')
            end
            clear s t
        elseif lindx == 3 %remove link after (out of) selected node
            s = pickpoint(xys,tidx,cumtimes,'source'); %source node
            if outdegree(G,s) == 1
                t = successors(G,s); %target node
                G = rmedge(G,s,t); %remove edge between these nodes
                updateflag = true;
            else
                disp('Node has multiple edges out, pick one of these to remove a link to instead')
            end
            clear s t
        end
    elseif strcmp(key,'n')
        %pick second point for link addition (always pick source then
        %target)
        if exist('s','var') == 1
            %pick second point for added tracking link
            t = pickpoint(xys,tidx,cumtimes,'target'); %target node
            %check if (1) the second point comes after the first in time; (2)
            %if the first point has more than one link out already, (3) if the
            %second point has more than 0 links in already
            if G.Nodes.frame(t) > G.Nodes.frame(s) && outdegree(G,s) < 2 && indegree(G,t) == 0
                EdgeTable = table([s t],1,1,'VariableNames',...
                    {'EndNodes','Weight','Cost'});
                G = addedge(G,EdgeTable);
                updateflag = true;
            else
                disp('Could not add link')
            end
        else
            disp('There is no source node chosen')
        end
    elseif strcmp(key,'c')
        %cancel link addition
        clear s
    elseif strcmp(key,'v')
        %validate a cell's track
        inds = find(~verified{pidx});
        [lindx, ~] = listdlg('ListString',cellstr(num2str(inds)),...
            'SelectionMode','single','PromptString','Choose a cell to validate');
        disp(inds(lindx))
        verified{pidx}(inds(lindx)) = true;
    elseif strcmp(key,'u')
        %unvalidate a cell's track
        inds = find(verified{pidx});
        [lindx, ~] = listdlg('ListString',cellstr(num2str(inds)),...
            'SelectionMode','single','PromptString','Choose a cell to unvalidate');
        disp(inds(lindx))
        verified{pidx}(inds(lindx)) = false;
    end
    
    if updateflag
        %if tracking digraph was modified, rebuild tracks and reinitialize
        %array of lines
        delete(p)
        xlabel('Rebuilding tracks')
        drawnow
        tracks = graphtracks(G, ntime);
        ntracks = length(tracks);
        
        [X,Y] = trackXY(tracks,allpoints);
        colors = make_colors(ntracks);
        p = initialize_lines(tracks, colors);
        xlabel('')
        
        updateflag = false; %reset update flag
    end
    
    if tidx > ntime
        tidx = ntime;
    elseif tidx < 1
        tidx = 1;
    end
    
    if tidxprev == ntime
        for ii = 1:numel(tt)
            delete(tt{ii})
        end
        clear tt
    end
    
    drawnow limitrate
end
close(f)

positions(pidx).G = G;
save(fullfile(dataDir,'validatedTracking.mat'),'verified');
save(fullfile(dataDir,'positions.mat'),'positions');

nverified = sum(verified{pidx});
ntotal = numel(verified{pidx});

fprintf('%d verified of %d = %.3g%%\n',nverified,ntotal,100*nverified/ntotal)

%% preliminary analysis of validated results
fields = {'XY','NCratio','nucLevel','cytLevel','nucArea','nucZ',...
    'nucMajorAxis','nucMinorAxis'};
histories = graphSignalingHistories(positions(pidx),G,fields,find(verified{pidx}));
% histories = graphSignalingHistories(positions(pidx),G,fields,1:numel(verified{pidx}));
start_times = cellfun(@(x) x(1), {histories.Time});
hists = histories(start_times==1);
fprintf('number of histories = %d\n',length(hists))
first_ids = cellfun(@(x) x(1), {hists.CellIdxs});
unique_ids = unique(first_ids);

ndaughters = zeros(size(unique_ids));
for ii = 1:length(unique_ids)
    ndaughters(ii) = sum(first_ids == unique_ids(ii));
end

close all
x = cellfun(@(x) x(end,1),{hists.XY});
y = cellfun(@(x) x(end,2),{hists.XY});
imshow(imgs(:,:,end))
hold on

colors = distinguishable_colors(length(unique_ids),{'k','w'});
for ii = 1:length(unique_ids)
    ids = first_ids == unique_ids(ii);
    scatter(x(ids),y(ids),50,colors(ii,:),'filled')
end

figure
histogram(ndaughters)

%% local functions
function colors = make_colors(ntracks)
%Define an array of colors for visualization
%Distuinguishable colors function makes at most 9000 different colors; if
%there are > 9000 tracks, start repeating colors
    if ntracks <= 9000
        colors = distinguishable_colors(ntracks,{'w','k'});
    else
        colors = distinguishable_colors(9000);
        numreps = ceil(ntracks/9000) - 1;
        remainder = mod(ntracks,9000); % = ntracks - numreps*9000
        colors = vertcat(repmat(colors,numreps), colors(1:remainder,:));
    end
end

function [X,Y] = trackXY(tracks,allpoints)
    X = cell(size(tracks));
    Y = cell(size(tracks));
    for ti = 1:length(tracks)
        X{ti} = NaN(size(tracks{ti}));
        Y{ti} = NaN(size(tracks{ti}));
        nonans = ~isnan(tracks{ti});
        X{ti}(nonans) = allpoints(tracks{ti}(nonans),1);
        Y{ti}(nonans) = allpoints(tracks{ti}(nonans),2);
    end
end

function p = initialize_lines(tracks, colors)
%initialize an array of Line objects to display tracking tails
    ntotal = sum(cellfun(@(x) size(x,2),tracks));
    ntracks = length(tracks);
    p = gobjects(ntotal,1);
    idx = 1;
    for ii = 1:ntracks
        for k = 1:size(tracks{ii},2)
            p(idx) = line(gca,1,1,'Color',colors(ii,:),'LineWidth',1.5);
            idx = idx + 1;
        end
    end
end

function allidx = pickpoint(xys,tidx,cumtimes,pointtype)

xlabel(strcat("Pick ", pointtype, " point"))
h = drawpoint;
xy = h.Position; %position of drawn point
delete(h)
%find the nearest cell centroid to the drawn point
D = pdist2(xys, xy(1,:));
[~, cellidx] = min(D);
%find the overall index of this point in the tracking digraph
allidx = cumtimes(tidx) + cellidx;
xlabel('')

end


