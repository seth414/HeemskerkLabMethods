function [barefname, pidx] = parseFilename(fname,dataDir)

[path,filename,ext] = fileparts(fname);

if ~exist('dataDir','var') || isempty(dataDir)
    dataDir = path;
end

fname = filename; %filename without extension
filename = [filename,ext]; %put the extension back on for comparisons below

if strcmp(ext,'.ims') || strcmp(ext,'.tif') || strcmp(ext,'.vsi')
    s = strsplit(filename,{'_',ext});
    teststr = s{end-1};
    if strcmp(teststr(1),'F')
        numstr = teststr(2:end);
        pidx = str2double(numstr);
        if isnan(pidx)
            pidx = [];
            suff = ext;
        else
            suff = ['_',teststr,ext];
        end
    else
        pidx = [];
        suff = ext;
    end
    
    s1 = strsplit(filename,suff);
    barefname = s1{1};
%     [splitfname, matches] = strsplit(filename,{'_F[0-9][0-9]*', ['\' ext]},...% F->. for any letter
%         'DelimiterType','RegularExpression',...
%         'CollapseDelimiters',false);
%     barefname = splitfname{1};
%     if length(matches) == 2
%         pidx = str2double(matches{1}(3:end));
%     elseif length(matches) == 1
%         pidx = [];
%     else
%         error('how could there be more than 2 matches?')
%     end
elseif strcmp(ext,'.nd2')
    r = bfGetReader(fullfile(dataDir,filename));
    ns = r.getSeriesCount;
    if ns == 1
        [splitfname, matches] = strsplit(filename,{['_[0-9][0-9][0-9][0-9]',ext]},...
            'DelimiterType','RegularExpression',...
            'CollapseDelimiters',false);
        if length(matches) == 1
            %if there is a matching number at the end of the filename
            pidx = str2double(matches{1}(2:end-(length(ext))));
            barefname = splitfname{1};
        elseif isempty(matches)
            %if there isn't a number matching the above pattern at the end
            %of the filename, it could either be first in the list of
            %numbered files or by itself (not part of a list)
            %look for the next file in the list
            barefname = fname;
            if exist(fullfile(dataDir,[fname,'_0001',ext]),'file')==2
                pidx = 0;
            else
                pidx = [];
            end
        end
    else
        %multiple series nd2 file
        barefname = fname;
        pidx = [];
    end
elseif strcmp(ext,'.lif')
    %multiple series lif file
    barefname = fname;
    pidx = [];
end



end