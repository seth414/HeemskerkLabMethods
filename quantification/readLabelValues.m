function [labels,lprop] = readLabelValues(LM,masks)

ul = unique(LM(LM > 0));
nl = length(ul);
ncells = length(masks);
labels = NaN(ncells,1);
lprop = NaN(ncells,nl);
for cii = 1:ncells
    vals = LM(masks(cii).nucmask);
    for li = 1:nl
        lprop(cii,li) = sum(vals == ul(li))/numel(vals);
    end
    labels(cii) = mode(vals);
end


end