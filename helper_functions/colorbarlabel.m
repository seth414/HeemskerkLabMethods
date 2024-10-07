function h = colorbarlabel(labelstring,limitflag)
%add a label to a colorbar defined by the string input labelstring
%limitflag is an optional input to determine whether or not to show
%numerical limits on the ends of the colorbar; default is true

if ~exist('limitflag','var')
    limitflag = true;
end

h = colorbar; h.Label.String = labelstring;
h.Ticks = {};

if limitflag
    hlabelpos = h.Label.Position;
    h.Ticks = h.Limits;
    h.Label.Position = hlabelpos;
end

end