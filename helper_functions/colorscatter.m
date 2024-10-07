function colorscatter(x,y,z,tol,pointsize,orderbycolor)

if ~exist('tol','var')
    tol = 0.02;
end

if ~exist('pointsize','var')
    pointsize = 15;
end

if ~exist('orderbycolor','var')
    orderbycolor = false;
end

if isstring(orderbycolor) || ischar(orderbycolor)
    orderdir = orderbycolor;
    orderbycolor = true;
elseif orderbycolor
    orderdir = 'ascend';
end

if orderbycolor
    [~,I] = sort(z,orderdir);
    x = x(I); y = y(I); z = z(I);
end

nnz = z(~isnan(z));
n = length(nnz);
zs = sort(nnz);
zmin = zs(max(ceil(n*tol),1));
zmax = zs(floor(n*(1-tol)));
z(z < zmin) = zmin; z(z > zmax) = zmax;

scatter(x,y,pointsize,z,'filled')
colormap turbo
colorbar
cleanSubplot
% axis square

end