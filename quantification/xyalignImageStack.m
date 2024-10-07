function IMa = xyalignImageStack(img,shiftyx)
m = size(img,1); n = size(img,2); nz = size(img,3); imclass = class(img);

IMa = zeros(m,n,nz,imclass);
for zi = 1:nz
    im = img(:,:,zi);
    ima = alignImage(im,im,shiftyx);    
    IMa(:,:,zi) = ima;
end

end
