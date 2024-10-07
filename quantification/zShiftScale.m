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