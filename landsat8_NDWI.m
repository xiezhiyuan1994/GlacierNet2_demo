function[im_NDWI]=landsat8_NDWI(im,band1,band2)


b1=double(im(:,:,band1));
b2=double(im(:,:,band2));

im_NDWI=(b1-b2)./(b1+b2);

imn=im_NDWI-min(im_NDWI(:));
im_NDWI=imn./max(imn(:));

im_NDWI(isnan(im_NDWI))=0;


end