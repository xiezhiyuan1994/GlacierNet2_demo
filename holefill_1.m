function[out]=holefill_1(imt2,im,T1,T2)
%
%
%



imt3=1-imt2;

c=bwconncomp(imt3);

cc=c.PixelIdxList;

ims=im(:,:,16);

for i=1:length(cc)
    
wr=ims(cc{i});
   
deg(i)=mean(abs(wr))./(2^16-1)*90;
wr_length(i)=length(wr);   

end

[maxv maxi]=max(wr_length);

hole_f=find(deg<T1);

hole_f2=find(wr_length<T2);

hole_f3=intersect(hole_f,hole_f2);

for q=1:length(hole_f3)
 
   imt2(cc{hole_f3(q)})=1;
   
   
end
 
    out=imt2;
   





end