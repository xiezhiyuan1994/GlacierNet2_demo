function[region]=relative_region_merge(a,b)

a=double(a);
b=double(b);

c=bwconncomp(a);
cc=c.PixelIdxList;

for i=1:size(cc,2)
    
  r=b(cc{i});
  if length(find(r==1))==0
    
      a(cc{i})=0;
    
      
  end
  
end

region=double(or(a,b));


end