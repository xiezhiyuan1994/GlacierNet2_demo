function lb7=map_SCG(DEMf,imt2,im_c)
% Input:
% DEMf: filled DEM
% imt2: DCG
% im_c: image stack


% Output:
% lb7: DCG and SCAZ

% Author: Zhiyuan Xie
% Year: 2022
% GlacierNet2: A Hybrid Multi-Model Learning Architecture for Alpine Glacier Mapping
% snow covered accumulation zone mapping function

addpath(genpath('D:\snow_covered_glacier_detect\topotoolbox-master\topotoolbox-master'));


% 
Z=double(DEMf.Z);

% to match the DEM with 30m resolution,
bwr2=imresize(imt2,0.5);

lb=bwlabel(bwr2);
lb_copy=lb;
% 1)	The surface reflectance snow index is used to identify all snow-covered pixels within the scene
[im_NDWI]=landsat8_NDWI(im_c,3,6);

imw=zeros(size(im_NDWI));

imw(im_NDWI>0.7)=1;

% 2) Snowy regions are removed if they are not adjacent to and overlapping with the ablation zone.
% 3)	The remaining snowy regions are merged with adjacent or overlapping ablation zone.
[im3]=relative_region_merge(imw,imt2);

% 4)	All gaps in the data are filled, representing the internal, 
% non-snow and non-glacier regions to assist with computational processing.
im3f=imfill(im3,'holes');

% reduce the size
im3f2=imresize(im3f,0.5,'nearest');

%make the pixel outside the merged region to NaN
lb2=lb;

DEMf2=DEMf;

DEMf2.Z=DEMf2.Z.*im3f2;

DEMf2.Z(DEMf2.Z==0)=NaN;

% compute FLOWobj

FD = FLOWobj(DEMf2);

FD2=FD;

%A = flowacc(FD);
% 5)	Ablation zones are defined by different code numbers.
c=bwconncomp(bwr2);

cc=c.PixelIdxList;
% 6)	Pixels are allocated to each target ablation zone, and an intermediate output G-1 was produced, 
% indicating the DBs of glaciers within the merged region. 
%Each DB subsequently inherits its corresponding ablation zone’s code number.
for i=1:length(cc)
    
    ind1=cc{i};% DCG index
    

        L = drainagebasins(FD,ind1);
        
        ind5=find(L.Z>0);
        
        ud=unique(lb(ind5));
        
        ud(ud==0)=[];
        
        if length(ud)>1
        
        lb(ismember(lb,ud))=i;
        
        end
        
        lb(ind5)=i;

        
    %end
   
    
    
end

% 7)	The DB for each target ablation zone is computed across the entire scene and also inherits the corresponding ablation zone. 

FD = FLOWobj(DEMf);

for i=1:length(cc)
ind1=cc{i};

L = drainagebasins(FD,ind1);

dba{i,1}=ind1;

dba{i,2}=find(L.Z>0);
    
end

dba2=dba;
% 
for i = 1:length(cc)

db_inter=cellfun(@intersect,repmat(dba(i,1),[length(cc),1]),dba(:,2),'UniformOutput',false);
nonempty=find(~cellfun(@isempty,db_inter));
nonempty=setdiff(nonempty,i);

if length(nonempty)>0
  
   
  for  p=1:length(nonempty)
    
  dba2(nonempty(p),2)={setdiff( dba{nonempty(p),2},  dba{i,2} )};
    
  end
    
end
    
end

for i=1:length(cc)
    lb2(dba2{i,2})=i;
    
   
end
% Then, the AND operator removes DB pixels located out of the merged region and obtains another intermediate output G-3.
lb2=lb2.*im3f2;

% 8)	The intermediate G-2 is computed, 
% representing the DB segment inside the merged region without a specified target area. 
db_all=drainagebasins(FD2);
% 9)	The difference between G-1 and G-3 is calculated
% 10)	The segments of G-2 corresponding to the differences are marked if it neither neighbored nor overlapped with G-1 DBs of the same code number.

db_all_z=db_all.Z;

db_all_z(db_all_z>0)=db_all_z(db_all_z>0)+max(lb(:));

db_all_z(lb>0)=lb(lb>0);

f_db_all_z = ordfilt2(double(db_all_z),1,ones(3,3));


aa=unique(db_all_z(f_db_all_z==lb2));

diff1=setdiff(unique(db_all_z) ,aa );

%db_all_z_copy=db_all.Z;

lb4=lb2; % lb4 
% 11)	The marked DB segments in G-2 are removed corresponding to the area in G-3.
lb4(ismember(db_all_z,diff1))=0;

lb4_2=bwareaopen(lb4,50);

lb4=lb4.*double(lb4_2);

% covert the db that index belong to lb to lb
db_all_z=db_all.Z;

gg=unique(db_all_z(lb>0));

gg(gg==0)=[];

lb5=double(ismember(db_all_z,gg)).*lb4;
%12)	The difference between G-1 and modified G-3 is recalculated. 
lb_diff=lb4-lb5; % 

bwx=boundarymask(lb_diff);




% 13)	The segments of G-2 corresponding to the recalculated differences are identified to locate each segment’s connecting borderlines. 
% These connections are part of the segment boundary directly touching the DBs from G-1, with the same segment code number. 
% The segments are then selected if the average elevation of its borderline was higher than the mean elevation of the entire segment boundary.

% junction point
tbp=(double(f_db_all_z==lb2)-double(lb>0)).*double(lb2>0)>0; %

db_all_z=db_all.Z;

lb_diff_db_b=double(bwx).*double(db_all_z);

lb_diff_db_jp=double(tbp).*double(db_all_z);


Zhigh=DEMf.Z;

Zhigh=Zhigh.*double(bwx);

lb_diff_b_ind=unique(lb_diff_db_b(:));

lb_diff_b_ind(lb_diff_b_ind==0)=[];


lb_diff_jp_ind=unique(lb_diff_db_jp(:));

lb_diff_jp_ind(lb_diff_jp_ind==0)=[];

clear kpp

ww=1;

for i=1:length(lb_diff_jp_ind)
    
    i
    
    mm=lb_diff_jp_ind(i);
    
    mm1=mean(Zhigh(find(lb_diff_db_jp==mm)));
    
    mm2=mean(Zhigh(find(lb_diff_db_b==mm)));
    
    %sort(Zhigh(find(lb_diff_db_b==mm)),'')
    
    if mm1<=mm2
       kpp(ww)=mm; 
        
       ww=ww+1; 
    end
    
    
    
    
end

lbb=double(ismember(db_all_z,kpp)).*lb4;
% lb4 was removed most negative case, now further check by elevation value.

% 14)	The selected segments in G-2 corresponding to the area in G-3 are again removed.
lb6=lb5;

lb6(lbb>0)=lbb(lbb>0);

% 15)	Morphological close operations are applied to merge the filled gaps from Step 4. Accordingly, all holes are restored.

se = strel('disk',6);

holes=xor(im3f,im3);

holes=imfill(holes,'holes');

holes = imclose(holes,se);


im3=im3f;

im3(holes)=0;

im3r2=imresize(im3,0.5,'nearest');



lb6=lb6.*double(im3r2>0);

indx=unique(lb6(:));

indx(indx==0)=[];

lb6_new=zeros(size(lb6));

% 16)	Isolated SCAZ regions caused by hole restoration are removed. 
% All isolated regions were independent pixel clusters, unconnected with the corresponding primary glacier region.

for i=1:size(indx,1)
    
  indx1=indx(i);
    
  imss=lb6==indx1;
  
 largestBlob = bwareafilt(imss, 1);
 largestArea = sum(largestBlob(:));
 lb6_new(bwareaopen(imss,largestArea-1)==1)=indx1;
    
    
end

lb6=lb6_new;



% in step 15) remove small holes

T=30;

indx=unique(lb6(:));

indx(indx==0)=[];

lb7=lb6;


for i=1:size(indx,1)
    
    indx1=indx(i);
    
    imss=double(lb6==indx1);
    
    imss2=1-bwareaopen(1-imss,T);
    
    lb7(logical(imss2-imss))=indx1;
    
end

end