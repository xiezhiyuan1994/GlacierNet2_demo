function[imt2]=terminus_improving(imtb,imtb2,im_c,trainedNet)
% terminus improving function
% Input:
% imtb: GlacierNet output
% imtb2: FusedNet output
% im_c: image stack
% trainedNet: the trained fused network

% Output:
% imt2: DCG estimations with improved termini

% Author: Zhiyuan Xie
% Year: 2022
% GlacierNet2: A Hybrid Multi-Model Learning Architecture for Alpine Glacier Mapping

T1=30;
T2=250;

imtb_final=imtb2;
%  give each DCG a unique codename(fused network output)
imtl2=bwlabel(imtb2);


fs=9;
im_ele=imfilter(im_c(:,:,12), ones(fs)/(fs^2), 'symmetric');

imtl = bwlabel(imtb);

im_zo=zeros(size(imtb));

%find termini
t=0.15;

c=bwconncomp(imtb);

cc=c.PixelIdxList;

for i=1:length(cc)
    
    ind1=cc{i};% DCG index
    
    [ele_s ele_ind]=sort(im_ele(ind1),'ascend');
    
    low_ele=ind1(ele_ind(1:round(length(ele_s)*t)));
    
    im_zo(low_ele)=1;
    
    
end




imtl_term=imtl.*im_zo;
% crop and then correct termini 
box_space=0.2;

IOU_T=0.7;

clear i_record
num=1;




for i=1:max(imtl(:))


[h l]=find(imtl_term==i);
% estimate the crop box 
h_dis=max(h)-min(h);

l_dis=max(l)-min(l);

h_min_space=min(h)-round(h_dis*box_space);

h_min_space=max(h_min_space,1);

l_min_space=min(l)-round(l_dis*box_space);

l_min_space=max(l_min_space,1);

l_max_space=min(l_min_space+round(l_dis*(1+2*box_space))-1,size(imtl_term,2));

h_max_space=min(h_min_space+round(h_dis*(1+2*box_space))-1,size(imtl_term,1));

box=[l_min_space h_min_space l_max_space-l_min_space  h_max_space-h_min_space];


imtl_term_c=imcrop(imtl,box);

imtl_term_c(~(imtl_term_c==i))=0;


imtl_term_c2=imcrop(imtl2,box);

i2=unique(imtl_term_c2(and(imtl_term_c>0,imtl_term_c2>0)));

if ~isempty(i2)

imtl_term_c2(~ismember(imtl_term_c2,i2))=0;


IOU=sum(sum(and(imtl_term_c,imtl_term_c2)))/sum(sum(or(imtl_term_c,imtl_term_c2)));
% check IOU to find disagreements
if IOU>0&IOU<IOU_T
    
    i_record(num)=i;
    num=num+1;
    
end

end



end
% crop the termini that may have misclassifications
winmin=512;

winmin=winmin-1;

count=1;

for i=i_record
    
 count/length(i_record)*100   
 
[h l]=find(imtl_term==i);

h_dis=max(h)-min(h);

l_dis=max(l)-min(l);

h_min_space=min(h)-round(h_dis*box_space);

h_min_space=max(h_min_space,1);

l_min_space=min(l)-round(l_dis*box_space);

l_min_space=max(l_min_space,1);

l_max_space=min(l_min_space+round(l_dis*(1+2*box_space))-1,size(imtl_term,2));

h_max_space=min(h_min_space+round(h_dis*(1+2*box_space))-1,size(imtl_term,1));

box=[l_min_space h_min_space l_max_space-l_min_space  h_max_space-h_min_space];

box_org=box;

if box(3)<winmin | box(4)<winmin
    
    if box(2)==1
        
       box(4)=winmin;
       
    end
    
    if box(1)==1
        
       box(3)=winmin;
       
    end
    
    if box(1)+box(3)==size(imtl_term,2)
        
        box(1)=size(imtl_term,2)-winmin;
        
        box(3)=winmin;
        
    end
    
    if box(2)+box(4)==size(imtl_term,1)
        
        box(2)=size(imtl_term,1)-winmin;
        
        box(4)=winmin;
    end
    
    if box(3)<winmin
        
       expan_dis=round((winmin-box(3))/2);
       
       box(1)=max(box(1)-expan_dis,1);
       
       box(3)=winmin;
      
        
        if box(1)+box(3)>size(imtl_term,2)

            box(1)=size(imtl_term,2)-winmin;

        end
       
        
    end
    
    if box(4)<winmin
        
       expan_dis=round((winmin-box(4))/2);
       
       box(2)=max(box(2)-expan_dis,1);
       
       box(4)=winmin;
        
       if box(2)+box(4)>size(imtl_term,1)

            box(2)=size(imtl_term,1)-winmin;

        end
       
        
    end
    
    %odd or even

    
    if rem(box(3), 2) == 0
        
       box(3)=box(3)-1; 
      
        
    end
    
    if rem(box(4), 2) == 0
        
       box(4)=box(4)-1; 
     
    end
    
    
    
end

clear im_term_allch
% cilp the subimage according to the box
for Q=1:size(im_c,3)

im_term_allch(:,:,Q)=imcrop(im_c(:,:,Q),box);

end


clear fmap
% extract feature map via GlacierNet
if size(im_term_allch,1)<=winmin+1 & size(im_term_allch,2)<=winmin+1


fmap=activations(trainedNet,im_term_allch,'encoder1_relu_3');

else
    
    for QQ=1:winmin+1:size(im_term_allch,1)
        
         for PP=1:winmin+1:size(im_term_allch,2)
             
        Mxend=min(QQ+winmin,size(im_term_allch,1));
        Mx=Mxend-winmin;
        
        Myend=min(PP+winmin,size(im_term_allch,2));
        My=Myend-winmin;
        
            
        M1=im_term_allch(Mx:Mxend,My:Myend,:);
        
        fmap(Mx:Mxend,My:Myend,:)=activations(trainedNet,M1,'encoder1_relu_3');
        
  
        
         end
        
        
    end
    

end

fmap=double(fmap)./double(mean(mean(fmap)));

fmap(isnan(fmap))=0;

% create reference for KNN
imtl_term_c=imcrop(imtl,box);

imtl_term_c(~(imtl_term_c==i))=0;

imtl_term_c2=imcrop(imtl2,box);

i2=unique(imtl_term_c2(and(imtl_term_c>0,imtl_term_c2>0)));

imtl_term_c2(~ismember(imtl_term_c2,i2))=0;

imtl_term_c_or=or(imtl_term_c>0,imtl_term_c2>0);

imtl_term_c_and=and(imtl_term_c>0,imtl_term_c2>0);

imtl_term_c_xor=xor(imtl_term_c>0,imtl_term_c2>0);

se = strel('disk',10); 

fse=double(se.Neighborhood);

im_near=xor(imfilter(imtl_term_c_or,fse)>0,imtl_term_c_or);




fmap_rs=reshape(fmap,[size(fmap,1)*size(fmap,2),size(fmap,3)]);

imtl_term_c_and_rs=reshape(imtl_term_c_and,[size(imtl_term_c_and,1)*size(imtl_term_c_and,2),1]);

data_pos=fmap_rs(imtl_term_c_and(:),:);

data_ng=fmap_rs(im_near(:),:);

data=[data_pos;data_ng];

data_label=categorical([ones(size(data_pos,1),1);zeros(size(data_ng,1),1)]);


% KNN

%apply KNN
Mdl=fitcknn(data,data_label,'NumNeighbors',5);

X=fmap_rs(imtl_term_c_xor(:),:);

Yfit = double(predict(Mdl,X))-1;




Y=zeros(size(imtl_term_c_xor));

Y(imtl_term_c_xor)=Yfit;



Y2=bwareaopen(Y,150);

se = strel('disk',10);

Y3 = imclose(Y2,se);

imtl_term_only_c=imcrop(imtl_term,box);

imtl_term_only_c(~(imtl_term_only_c==i))=0;

Y3=and(Y3,imtl_term_only_c);
% restore the corrected terminus back
imtb_final(box(2):box(2)+box(4),box(1):box(1)+box(3))=or(imtb_final(box(2):box(2)+box(4),box(1):box(1)+box(3)),Y3);

count=count+1;

end
    

NDVI=(double(im_c(:,:,5))-double(im_c(:,:,4)))./(double(im_c(:,:,5))+double(im_c(:,:,4)));

NDVI_n=NDVI-min(NDVI(:));

NDVI_n=NDVI_n./max(NDVI_n(:));

NDVI_t=NDVI_n>0.7;

NDVI_t_imt=and(NDVI_t,imtb_final);
            
 
imtb_final2=imtb_final;

imtb_final2(NDVI_t_imt==1)=0;

se = strel('disk',3);

imtb_final2= ~imclose(~imtb_final2,se);

imtb_final2=holefill_1(imtb_final2,im_c,T1,T2);
            

% adjust snowline
imt_diff=zeros(size(imtb));

imtl2=bwlabel(imtb_final2);

t=0.2;

for R=1:max(imtl2(:))
    
    ind=find(imtl2==R);
    
    imtl_uni=unique(imtl(ind));
    
    imtl_uni(imtl_uni==0)=[];
    
    dcg1=ismember(imtl,imtl_uni);
    
    dcg2=ismember(imtl2,R);
    
    dcg2_diff=and(xor(dcg1,dcg2),dcg1);
    
    c=bwconncomp(dcg2_diff);

    cc=c.PixelIdxList;
    
    [ele_s ele_ind]=sort(im_ele(ind),'descend');
    
    high_ele_T=ele_s(round(length(ele_s)*t));
    
    for R2=1:length(cc)
        
        ind2=cc{R2};
        
       if mean(im_ele(ind2))>=high_ele_T
           
           imt_diff(ind2)=1;
           
       elseif length(ind2)>1250
           
          if sum(im_zo(ind2))==0
                 
           imt_diff(ind2)=1;
           
          end
       
           
       end
       

        
    end
   
end

imt2=or(imtb_final2,imt_diff);

end

