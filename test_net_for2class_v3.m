function [testing_output]=test_net_for2class_v3(test_im,trainedNet,s,w)

% terminus improving function
% Input:
% test_im: input raster/image
% trainedNet: Trained Network
% s: window shifting stride
% w: window size

% Output:
% testing_output: network prediction of  test_im

% Author: Zhiyuan Xie
% Year: 2022
% GlacierNet2: A Hybrid Multi-Model Learning Architecture for Alpine Glacier Mapping



a=size(test_im,1);
b=size(test_im,2);

a1=ceil(a/2)*2;
b1=ceil(b/2)*2;

p1=uint16(zeros(a1,b1,size(test_im,3)));

p1(1:a,1:b,:)=test_im;

clear test_im

a2=ceil(a1/s)*s;
b2=ceil(b1/s)*s;

%imc_p=padarray(p1,[(a2-a1)/2,(b2-b1)/2],0);

imc_p=padarray(p1,[(a2-a1)/2,(b2-b1)/2],'symmetric');

clear p1;

%w=256;

h=single(ones(w,w));
%im_bg_all= zeros(size(imc_p,1),size(imc_p,2));
im_gl_all= single(zeros(size(imc_p,1),size(imc_p,2)));
% im_sn_all= zeros(size(imc_p,1),size(imc_p,2));
% im_cl_all= zeros(size(imc_p,1),size(imc_p,2));
H=single(zeros(size(imc_p,1),size(imc_p,2)));

PC=((size(imc_p,1)-w)/(s)+1)*((size(imc_p,2)-w)/(s)+1);

for i=1:((size(imc_p,1)-w)/(s)+1)

for p = 1:((size(imc_p,2)-w)/(s)+1)
    
    seg_im_sub=  single(semanticseg(  imc_p((i-1)*s+1:(i-1)*s+w,(p-1)*s+1:(p-1)*s+w,:), trainedNet));



    im_gl=zeros(size(seg_im_sub));
    im_gl(seg_im_sub==1)=1;
    im_gl_all((i-1)*s+1:(i-1)*s+w,(p-1)*s+1:(p-1)*s+w)=im_gl_all((i-1)*s+1:(i-1)*s+w,(p-1)*s+1:(p-1)*s+w)+im_gl;
    
    
    
%     im_bg=zeros(size(seg_im_sub));
%     im_bg(seg_im_sub==2)=1;
%     im_bg_all((i-1)*s+1:(i-1)*s+w,(p-1)*s+1:(p-1)*s+w)=im_bg_all((i-1)*s+1:(i-1)*s+w,(p-1)*s+1:(p-1)*s+w)+im_bg;


    
    H((i-1)*s+1:(i-1)*s+w,(p-1)*s+1:(p-1)*s+w)=H((i-1)*s+1:(i-1)*s+w,(p-1)*s+1:(p-1)*s+w)+h;


end

pecent=i*p/PC*100


end


im_gl_p=im_gl_all./H;
%im_sn_p=im_sn_all./H;
%im_cl_p=im_cl_all./H;
%im_bg_p=im_bg_all./H;

clear im_gl_all
clear im_bg_all
clear im_gl
clear im_bg

testing_output=im_gl_p((a2-a1)/2+1:(a2-a1)/2+a,(b2-b1)/2+1:(b2-b1)/2+b);


%imc_p=imc_p((a2-a1)/2+1:(a2-a1)/2+a,(b2-b1)/2+1:(b2-b1)/2+b,:);








end