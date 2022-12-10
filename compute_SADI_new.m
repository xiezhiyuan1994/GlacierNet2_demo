function[geom]=compute_SADI_new(Z,r,cellsize)

% x: left to right
% y: top to bot

fz=ones(r,r);

fzx=repmat(([1:r]-0.5).*cellsize,[r,1])-r/2*cellsize;

fzy= -(repmat(([1:r]'-0.5).*cellsize,[1,r])-r/2*cellsize);

fzxy=(([1:r]'-0.5).*cellsize-r/2*cellsize)*(([1:r]-0.5).*cellsize-r/2*cellsize);

fzx2=fzx.^2;

fzy2=fzy.^2;

s_x2y2=sum(sum(fzx2.*fzy2));

s_x2=sum(sum(fzx2));

s_x4=sum(sum(fzx2.^2));

s_z=imfilter(Z,fz,'replicate');

s_zx=imfilter(Z,fzx,'replicate');

s_zy=imfilter(Z,fzy,'replicate');

s_zxy=imfilter(Z,fzxy,'replicate');

s_zx2=imfilter(Z,fzx2,'replicate');

s_zy2=imfilter(Z,fzy2,'replicate');



M=zeros(6,6);

M(1,1)=s_x4;

M(1,2)=s_x2y2;

M(1,6)=s_x2;

M(2,1)=s_x2y2;

M(2,2)=s_x4;

M(2,6)=s_x2;

M(3,3)=s_x2y2;

M(4,4)=s_x2;

M(5,5)=s_x2;

M(6,1)=s_x2;

M(6,2)=s_x2;

M(6,6,:)=r^2;

clear s_x4 s_x2y2 s_x2

Q=zeros(6,size(Z,1)*size(Z,2));

Q(1,:)=s_zx2(:);

Q(2,:)=s_zy2(:);

Q(3,:)=s_zxy(:);

Q(4,:)=s_zx(:);

Q(5,:)=s_zy(:);

Q(6,:)=s_z(:);

clear s_zx2 s_zy2 s_zxy s_zx s_zy s_z

[L,U,P]=lu(M);

y=L\(P*Q);

c=U\y;

p=c(4,:);

q=c(5,:);

s=c(3,:);

r=c(1,:);

t=c(2,:);

clear Q M y c

aspect3=-90.*(1-sign(q)).*(1-abs(sign(p)))+180.*(1+sign(p))-180.*sign(p).*acos(-q./ ((p.^2+q.^2).^0.5) )./pi;

clear pq1 qp a pq2

aspect3=reshape(aspect3,size(Z));

geom.aspect=aspect3;

f2(:,:,1)=[0 1 0;0 0 0;0 0 0];

f2(:,:,2)=[0 0 1;0 0 0;0 0 0];

f2(:,:,3)=[0 0 0;0 0 1;0 0 0];

f2(:,:,4)=[0 0 0;0 0 0;0 0 1];

f2(:,:,5)=[0 0 0;0 0 0;0 1 0];

f2(:,:,6)=[0 0 0;0 0 0;1 0 0];

f2(:,:,7)=[0 0 0;1 0 0;0 0 0];

f2(:,:,8)=[1 0 0;0 0 0;0 0 0];

for i=1:size(f2,3)
    
    fazi(:,:,i)=imfilter(double(aspect3),f2(:,:,i));
    
end

ab=[0:45:360-45];

ab=reshape(ab,[1 1 length(ab)]);

geom.SADI=1-sum(abs(cosd(fazi)-cosd(ab))+abs(sind(fazi)-sind(ab)) , 3)/4/size(fazi,3);

end