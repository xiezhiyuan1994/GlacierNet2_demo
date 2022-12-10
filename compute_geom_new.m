function[geom]=compute_geom_new(Z,r,cellsize)


% x: left to right
% y: top to bot

fz=ones(r,r);

fzx=repmat(([1:r]-0.5).*cellsize,[r,1])-r/2*cellsize;

fzy= repmat(([1:r]'-0.5).*cellsize,[1,r])-r/2*cellsize;

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


% s(:,:,1)=s_z;
% s(:,:,2)=s_zx;
% s(:,:,3)=s_zy;
% s(:,:,4)=s_zxy;
% s(:,:,5)=s_zx2;
% s(:,:,6)=s_zy2;
% s(:,:,7)=s_x2y2;
% s(:,:,8)=s_x2;
% s(:,:,9)=s_x4;
% 
% s2=reshape(s,[size(s,1)*size(s,2) size(s,3)]);

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

pq1 = 1.0 + p.*p + q.*q;
qp = sqrt((1.0 + q.*q) ./ (1.0 + p.*p));

a = r .* qp - t ./ qp;
b = p.*q.*r .* qp ...
 - 2.0.* s.* sqrt((1.0 + q.*q) .* (1.0 + p.*p)) ...
 + p.*q.*t .* qp;

pq2 = p.*p + q.*q;
    
curv_tangential=-(q.*q .* r - 2.0 .* p.*q.*s + p.*p .* t) ./ (pq2 .* sqrt(pq1));

curv_profile=-(p.*p .* r + 2.0 .* p.*q.*r.*s + q.*q .* t) ./...
        ((p.*p + q.*q) .* sqrt(pq1 .* pq1 .* pq1));
    
slope=atand(sqrt(p.*p + q.*q));

uns= -sqrt(a.*a ./ pq1 + b.*b) ./ (2.0 .* sqrt(pq1.*pq1.*pq1));

clear pq1 qp a pq2

unsphericity=reshape(uns,size(Z));

curv_tangential=reshape(curv_tangential,size(Z));

curv_profile=reshape(curv_profile,size(Z));

slope=reshape(slope,size(Z));


geom.unsphericity=unsphericity;

clear unsphericity

geom.curv_tangential=curv_tangential;

clear curv_tangential

geom.curv_profile=curv_profile;

clear curv_profile

geom.slope=slope;

clear slope Z


end

