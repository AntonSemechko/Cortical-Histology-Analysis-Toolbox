function [F,P]=step3_ParameterizeContour(C,idx)
% Map closed contour to a trapezoid of unit area using area-preserving
% parameterization.
%
% INPUT:
%   - C,idx     : outputs of 'step2_ConstructClosedContour' function.
%
% OUTPUT:
%   - F         : M-by-3 face-vertex connectivity list, where M is the
%                 total numer of triangles. 
%   - P         : (N+1)-by-2 array of vertices in the parameter domain, 
%                 so that C(i,:) corresponds to P(i,:). 
%
% ******************************* IMPORTANT *******************************  
% Pia (i.e., C(idx(3):idx(4),:)) will be mapped to a vertical line at x=0 
% and white matter boundary (i.e., C(idx(1):idx(2),:)) will be mapped to a
% vertical line at x=1.
% *************************************************************************
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: Mar.2014
%


% Triangulate the contour
C(end,:)=[];
N=size(C,1);
c=[1:N-1;2:N]'; % edges
TR=delaunayTriangulation(C,c); % constrained Delaunay triangulation

% Only retain triangles whose centroids are inside the contour
F=TR.ConnectivityList; 
X=(C(F(:,1),:)+C(F(:,2),:)+C(F(:,3),:))/3; % triangle centroids

idx_in=inpoly(X,[C;C(1,:)]);
F=F(idx_in,:);

% Find triangles that only connect to a single boundary (either outer or
% inner cortex boundary) and remove them
chk1=sum(F>=idx(1) & F<=idx(2),2)==3;
chk2=sum(F>=idx(3) & F<=idx(4),2)==3;
F(chk1|chk2,:)=[];

% Re-triangulate the region enclosed by the contour
N=size(F,1);
Vidx=deal(zeros(N,1));
n=zeros(N,2);
for i=1:N
    
    E=F(i,:)';
    E=[E,circshift(E,[-1 0])]; %#ok<*AGROW>
    
    D=sum((C(E(:,1),:)-C(E(:,2),:)).^2,2);
    D=D/max(D);
    [Dmin,j]=min(D);
    
    if Dmin>0.3
        chk=F(i,:)>=idx(1) & F(i,:)<=idx(2);
        if sum(chk)~=2, chk=~chk; end
        n(i,:)=F(i,chk);
        Vidx(i)=F(i,~chk);
    else
        n(i,:)=E(j,:);
        if j==1
            Vidx(i)=F(i,3);
        elseif j==2
            Vidx(i)=F(i,1);
        else
            Vidx(i)=F(i,2);
        end
    end
    
end

n=sort(n,2);
dn=n(:,2)-n(:,1);
F(dn~=1,:)=[];

Vidx(dn==1)=[];
n(dn==1,:)=[];
for i=1:size(n,1)
    Fi=[n(i,1):n(i,2)-1;...
       n(i,1)+1:n(i,2)];
    Fi(3,:)=Vidx(i);
    F=[F;Fi'];
end

% Map triangulated contour to a trapezoid with unit area ------------------

% Find triangles attached to the "interior" edges (as opposed to the edges that comprise the boundary)
TR=triangulation(F,C);
E=edges(TR);
EA=edgeAttachments(TR,E);

N=size(E,1);
flag=false(N,1);
for i=1:N
    if numel(EA{i})==2, flag(i)=true; end
end
EA=EA(flag);
EA=cell2mat(EA);

% Find two triangles at the opposite ends
i1=find(sum(F==idx(1) | F==idx(4),2)==2,1);
i2=find(sum(F==idx(2) | F==idx(3),2)==2,1);

% Order the triangles
N=size(F,1);
f_srt=zeros(N,1);
f_srt([1 N])=[i1 i2];
for i=2:(N-1)
    ind=sum(EA==f_srt(i-1),2)>0;
    f=EA(ind,:);
    EA(ind,:)=[];
    f(f==f_srt(i-1))=[];
    f_srt(i)=f;
end   
F=F(f_srt,:);
ind=find(F(1,:)==idx(4))-1;
F(1,:)=circshift(F(1,:),[0 -ind]);

% Assign labels to the contour vertices (based on which boundary they belong to)
X=ones(size(C,1),1);
X(idx(3):idx(4))=0;

% Compute triangle areas
A=TriangleAreas({F C});
A=A/sum(A);

% Map the mesh to a trapezoid of unit area
C_new=C;
C_new(idx(4),:)=[0 0];
C_new(1,:)=[1 0];

C_new(F(1,3),:)=[X(F(1,3)) 2*A(1)];

F_flg=F;
F_flg(1,1:2)=NaN;
F_flg(F==F(1,1))=NaN;
F_flg(F==F(1,2))=NaN;
F_flg(F==F(1,3))=NaN;

N=size(F,1);
for i=2:N
    
    f=F(i,:);
    chk=isnan(F_flg(i,:));
       
    y=C_new(f(chk),2);
    x=X(f(chk));
    
    x3=X(f(~chk));
    if x3==x(1)
        y3=y(1)+2*A(i);
    else
        y3=y(2)+2*A(i);
    end    
    C_new(f(~chk),:)=[x3 y3];
    F_flg(F==f(~chk))=NaN;
        
end
P=C_new;
P=[P;P(1,:)];

close all
figure, triplot(triangulation(F,C_new)), axis equal


function TA=TriangleAreas(TR)
% Calculate area of individual triangles in a triangular mesh.
%
% INPUT:
%   - TR    : 1-by-2 cell, TR={Tri,V}, where Tri is an M-by-3 array of 
%             faces and V is an N-by-3 array of vertex coordinates.
%   - Fidx  : index list of faces for which to compute the areas.
%
% OUTPUT:
%   - TA    : M-by-1 array of triangle areas where, M is the number of 
%             triangular faces. 
%


[Tri,V]=deal(TR{1},TR{2});  
if size(V,2)==2, V(:,3)=0; end

V1=V(Tri(:,1),:);
V2=V(Tri(:,2),:);
V3=V(Tri(:,3),:);

% Direction vectors
D1=V2-V1;
D2=V3-V1;

% Areas
TA=cross(D1,D2,2);
TA=sqrt(sum(TA.^2,2))/2;

