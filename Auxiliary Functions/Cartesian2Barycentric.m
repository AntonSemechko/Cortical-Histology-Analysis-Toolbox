function [B,f]=Cartesian2Barycentric(F,X,P,f)
% Given a set of 2D points, compute their barycentric coordinates with
% respect to the underlying triangulation. 
%
%   - F     : triangulation.
%   - X     : Cartesian coordinates of the mesh vertices.
%   - P     : N-by-2 array of Cartesian point coordinates.
%   - f     : N-by-1 array of triangle indices such that f(i) and P(i,:)
%             are corresponding. If f is not know, specify it as an empty 
%             array and f will be determined automatically.  
%
%   - B     : N-by-2 array of barycentric point coordinates.
%   - f     : N-by-1 array of face indices specifying the triangles with 
%             respect to which the barycentric coordinates are defined.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: June.2014
%


if nargin<4, f=[]; end

n=size(P,1);
if ~isempty(f)
    f=f(:);
    if numel(f)~=n
        error('Number of face indices does not match the number of query points')
    end
    F=F(f,:); 
end

% Triangle edge direction vectors
D12=X(F(:,2),:)-X(F(:,1),:);
D13=X(F(:,3),:)-X(F(:,1),:);

% Determinants
DT=D12(:,1).*D13(:,2)-D13(:,1).*D12(:,2);

% Find barycentric co-ordinates: f is known
if ~isempty(f)
    d=bsxfun(@minus,P,X(F(:,1),:));
    u=(d(:,1).*D13(:,2)-D13(:,1).*d(:,2))./DT;
    v=(d(:,2).*D12(:,1)-D12(:,2).*d(:,1))./DT;
    B=[u,v];
    return
end

% Find barycentric co-ordinates: f is unknown
[U,V,f]=deal(nan(n,1));
idx=1:n;
tol=1E-8;
for i=1:size(F,1)
    
    d=bsxfun(@minus,P,X(F(i,1),:));
    u=(d(:,1).*D13(i,2)-D13(i,1).*d(:,2))/DT(i);
    v=(d(:,2).*D12(i,1)-D12(i,2).*d(:,1))/DT(i);
    
    chk=(u>-tol & u<(1+tol)) & (v>-tol & v<(1+tol)) & (u+v)<(1+tol);
    
    if sum(chk)==0, continue; end
        
    U(idx(chk))=u(chk);
    V(idx(chk))=v(chk);
    f(idx(chk))=i;
    
    P(chk,:)=[];
    idx(chk)=[];
    if isempty(idx), break; end
    
end
B=[U,V];

