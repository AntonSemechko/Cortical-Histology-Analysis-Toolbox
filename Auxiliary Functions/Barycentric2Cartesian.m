function P=Barycentric2Cartesian(F,X,B,f)
% Convert from barycentric back to Cartesian coordinates. 
%
%   - F     : triangulation.
%   - X     : Cartesian coordinates of the mesh vertices.
%   - B     : N-by-2 array of barycentric point coordinates.
%   - f     : N-by-1 array of face indices specifying the triangles with 
%             respect to which the barycentric coordinates are defined.
%
%   - P     : N-by-2 array of Cartesian point coordinates.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: June.2014
%


n=size(B,1);
f=f(:);
if numel(f)~=n
    error('Number of face indices does not match the number of query points')
end
F=F(f,:);
    
% Triangle edge direction vectors
D12=X(F(:,2),:)-X(F(:,1),:);
D13=X(F(:,3),:)-X(F(:,1),:);

% Cartesian coordinates
P=bsxfun(@times,B(:,1),D12)+bsxfun(@times,B(:,2),D13);
P=P+X(F(:,1),:);
