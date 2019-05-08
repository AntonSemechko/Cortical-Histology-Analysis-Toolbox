function [D,P]=Pt2ContourDistance(C,X)
% Compute closest distance from a set of points to a 2D (open or closed)
% contour.
%
% INPUT:
%   - C     : N-by-2 array of contour coordinates; x- and y-coordinates are
%             contained in the first and second columns of C, respectively.
%   - X     : M-by-2 array of 2D point coordinates.
%
% OUTPUT:
%   - D     : M-by-1 array of closest distances.
%   - P     : M-by-2 array of closest points on the contour, such that
%             D(i)=norm(X(i,:)-P(i,:)). It is important to note that the
%             closest point may not be unique as there may be other points
%             along the contour that have the same distance to the query
%             point in X. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: Mar.2014
%


% Is the contour open or closed?
flag=false;
if norm(C(end,:)-C(1,:))<1E-6
    flag=true;
    C(end,:)=[];
end

% Edge length
D12=circshift(C,[-1 0])-C;    
if ~flag, D12(end,:)=[]; end

L=sqrt(sum(D12.^2,2));
N12=bsxfun(@rdivide,D12,L);

% Loop through the edges and vertices
DF=Inf*ones(size(X,1),1);
P=zeros(size(X,1),2);
for i=1:size(C,1)
    
    % Compute the distances from all points to the current vertex
    D=bsxfun(@minus,X,C(i,:));
    d2=sum(D.^2,2);
    
    % Compare the distances
    idx=d2<DF;
    d2=d2(idx);
    DF(idx)=d2;
    P(idx,1)=C(i,1);
    P(idx,2)=C(i,2);
    if i==size(C,1), break; end
    
    % Find points inside the edge characteristic
    t=bsxfun(@minus,X,C(i,:))*N12(i,:)';
    chk=t>=0 & t<=L(i);

    % Compute distances to points inside characteristic
    t=t(chk);
    if ~isempty(t)
    
        D=bsxfun(@times,t,N12(i,:));
        D=bsxfun(@plus,D,C(i,:));
        d2=sum((X(chk,:)-D).^2,2);
        
        % Compare to the existing distances
        d2_old=DF(chk,:);
        idx=d2<d2_old;
        d2=d2(idx);
        if ~isempty(d2)
            chk(chk)=idx;
            DF(chk)=d2;
            P(chk,1)=D(idx,1);
            P(chk,2)=D(idx,2);
        end
        
    end
         
end
D=sqrt(DF);

