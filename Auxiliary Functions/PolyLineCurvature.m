function [Ka,K]=PolyLineCurvature(C)
% Estimate average absolute curvature of a 2D piecewise linear contour.
%
% INPUT:
%   - C     : N-by-2 array of contour coordinates, where N is the number
%             of contour vertices. For closed contours, C(1,:) must be
%             equal to C(end,:). If C(1,:)~=C(end,:), contour will be 
%             treated as an open contour. 
%
% OUTPUT:
%   - Ka    : average absolute curvature of C. For example, open contour 
%             composed of collinear line segments has Ka=0. You can also 
%             verify that a unit circle sampled at equal arc length 
%             intervals will have have Ka=1, as expected:
%   
%             figure('color','w')
%             hold on
%             N=10:10:500;
%             Ka=zeros(size(N));
%             for i=1:numel(N)
%                 t=linspace(0,1,N(i)+1);
%                 t(end)=0;
%                 t=2*pi*t(:);
%                 Ka(i)=PolyLineCurvature([cos(t) sin(t)]);
%                 plot(N(i),Ka(i),'.','MarkerSize',15)
%             end             
%             xlabel('# point samples','FontSize',25)
%             ylabel('mean abs. curvature','FontSize',25)
%             drawnow
%
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: Mar.2014
%


% Check the input
if ~isnumeric(C) || ~ismatrix(C) || size(C,2)~=2 || sum(~isfinite(C(:)))>0
    error('Contour specified incorrectly')
end

N=size(C,1);
if N==2, Ka=0; return; end

% Triangles
chk_open=false;
if norm(C(1,:)-C(end,:))<=1E-6 % closed 
    N=N-1;    
    if N==2, Ka=NaN; return; end
    v_id=1:N;
    F=cat(1,v_id,circshift(v_id,-1),circshift(v_id,-2))';
    C(end,:)=[];
else % open
    v_id=1:(N-2);
    F=cat(1,v_id,v_id+1,v_id+2)';
    chk_open=true;
end

% Curvature at the vertices
a=C(F(:,2),:)-C(F(:,1),:); 
b=C(F(:,3),:)-C(F(:,1),:);
c=C(F(:,3),:)-C(F(:,2),:);

A=abs(a(:,1).*b(:,2)-a(:,2).*b(:,1))/2;
a=sqrt(sum(a.^2,2));
b=sqrt(sum(b.^2,2));
c=sqrt(sum(c.^2,2));

K=4*A./(a.*b.*c);

% Vertex weights
W=a+c;

% Remove degenerate triangles, if any
idx=~isfinite(K);
Ko=K;
K(idx)=0;
W(idx)=0;

% Mean absolute curvature
Ka=(W'*K)/sum(W);

if nargout>1
    K=Ko;
    if chk_open
        K=cat(1,NaN,K,NaN);
    else
        K=cat(1,K,K(1));
    end
end
