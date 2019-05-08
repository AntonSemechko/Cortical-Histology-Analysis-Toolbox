function [De,D,Dn,f]=DistBetween2Contours(C1,C2)
% Compute mean Euclidean and mean normal distances from one 2D piecewise 
% linear contour (C1) to another (C2).
%
% INPUT:
%   - C1, C2    : N1-by-2 and N2-by-2 arrays of point coordinates
%                 representing the contours. Mean distances will be
%                 computed from C1 to C2.
%
% OUTPUT:
%   - De        : mean Euclidean distance from C1 to C2.
%   - D         : vectors of point to contour distances
%   - Dn        : mean normal distance from C1 to C2. 
%   - f         : in computing Dn, some of the normal lines passing through 
%                 the vertices of C1 may not intersect C2 at any point,
%                 thus, the distances at these locations is undefined.
%                 f is a number between 0 and 1, indicating the fraction
%                 of C1 length for which Dn could be computed. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: Mar.2014
%


% Are the contours open or closed?
flag1=false; 
if norm(C1(1,:)-C1(end,:))<1E-6
    flag1=true; 
    C1(end,:)=[];
end

flag2=false;
if norm(C2(1,:)-C2(end,:))<1E-6
    flag2=true; 
    C2(end,:)=[];
end


% MEAN EUCLIDEAN DISTANCE -------------------------------------------------

% Distance from the vertices of C1 to C2
D=Pt2ContourDistance(C2,C1); 

% Edge lengths of C1
E=circshift(C1,[-1 0])-C1;
L=sqrt(sum(E.^2,2));
if ~flag1
    E(end,:)=[];
    L(end)=[]; 
end

% Estimate mean Euclidean distance from C1 to C2 using 1st order approx
D=D(:);
De=(D+circshift(D,[-1 0]))/2;
if ~flag1, De(end)=[]; end
De=sum(De.*L)/sum(L);


% MEAN NORMAL DISTANCE ----------------------------------------------------
if nargout<3, return; end

% Compute normal vectors for every vertex in C1. For open contours, the
% end-points are assigned the normals of their respective line segments.  
E=bsxfun(@rdivide,E,L);
Ne=[E(:,2),-E(:,1)];    % edge normals
Nv=Ne+circshift(Ne,[1 0]);
Nv=bsxfun(@rdivide,Nv,sqrt(sum(Nv.^2,2)));
if ~flag1
    Nv(1,:)=Ne(1,:);
    Nv=[Nv;Ne(end,:)];
end

% Define lines along the vertex normals of C1 as cos(t)*x + sin(t)*y = r; 
% [t,r] are the parameters. The normal distances between C1 and C2 will be
% measured along these lines. Main advantage of the parameterization used 
% above is that it enables a very simple way of computing points of 
% intersection (with another line). 
t=atan2(Nv(:,2),Nv(:,1));
t(t<0)=t(t<0)+2*pi; 

% t above is the angle defining the orientation of the vertex normal.
% However, the t we seek will be +90 or -90 degrees offset from this t
% depending on the relative orientation betweent the position vector of the
% vertex and the vertex normal. 
chk=sign(sum([-Nv(:,2),Nv(:,1)].*C1,2));
chk(chk==0)=1;
t(chk>0)=t(chk>0)+pi/2;
t(chk<0)=t(chk<0)-pi/2;

t(t>2*pi)=t(t>2*pi)-2*pi;
t(t<0)=t(t<0)+2*pi; 

% Now can compute the r parameter; which equals the shortest distance
% between the origin and the normal line passing through the vertex.
r=abs(cos(t).*C1(:,1)+sin(t).*C1(:,2)); 

% Finaly, can start solving for points of intersection between normal lines 
% through vertices of C1 and line segments making-up C2. 
Va=C2;
Vb=circshift(C2,[-1 0]);
if ~flag2
    Va(end,:)=[];
    Vb(end,:)=[];
end
D=Vb-Va;

m=size(C1,1);
Dn=nan(m,1);
for i=1:m
    
    % Points of intersection between the i-th vertex normal (on C1) and C2
    Num=r(i)-(cos(t(i))*Va(:,1)+sin(t(i))*Va(:,2));
    Den=cos(t(i))*D(:,1)+sin(t(i))*D(:,2);
    s=Num./Den;
    
    idx=s>=0 & s<=1;
    
    if sum(idx)==0, continue; end % there are no points of intersection 
    
    s=s(idx);
    P=Va(idx,:)+bsxfun(@times,s,D(idx,:)); % actual points on intersection
    
    % Compute (squared) distance and select the shortest
    D2=bsxfun(@minus,P,C1(i,:));
    D2=sum(D2.^2,2);
    Dn(i)=min(D2);
    
end
Dn=sqrt(Dn);

% Estimate mean normal distance from C1 to C2 using 1st order
% approximation. Keep record of how many edges were actually usable,
% howver. 
Dn=Dn(:);
Dn=(Dn+circshift(Dn,[-1 0]))/2;
if ~flag1, Dn(end)=[]; end

idx=isnan(Dn);
Dn=sum(Dn(~idx).*L(~idx))/sum(L(~idx));
f=sum(L(~idx))/sum(L);
if isempty(f)
    f=0;
    Dn=NaN;
end



% Visualize normal lines through vertices of C1 to make sure they are
% correct -----------------------------------------------------------------

% Current domain
X_min=min(min(C1,[],1),min(C2,[],1));
X_max=max(max(C1,[],1),max(C2,[],1));

% Expand the domain by 10% in all directions
D=X_max-X_min;
X_min=X_min-0.05*D;
X_max=X_max+0.05*D;

a=cos(t);
b=sin(t);

X=repmat([X_min(1) X_max(1)],[m 1]);
Y=repmat([X_min(2) X_max(2)],[m 1]);


idx=abs(a)>abs(b);
if sum(idx)>0
    X(idx,1)=(r(idx)-b(idx).*Y(idx,1))./a(idx);
    X(idx,2)=(r(idx)-b(idx).*Y(idx,2))./a(idx);
end

idx=~idx;
if sum(idx)>0
    Y(idx,1)=(r(idx)-a(idx).*X(idx,1))./b(idx);
    Y(idx,2)=(r(idx)-a(idx).*X(idx,2))./b(idx);
end

close all
figure('color','w')
plot(C1(:,1),C1(:,2),'-b','LineWidth',2), hold on, axis equal
plot(C2(:,1),C2(:,2),'-g','LineWidth',2), hold on, axis equal
for i=1:30:m
    plot(X(i,:),Y(i,:),'-k')
    plot(C1(i,1),C1(i,2),'.b','MarkerSize',20)
end
drawnow


