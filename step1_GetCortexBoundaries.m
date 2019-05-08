function bw=step1_GetCortexBoundaries(IM)
% Extract binary image containing two polylines defining the inner and 
% outer boundaries of the cortex.
%
% INPUT:
%   - IM    : RGB image where the boundaries of the cortex are represented
%             by BLACK polylines.
%
% OUPUT:
%   - bw    : binary image (same size as IM) containing two boundary 
%             segments.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: Feb.2014
%


% Convert to binary
bw=sum(double(IM),3)<50;

% Find largest foreground objects
bw=bwareaopen(bw,11);
L=bwlabel(bw);
m=max(L(:));
if m<2
    error('Unable to find cortex boundaries')
end

% Merge disconnected line segments until only two lines remain
bw=MergeLineSegments(bw);


function bw=MergeLineSegments(bw)

% Thin the lines to one pixel
bw=bwmorph(bw,'dilate');
bw=bwmorph(bw,'thin',Inf);

% Get the end-points
bw_end=bwmorph(bw,'endpoints');
idx=find(bw_end(:));
[y,x]=ind2sub(size(bw),idx);
n=numel(y);
if n==4, return; end

% Remove spurious endpoints (if any)
idx=false(n,1);
for i=1:n
    
    % Crop a small window around the end-point
    y_lim=[y(i)-10 y(i)+10];
    x_lim=[x(i)-10 x(i)+10];
    if y_lim(1)<1, y_lim(1)=1; end
    if y_lim(2)>size(bw,1), y_lim(2)=size(bw,1); end
    if x_lim(1)<1, x_lim(1)=1; end
    if x_lim(2)>size(bw,2), x_lim(2)=size(bw,2); end
    W=bw(y_lim(1):y_lim(2),x_lim(1):x_lim(2));
    m=sum(W(:));
    
    % Remove spurs
    W=bwmorph(W,'spur',5);
    if (m-sum(W(:)))<5
        bw(y_lim(1):y_lim(2),x_lim(1):x_lim(2))=W;
        idx(i)=true;
    end

end
y=y(:); y(idx)=[];
x=x(:); x(idx)=[];
n=numel(y);
if n==4, return; end

% Compute distances between the endpoints
X1=permute([y x],[1 3 2]);
X2=permute([y x],[3 1 2]);
D=bsxfun(@minus,X1,X2);
D=sqrt(sum(D.^2,3));
D(1:(n+1):end)=Inf;

% Connect lines segments until only two are left
while size(D,1)>4
    
    % Find a pair of closest end-points and connect them
    [D_min,idx]=min(D(:));
    [idx1,idx2]=ind2sub(size(D),idx);
    
    % Connect the endpoints by a straight line
    X1=[y(idx1) x(idx1)];
    X2=[y(idx2) x(idx2)];

    d=X2-X1; 
    t=linspace(0,D_min,3*D_min);    
    l=bsxfun(@times,d/D_min,t(:));
    l=bsxfun(@plus,X1,l);

    % Add the connecting line segment to the image 
    l=round(l); 
    l=sub2ind(size(bw),l(:,1),l(:,2));
    bw(l)=true;

    % Remove the connected end-points from the list
    y([idx1,idx2])=[];   x([idx1,idx2])=[];
    D([idx1,idx2],:)=[]; D(:,[idx1,idx2])=[];


end

if size(D,1)~=4
    error('Unable to extract cortex boundaries')
end

