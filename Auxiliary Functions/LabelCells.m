function L=LabelCells(bw)
% Given a binary image containing cell outlines, generate a label image,
% where pixels with the same value belong to the same cell. While
% this step may seem trivial, in reality special care must be taken as 
% some of the boundaries of the segmented cells are either touching or
% overlapping. 
%
%   - bw    : one of binary images output by 'step5_CleanUpSegmentedImages'
%             function.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: Jun.2014
%

% Find 8-connected objects contained inside the segmentation contours
bw_prm=bw>0;
bw_fil=imfill(bw_prm,'holes');
bw_in=bw_fil & ~bw_prm;

bw_prm=imdilate(bw_in,true(3)) & bw_prm;
bw_fil=imfill(bw_prm,'holes');
clear bw_prm bw


% If the total number of 8-connected objects in the original image equals
% the number in bw_in, then were are done
L=bwlabel(bw_fil);
L_in=bwlabel(bw_in);

N=max(L(:));
N_in=max(L_in(:));
if N==N_in, return; end

% Overlap objects in L with L_in to find those that contain more than one 
% label (not including background)
P=regionprops(L,'PixelIdxList'); %#ok<*MRPBW>
idx=false(N,1); idx_in=false(N_in,1);
idx_del=false(N,1);
for i=1:N
    f=L_in(P(i).PixelIdxList);
    l=unique(f);
    if l(1)==0, l(1)=[]; end
    m=numel(l);
    if m>1
        if m==2 && (sum(f==l(1))<5 || sum(f==l(2))<5)
           continue
        end
        if m==3
            a=[sum(f==l(1)) sum(f==l(2)) sum(f==l(3))]<5;
            if sum(a)>1, continue; end
        end
        if m>5 % objects with m>5 are considered outliers
            idx_del(i)=true;
            continue
        end        
        idx(i)=true;
        idx_in(l)=true;
    end
end
clear P
if sum(idx)==0, return; end

if sum(idx_del)>0
    bw_del=ismember(L,find(idx_del));  
else
    bw_del=[];
end

% Remove outliers (if any)
bw_in=ismember(L,find(idx)) & bw_in;
if ~isempty(bw_del)
    bw_fil(bw_del)=0;
    bw_in(bw_del)=0;
end
clear bw_del
L_in=bwlabel(bw_in);
L=bwlabel(bw_fil);


% "Split-up" the label of touching objects
idx=unique(L(bw_in));
P=regionprops(L,{'PixelIdxList','BoundingBox'}); P=P(idx);
N=max(L(:));
for i=1:numel(idx)
    
    f=L_in(P(i).PixelIdxList);
    l=unique(f);
    if l(1)==0, l(1)=[]; end
    
    n=zeros(size(l));
    for j=1:numel(l), n(j)=sum(f==l(j)); end
    [~,srt]=sort(n,'descend');
    l=l(srt);
    
    uc=ceil(P(i).BoundingBox(1:2));
    dx=P(i).BoundingBox(3:4);
    
    x1=uc(1)-2;
    x2=uc(1)+dx(1)+1;

    y1=uc(2)-2;
    y2=uc(2)+dx(2)+1;
    
    if x1<1, x1=1; end
    if x2>size(L,2), x2=size(L,2); end
    
    if y1<1, y1=1; end
    if y2>size(L,1), y2=size(L,1); end
    
    Li=L(y1:y2,x1:x2);
    bw=Li==idx(i);
    L_in_i=L_in(y1:y2,x1:x2);
    for j=1:numel(l)
        bw_j=imdilate(L_in_i==l(j),true(3)) & bw;
        if j==1
            Li(bw_j)=idx(i);
        else
            N=N+1; 
            Li(bw_j)=N;
        end
        bw=bw&~bw_j;
    end
    Li(bw)=0;
    L(y1:y2,x1:x2)=Li;
end
clear L_in bw_in bw_fil P

% Find small objects and append them to large ones if they are adjacent,
% and delete them otherwise.
S=regionprops(L,'Area');
A=zeros(N,1);
for i=1:N, A(i)=S(i).Area; end

a=sort(A);
idx=find(A<a(floor(0.005*N)));
if isempty(idx), return; end

idx=sort(idx,'descend');
P=regionprops(L,{'BoundingBox'}); P=P(idx);
for i=1:numel(idx)
    
    uc=ceil(P(i).BoundingBox(1:2));
    dx=P(i).BoundingBox(3:4);
    
    x1=uc(1)-2;
    x2=uc(1)+dx(1)+1;

    y1=uc(2)-2;
    y2=uc(2)+dx(2)+1;
    
    if x1<1, x1=1; end
    if x2>size(L,2), x2=size(L,2); end
    
    if y1<1, y1=1; end
    if y2>size(L,1), y2=size(L,1); end
    
    Li=L(y1:y2,x1:x2);
    bw=imdilate(Li==idx(i),true(3));
    f=Li(bw);
    l=unique(f);
    l(l==0)=[];
    l(l==idx(i))=[];
    
    if ~isempty(l)
        chk=ismember(l,idx);
        l(chk)=[];
    end
    
    % If on its own and less than 25 pixels, delete it; keep otherwise.    
    if isempty(l)
        if sum(f==idx(i))<25
            Li(Li==idx(i))=0;
            L(y1:y2,x1:x2)=Li;
            L(L>idx(i))=L(L>idx(i))-1;
        end
        continue
    end
    
    if numel(l)==1 % one adjecent object
        Li(Li==idx(i))=l;
        L(y1:y2,x1:x2)=Li;
        L(L>idx(i))=L(L>idx(i))-1;
        continue
    end
        
    % Assign to the largest adjacent object
    n=zeros(size(l));
    for j=1:numel(l)
        n(j)=sum(L(:)==l(j));
    end
    [~,n]=max(n);
    Li(Li==idx(i))=l(n);
    L(y1:y2,x1:x2)=Li;
    L(L>idx(i))=L(L>idx(i))-1;
end



