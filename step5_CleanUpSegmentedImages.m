function [R,G,D]=step5_CleanUpSegmentedImages(C,idx,BW,R,G,D,s)
% Clean-up segmented 'red', 'green' and DAPI stained images by removing 
% cells that fall outside the region enclosed by the pia and white matter
% boundaries. Additionally, remove cells that are within user-defined 
% distance away from pia.
%
% INPUT:
%   - C, idx, BW    : outputs of 'step2_ConstructClosedContour' function.
%   - R, G, D       : segmentations of 'red', 'green' and DAPI stained 
%                     images.
%   - s             : distance away from pia threshold. s=10 (pixels) is
%                     the default setting.
%
% OUTPUT:
%  - R, G, D        : processed R, G and D images in binary form.  
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: Mar.2014
%


if nargin<5, G=[]; end
if nargin<6, D=[]; end
if nargin<7 || isempty(s), s=10; end

% Convert segmented images to binary
RGD={sum(R,3)<50};
if ~isempty(G)
    RGD{2}=sum(G,3)<50;
    if ~isequal(size(RGD{1}),size(RGD{2}))
        error('Size of the ''red'' image does not match that of the ''green'' image')
    end
end
if ~isempty(D)
    RGD{3}=sum(D,3)<50;
    if ~isequal(size(RGD{1}),size(RGD{3}))
        error('Size of the ''red'' image does not match that of the DAPI image')
    end
end
clear R G D

% Remove cells that are outside the contour
BW=~imfill(BW,'holes');
M=numel(RGD);
for i=1:M
    bw=BW&RGD{i};
    bw=imreconstruct(bw,RGD{i});
    RGD{i}=RGD{i}&(~bw);
end

% Identify cells that are located approximately within r units of the 
% contour boundaries
se=strel('disk',ceil(s),0);
BW_r=imdilate(BW,se);
BW_r=(~BW)&BW_r;

RGD_r=RGD;
for i=1:M
    bw=BW_r&RGD{i};
    RGD_r{i}=imreconstruct(bw,RGD{i});
end
clear BW BW_r

% Remove cells that are within r units of the pia
C_pia=C(idx(3):idx(4),:);
for i=1:M
    
    % Get the foreground pixel co-ordinates
    L=bwlabel(RGD_r{i});
    P=regionprops(L,'PixelList'); %#ok<*MRPBW>
    n=max(L(:));
    if n==0, continue; end
    
    % Check the distance to pia
    idx=false(n,1);
    for j=1:n
        Dj=Pt2ContourDistance(C_pia,P(j).PixelList);
        if sum(Dj<=s), idx(j)=true; end
    end
    idx=find(idx);
    
    % Remove cells within r units away from pia
    if isempty(idx), continue; end
    L=ismember(L,idx);
    RGD{i}=RGD{i}&(~L);
    
end
clear RGD_r L P

R=RGD{1}; 
if M>1, G=RGD{2}; end 
if M>2, D=RGD{3}; end

