function [Data,im]=step6_GetSampleData(IM1,IM2,R,G,D)
% Extract data from a single histology slice.
%
% INPUT:
%   - IM1       : image containing the inner and outer boundaries of the 
%                 cortex.
%   - IM2       : image containing polylines partitioning cortex into 
%                 layers.
%   - R,G,D     : segmented images (of the same size) corresponding to
%                 the "red", "green" and DAPI stains, respectively.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: June.2014
%

% PRELIMINARY PRE-PROCESSING ----------------------------------------------

% Extract pia and white matter lines
bw=step1_GetCortexBoundaries(IM1);
fprintf('.')

% Extract the coordinates of the lines and produce a closed contour
[C,idx,BW]=step2_ConstructClosedContour(bw);
fprintf('.')

% Map the contour to a trapezoid with unit area
[F,P]=step3_ParameterizeContour(C,idx);
fprintf('.')

% Get the polylines 
X=step4_GetCortexLayers(IM2,C,idx);
[C_sub,A_sub,idx_sub]=SubregionContours(C,idx,X);
fprintf('.')

% Clean-up and label segmented images
Do=D;
[R,G,D]=step5_CleanUpSegmentedImages(C,idx,BW,R,G,D);
R=LabelCells(R);
G=LabelCells(G);
D=LabelCells(D);
clear BW bw
close all
fprintf('.')

% Visualize the extracted contour and polylines
close all
hf=figure('color','w','name','DO NOT MINIMIZE OR CLOSE! THIS FIGURE IS CURRENTLY BEING SAVED ...');
plot(C(:,1),C(:,2),'-k','LineWidth',2), hold on
axis equal off
h1=plot(C(idx(1):idx(2),1),C(idx(1):idx(2),2),'-g','LineWidth',3);
h2=plot(C(idx(3):idx(4),1),C(idx(3):idx(4),2),'-b','LineWidth',3);
for i=1:numel(X)
    plot(X{i}(:,1),X{i}(:,2),'--r','LineWidth',2)
end
XLim=get(gca,'XLim'); XLim=XLim+[-50 50];
YLim=get(gca,'YLim'); YLim=YLim+[-50 50];
set(gca,'XLim',XLim,'YLim',YLim,'Box','off','Xtick',[],'YTick',[],'YDir','reverse')
hl=legend([h1 h2],{'wtm' 'pia'});
warning('off')
maximize_fig(hf)
set(hl,'Location','SouthOutside','Orientation','horizontal','FontSize',40) 
im=export_fig('-r300','-a2',hf);
warning('on')
im=padarray(im,[10 10 0],255,'both');
close all
pause(0.5)

% WIDTH OF THE CORTEX & CORTEX LAYERS -------------------------------------

% Length of the wtm and pia polylines
C1=C(idx(1):idx(2),:); 
C2=C(idx(3):idx(4),:); 

L1=C1(2:end,:)-C1(1:end-1,:); L1=sum(sqrt(sum(L1.^2,2)));
L2=C2(2:end,:)-C2(1:end-1,:); L2=sum(sqrt(sum(L2.^2,2)));

% Distance from wtm line to pia and vice versa
W12=DistBetween2Contours(C1,C2);
W21=DistBetween2Contours(C2,C1);

% Average width of the cortex
W_ave=2*sum(A_sub)/sum(L1+L2);

% Average widths of the sub-regions
n=numel(C_sub);
W_sub=zeros(n,1);
for i=1:n
    
    % Lengths of the sub-region boundaries  
    L1=C_sub{i}(idx_sub(i,1):idx_sub(i,2),:);
    L2=C_sub{i}(idx_sub(i,3):idx_sub(i,4),:);
    
    L1=L1(2:end,:)-L1(1:end-1,:); L1=sum(sqrt(sum(L1.^2,2)));
    L2=L2(2:end,:)-L2(1:end-1,:); L2=sum(sqrt(sum(L2.^2,2)));
    
    % Average width
    W_sub(i)=2*A_sub(i)/(L1+L2);
    
end

Data.dist_wtm2pia=W12;
Data.dist_pia2wtm=W21;
Data.ave_width_cortex=W_ave;
Data.ave_width_layers=W_sub;
Data.area_layers=A_sub(:)';

% CELL DISTRIBUTION RELATED MEASUREMENTS ----------------------------------

% Measurements without DAPI overlap

% ZNF
[CellCnt_sub_R0,CellArea_sub_R0,RDP_R0]=CellMeasurements(R);
fprintf('.')

% CUX
[CellCnt_sub_G0,CellArea_sub_G0,RDP_G0]=CellMeasurements(G);
fprintf('.')

% DAPI
[CellCnt_sub_D,CellArea_sub_D,RDP_D]=CellMeasurements(D);
fprintf('.')

Data.DAPI.CellCnt=CellCnt_sub_D;
Data.DAPI.CellArea=CellArea_sub_D;
Data.DAPI.RelativeDistance=RDP_D;
clear CellCnt_sub_D CellArea_sub_D RDP_D

Data.R0.CellCnt=CellCnt_sub_R0;
Data.R0.CellArea=CellArea_sub_R0;
Data.R0.RelativeDistance=RDP_R0;
clear CellCnt_sub_R0 CellArea_sub_R0 RDP_R0

Data.G0.CellCnt=CellCnt_sub_G0;
Data.G0.CellArea=CellArea_sub_G0;
Data.G0.RelativeDistance=RDP_G0;
clear CellCnt_sub_G0 CellArea_sub_G0 RDP_G0

% Find 'green' - 'red' label image
bw=(G>0) & (R>0);
bw=imreconstruct(bw,G>0);

GR=G;
GR(bw)=0;
GR=RelabelImage(GR);
clear bw

% CUX-ZNF
[CellCnt_sub_GR0,CellArea_sub_GR0,RDP_GR0]=CellMeasurements(GR);
fprintf('.')

Data.GR0.CellCnt=CellCnt_sub_GR0;
Data.GR0.CellArea=CellArea_sub_GR0;
Data.GR0.RelativeDistance=RDP_GR0;
clear CellCnt_sub_GR0 CellArea_sub_GR0 RDP_GR0

% DAPI centroids (from the raw image)
Do=sum(Do,3);
Do=imfill(Do,'holes') & ~Do;
Do=bwareaopen(Do,6);
x=regionprops(Do,'Centroid');
DC=zeros(numel(x),2);
for i=1:numel(x), DC(i,:)=x(i).Centroid; end
DC=round(DC);
DC=sub2ind(size(R),DC(:,2),DC(:,1));
D=false(size(R));
D(DC)=true;
clear Do DC

% Remove objects in R and G that do not overlap with the centroids of 
% labelled objects in D
bw=imreconstruct(D&(R>0),R>0);
R(~bw)=0;
R=RelabelImage(R);

bw=imreconstruct(D&(G>0),G>0);
G(~bw)=0;
G=RelabelImage(G);
clear bw D

% Measurements with DAPI overlap
% ZNF & DAPI
[CellCnt_sub_R1,CellArea_sub_R1,RDP_R1]=CellMeasurements(R);
fprintf('.')

% CUX & DAPI
[CellCnt_sub_G1,CellArea_sub_G1,RDP_G1]=CellMeasurements(G);
fprintf('.')

Data.R1.CellCnt=CellCnt_sub_R1;
Data.R1.CellArea=CellArea_sub_R1;
Data.R1.RelativeDistance=RDP_R1;
clear CellCnt_sub_R1 CellArea_sub_R1 RDP_R1

Data.G1.CellCnt=CellCnt_sub_G1;
Data.G1.CellArea=CellArea_sub_G1;
Data.G1.RelativeDistance=RDP_G1;
clear CellCnt_sub_G1 CellArea_sub_G1 RDP_G1

% Find 'green' - 'red' label image
bw=(G>0) & (R>0);
bw=imreconstruct(bw,G>0);

GR=G;
GR(bw)=0;
GR=RelabelImage(GR);
clear bw

% (CUX & DAPI) - (ZNF & DAPI)
[CellCnt_sub_GR1,CellArea_sub_GR1,RDP_GR1]=CellMeasurements(GR);
fprintf('.')

Data.GR1.CellCnt=CellCnt_sub_GR1;
Data.GR1.CellArea=CellArea_sub_GR1;
Data.GR1.RelativeDistance=RDP_GR1;
clear CellCnt_sub_GR1 CellArea_sub_GR1 RDP_GR1
fprintf('.')


    function [CellCnt_sub,CellArea_sub,RDP]=CellMeasurements(L)
    % Measure cell size, cell density and relative distance of cells to 
    % pia in a labelled image L.
    %    
    
    % Cell centroids and sizes
    rp=regionprops(L,'Centroid','Area');
    clear L
    m=numel(rp);
    CC=zeros(m,2);
    CA=zeros(m,1);
    for k=1:m
        CC(k,:)=rp(k).Centroid;
        CA(k)=rp(k).Area;
    end
    clear rp
    
    % For different cortex layers find: list of cell sizes, number of
    % cells and relative distance of cells to pia
    m=numel(C_sub);
    CellCnt_sub=zeros(1,m);
    [CellArea_sub,RDP]=deal(cell(1,m));
    CHK=false(size(CC,1),1);
    for k=1:m
        
        % Centroids in the k-th layer
        chk=inpoly(CC,C_sub{k});
        CHK=CHK|chk;
        
        % Number of cells in the k-th layer and their sizes
        CellCnt_sub(k)=sum(chk);
        CellArea_sub{k}=CA(chk);
        
        % Normalized distance to pia
        [B,f]=Cartesian2Barycentric(F,C,CC(chk,:));
        f=Barycentric2Cartesian(F,P,B,f);
        RDP{k}=f(:,1);
        clear B f chk
        
    end
        
    % CHK is used in case some cell centroids fall on the boundaries between
    % different layers. While this is highly improbable, CHK is used as
    % a precaution to enforce the condition that the total number of cells
    % enclosed by the pia and white matter boundaries must equal the sum 
    % of the cells contained in individual sub-regions of the cortex. --> 
    % This condition may not hold if there are cells on the boundaries, so
    % these cells are removed. 
    
    end

end


function L=RelabelImage(Lo)
% Given a labelled image where some of the labels have been removed, 
% reassign the labels so that max(L(:)) equals the number of unique labels
% in Lo.


siz=size(Lo);

[L,srt_fwd]=sort(Lo(:));
clear Lo
[l,~,idx_unq]=unique(L);

l_new=(0:numel(l)-1)';
L=l_new(idx_unq);
clear idx_unq

srt_bck=(1:numel(L))';
srt_bck=srt_bck(srt_fwd);
[~,srt_bck]=sort(srt_bck,'ascend');
clear srt_fwd

L=L(srt_bck);
clear srt_bck
L=reshape(L,siz);

end

