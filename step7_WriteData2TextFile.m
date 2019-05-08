function step7_WriteData2TextFile(Data,FileName)
% Write data extracted in step 6 to a text file.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: June.2014
%


[PD,~,Ext]=fileparts(FileName);
if isempty(Ext)
    FileName=strcat(FileName,'.txt');
end

fid=fopen(FileName,'w');
if fid<3
    error('Unable to open file: %s',FileName)
end

fprintf(fid,'%s\n%s\n%s\n\n','*** VIEW DOCUMENT CONTENTS USING WORDPAD ***',....
                             '*** ALL UNITS ARE IN PIXELS ***',...
                             '*** LAYER 1 STARTS AT PIA (see cortex_boundaries.tif for verification) ***' );

fprintf(fid,'%s ','Data folder:');
fprintf(fid,'%s\n\n',PD);

fprintf(fid,'%s ','Date & time of analysis:');
fprintf(fid,'%s\n\n',datestr(clock));

fprintf(fid,'%s\n\n','**************************************************');

fprintf(fid,'%s\n','White matter to pia (average distance):');
fprintf(fid,'\t%.5f\n\n',Data.dist_pia2wtm);

fprintf(fid,'%s\n','Pia to white matter (average distance):');
fprintf(fid,'\t%.5f\n\n',Data.dist_wtm2pia);

fprintf(fid,'%s\n','Average width of the cortex:');
fprintf(fid,'\t%.5f\n\n',Data.ave_width_cortex);

n=numel(Data.ave_width_layers);
str=sprintf('Average width of the cortex layers (1 to %u):',n);
fprintf(fid,'%s\n',str);
fprintf(fid,'\t%.5f\n',Data.ave_width_layers(:)');
fprintf(fid,'\n');

str=sprintf('Area of the cortex layers (1 to %u):',n);
fprintf(fid,'%s\n',str);
fprintf(fid,'\t%.5f\n',Data.area_layers(:)');

Fields={'DAPI' 'R0' 'G0' 'GR0' 'R1' 'G1' 'GR1'};
Title={'DAPI',...
       'ZNF312 (original)',...
       'CUX2 (original)', ...
       'CUX2 - ZNF312 (original)',...
       'ZNF312 (after intersection with DAPI)',...
       'CUX2 (after intersection with DAPI)', ...
       'CUX2 - ZNF312 (after intersection with DAPI)'};
       
   
for i=1:numel(Fields)
    
    if ~isfield(Data,Fields{i}), continue; end
    Data_i=getfield(Data,Fields{i}); %#ok<GFLD>
    n=numel(Data_i.CellCnt);
    
    fprintf(fid,'\n\n%s\n',Title{i});
    fprintf(fid,'%s\n','=================================================');

    fprintf(fid,'%s\n','# of cells per layer:');
    fprintf(fid,'\t%u\n',Data_i.CellCnt);
    
    
    fprintf(fid,'\n\n%s\n','Relative distance from pia (0 to 1 range):');
    fprintf(fid,'%s\n','-------------------------------------------------');
    for j=1:n
        if j==1
            fprintf(fid,'%s\n',sprintf('layer %u',j));
        else
            fprintf(fid,'\n\n%s\n',sprintf('layer %u',j));
        end
        fprintf(fid,'%.5f ',Data_i.RelativeDistance{j}');
    end

    
    fprintf(fid,'\n\n%s\n','Cell area (pixels^2):');
    fprintf(fid,'%s\n','-------------------------------------------------');
    for j=1:n
        if j==1
            fprintf(fid,'%s\n',sprintf('layer %u',j));
        else
            fprintf(fid,'\n\n%s\n',sprintf('layer %u',j));
        end
        fprintf(fid,'%u ',Data_i.CellArea{j}');
    end
    fprintf('.')
end

flag=fclose(fid);
if flag~=0
    error('Unable to write data to file %s',FileName)
end

