function [ cyto_out ] = cytoscape_transform_MSCtree_mixsimDis_bkc( link_table,label_table,savename,node_color_list,level_index,pos )


switch_color=1;

index_save=1;
if (nargin < 3)
index_save=0;    
end

if (nargin < 3)
savefolder='';    
end

if (nargin < 4)
switch_color=0;
end

if (nargin < 5)
level_index=100000;
end

pos_index=1;
if (nargin < 6)
pos_index=0;
end


link_table=link_table(find(cell2mat(link_table(:,5))<=level_index),:);


if length(link_table) > 0
IDlist=unique([label_table(:,1)',link_table(:,1)',link_table(:,2)']');
else
IDlist=unique(label_table(:,1));
end

N=length(IDlist);
meta=size(link_table);


% ¨º­ÓÃC¦â¬O¶Â¬õ¶ÀºñÂÅÀQµµ from BK
color_list_level{1}='#FF0000';
color_list_level{2}='#FFFF00';
color_list_level{3}='#00FF00';
color_list_level{4}='#0000FF';
color_list_level{5}='#00FFFF';
color_list_level{6}='#FF00FF';



cyto_out{1}= '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>';
cyto_out{2}= '<graph id="8" label="64" directed="1" cy:documentVersion="3.0" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">';

% get nodes information
r=100;
node_color_in  ='#ccccff';
node_color_out ='#000000';
node_type ='ELLIPSE';
node_type_ideG='Octagon';

NL=9;    
[aa bb]=size(IDlist);    
for i=1:aa
    
    if switch_color ==1
        node_color_in  = node_color_list{find(strcmp(node_color_list(:,1),IDlist{i})==1),2};
        %node_color_out ='#000000';
        node_color_out =node_color_in ;
        
    end    
    
    if pos_index==0
    px=num2str(r*sin((i-1)*(2*pi)/N));
    py=num2str(r*cos((i-1)*(2*pi)/N));
    else
    meta_ii=strcmp(pos(:,1),IDlist{i});
    px=num2str(pos{meta_ii,2});
    py=num2str(pos{meta_ii,3});   
    end

    cyto_out{2+NL*(i-1)+1}= ['<node id="',num2str(i),'" label="',IDlist{i},'">'];
    cyto_out{2+NL*(i-1)+2}= ['  <att name="shared name" value="',IDlist{i},'" type="string"/>'];
    cyto_out{2+NL*(i-1)+3}= ['  <att name="selected" value="0" type="boolean"/>'];
    cyto_out{2+NL*(i-1)+4}= ['  <att name="name" value="',label_table{find(strcmp(label_table(:,1),IDlist{i})==1),2},'" type="string"/>'];  
        if  strcmp(IDlist{i}(length(IDlist{i})),')') ==0 
    cyto_out{2+NL*(i-1)+5}= ['<graphics y="',py,'" fill="',node_color_in,'" x="',px,'" z="0.0" h="70.0" type="',node_type,'" width="3.0" outline="',node_color_out,'" w="70.0">'];
        else
    cyto_out{2+NL*(i-1)+5}= ['<graphics y="',py,'" fill="',node_color_in,'" x="',px,'" z="0.0" h="90.0" type="',node_type_ideG,'" width="3.0" outline="',node_color_out,'" w="90.0">'];            
        end
    cyto_out{2+NL*(i-1)+6}= ['<att name="NODE_LABEL_COLOR" value="#333333" type="string"/>'];
    cyto_out{2+NL*(i-1)+7}= ['<att name="LABEL_FONT_SIZE" value="36" type="string"/>'];
    
%    cyto_out{2+7*(i-1)+5}= ['<graphics fill="',node_color_in,'" w="40.0" width="0.0" h="40.0" z="0.0" y="',py,'" outline="',node_color_out,'" type="',node_type,'" x="',px,'">']       
%    cyto_out{2+8*(i-1)+5}= ['<graphics outline=""$color"" type="ELLIPSE" y=""$yp"" h="5.0" width="0.5" w="5.0" x=""$xp"" fill=""$color"" z="0.0">']   
    cyto_out{2+NL*(i-1)+8}= ['</graphics>'];  
%    cyto_out{2+8*(i-1)+7}= ['  \n'];   
    cyto_out{2+NL*(i-1)+9}= ['</node>'];   
end

index_L1=length(cyto_out);

% get edges information
edge_type_s ='SOLID';
edge_type_di ='SOLID';
edge_type_zero ='SOLID';
edge_type_out ='SINEWAVE';%'SOLID';


arrow_type_s_s ='NONE';
arrow_type_s_di ='Arrow';
arrow_type_s_zero ='NONE';
arrow_type_s_out ='NONE';


arrow_type_t_s ='Arrow';
arrow_type_t_di ='Arrow';
arrow_type_t_zero ='NONE';
arrow_type_t_out ='NONE';


edge_color_s ='#000000';
edge_color_di ='#000000';
edge_color_zero ='#999933';
edge_color_out ='#0000FF';


edge_width_s='6';
edge_width_di='6';
edge_width_zero='6';
edge_width_out='6';

CL= 16;

if length(link_table) >0

[aa bb]=size(link_table);    
for i=1:aa
    index_eee=1;
    %1
    sc_str=sprintf('%1.1e',link_table{i,3});
    if strcmp(sc_str(6),'0')==1
        sc_str=sc_str([1:5,7]);
    end
    cyto_out{index_L1+CL*(i-1)+index_eee}  = ['<edge id="" label=" ',link_table{i,1},' (',sc_str,') ',link_table{i,2},'" source="',num2str(find(strcmp(IDlist,link_table{i,1})==1)),'" target="',num2str(find(strcmp(IDlist,link_table{i,2})==1)),'" cy:directed="1">'];
    index_eee=index_eee+1;
    %2
    cyto_out{index_L1+CL*(i-1)+index_eee}  = ['    <att name="interaction" value="',sc_str,'" type="string"/>'];
    index_eee=index_eee+1;
    %3
    cyto_out{index_L1+CL*(i-1)+index_eee}  = ['    <att name="shared name" value="',link_table{i,1},' (',sc_str,') ',link_table{i,2},'" type="string"/>'];
    index_eee=index_eee+1;
    %4
    cyto_out{index_L1+CL*(i-1)+index_eee}  = ['    <att name="selected" value="0" type="boolean"/>'];
    index_eee=index_eee+1;
%    cyto_out{index_L1+9*(i-1)+5}  = ['    <att name="name" value="',link_table{i,1},' (',num2str(i),') ',link_table{i,2},'" type="string"/>'];
    %5
    cyto_out{index_L1+CL*(i-1)+index_eee}  = ['    <att name="name" value="',link_table{i,1},' (',') ',link_table{i,2},'" type="string"/>'];
    index_eee=index_eee+1;
    %6
    cyto_out{index_L1+CL*(i-1)+index_eee}  = ['    <att name="shared interaction" value="',num2str(i),'" type="string"/>'];
    index_eee=index_eee+1;


    
if link_table{i,4} == 0.6

    edge_type    = edge_type_zero;
    edge_color   = edge_color_zero;
    edge_width   = edge_width_zero;
    arrow_type_s = arrow_type_s_zero;
    arrow_type_t = arrow_type_t_zero;

end    

    % high level line

edge_type  ='SOLID';

if link_table{i,5} ==0
edge_width = '6';   
edge_color_s ='#000000';
edge_color_di ='#000000';
end

if link_table{i,5} ==1
edge_width = '6';   
edge_color_s =color_list_level{link_table{i,5}};
edge_color_di =color_list_level{link_table{i,5}};
end

if link_table{i,5} ==2
edge_width = '6';   
edge_type  ='SOLID';
edge_color_s =color_list_level{link_table{i,5}};
edge_color_di =color_list_level{link_table{i,5}};
end

if link_table{i,5} ==3
edge_width = '6';   
edge_type  ='SOLID';
edge_color_s =color_list_level{link_table{i,5}};
edge_color_di =color_list_level{link_table{i,5}};
end

if link_table{i,5} ==4
edge_width = '6';   
edge_type  ='SOLID';
edge_color_s =color_list_level{link_table{i,5}};
edge_color_di =color_list_level{link_table{i,5}};
end

if link_table{i,5} ==5
edge_width = '6';   
edge_type  ='SOLID';
edge_color_s =color_list_level{link_table{i,5}};
edge_color_di =color_list_level{link_table{i,5}};
end

if link_table{i,5} ==6
edge_width = '6';   
edge_type  ='SOLID';
edge_color_s =color_list_level{link_table{i,5}};
edge_color_di =color_list_level{link_table{i,5}};
end

    %



  


if link_table{i,4} ==1.0
%     edge_type  =edge_type_s;
    edge_color = edge_color_s; %edge_color_list{linktable{i,4}}; %edge_color_list;
    edge_width = edge_width_s;
    arrow_type_s = arrow_type_s_s;
    arrow_type_t = arrow_type_t_s;
    
end    

if link_table{i,4} ==2.0
%     edge_type  = edge_type_di;
    edge_color = edge_color_di;
    edge_width = edge_width_di;
    arrow_type_s = arrow_type_s_di;
    arrow_type_t = arrow_type_t_di;    
end   

if link_table{i,4} ==2.5
%     edge_type  = edge_type_di;
    edge_color = edge_color_zero;
    edge_width = edge_width_di;
    arrow_type_s = arrow_type_s_di;
    arrow_type_t = arrow_type_t_di;
end   



if link_table{i,5} ~=0
 %edge_type  ='LONG_DASH';    
end
if link_table{i,5} >5
edge_type  ='SOLID';   
end

if link_table{i,4} ==-1
        % if link_table{i,5} >2 % for last level
    edge_type  = edge_type_out;
    edge_color = edge_color_out;
    edge_width = edge_width_out;
    arrow_type_s = arrow_type_s_out;
    arrow_type_t = arrow_type_t_out;
        % end
end



     %7
     cyto_out{index_L1+CL*(i-1)+index_eee}  = ['    <graphics width="',edge_width,'" fill="',edge_color,'">'];    
     index_eee=index_eee+1;
     %8
     cyto_out{index_L1+CL*(i-1)+index_eee}  = ['      <att name="EDGE_LINE_TYPE" value="',edge_type,'" type="string"/>'];
      index_eee=index_eee+1;
      %9
     cyto_out{index_L1+CL*(i-1)+index_eee}  = ['      <att name="EDGE_SOURCE_ARROW_SHAPE" value="',arrow_type_s,'" type="string"/>'];
      index_eee=index_eee+1;
      %10
     cyto_out{index_L1+CL*(i-1)+index_eee}  = ['      <att name="EDGE_TARGET_ARROW_SHAPE" value="',arrow_type_t,'" type="string"/>'];
      index_eee=index_eee+1;
     %11
     cyto_out{index_L1+CL*(i-1)+index_eee}  = ['      <att name="EDGE_TARGET_ARROW_SELECTED_PAINT" value="#333333" type="string"/>'];
       index_eee=index_eee+1;
       %12
     cyto_out{index_L1+CL*(i-1)+index_eee}  = ['      <att name="EDGE_SOURCE_ARROW_SELECTED_PAINT" value="#333333" type="string"/>'];    
      index_eee=index_eee+1;
      %13
     cyto_out{index_L1+CL*(i-1)+index_eee}  = ['      <att name="EDGE_TARGET_ARROW_UNSELECTED_PAINT" value="',edge_color,'" type="string"/>'];
       index_eee=index_eee+1;
       %14
     cyto_out{index_L1+CL*(i-1)+index_eee}  = ['      <att name="EDGE_SOURCE_ARROW_UNSELECTED_PAINT" value="',edge_color,'" type="string"/>'];    
      index_eee=index_eee+1;


    %15  
    cyto_out{index_L1+CL*(i-1)+index_eee}  = ['    </graphics>'];        
    index_eee=index_eee+1;
    %16
    cyto_out{index_L1+CL*(i-1)+index_eee}  = ['  </edge>'];
    index_eee=index_eee+1;    
end
end

index_L2= length(cyto_out);

cyto_out{index_L2+1}='</graph>';

% save file
 if index_save==1
  fid=fopen([savename,'.xml'],'w+');
  for i=1:length(cyto_out)
  fprintf(fid,'%s\n',cyto_out{i});
  end
  fclose(fid);
 end


end

