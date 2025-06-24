function [ Network_pos ] = get_cyto_pos( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


clear Network_pos
index_o=0;
index_t=0;
index_s=1;
index_i=1;

IDname_meta=[];

 fid=fopen(filename);
 tic
  while index_s >0
      clear meta
  
      % time line
      index_t=index_t+1;
  
  
      
      % read a line in file for one times 
      ooxx=fgetl(fid);
      meta_aa=findstr(ooxx,'"');
      
      if length(findstr(ooxx,'<att name="shared name" value="')) ==1
          index_o=index_o+1;
          Network_pos{index_o,1}=ooxx(meta_aa(3)+1:meta_aa(4)-1);
      end

      if length(findstr(ooxx,' y="')) ==1
      if length(findstr(ooxx,' x="')) ==1
          ooxx;
         meta_strx=ooxx(findstr(ooxx,' x="'):min(findstr(ooxx,' x="')+50,length(ooxx)));
         meta_aax=findstr(meta_strx,'"');
          Network_pos{index_o,2}=str2num(meta_strx(meta_aax(1)+1:meta_aax(2)-1));

         meta_stry=ooxx(findstr(ooxx,' y="'):min(findstr(ooxx,' y="')+50,length(ooxx)));
         meta_aay=findstr(meta_stry,'"');          
          Network_pos{index_o,3}=str2num(meta_stry(meta_aay(1)+1:meta_aay(2)-1));
      end
      end
           
      
      
      % stop loop condition
      if length (ooxx) ==1
      if  ooxx == -1
          index_s =-1;
      end    
      end
      
  end
  
  toc
  % close fid
  fclose(fid);


Network_pos=Network_pos(find(cellfun('length',Network_pos(:,2))==1),:);
  
end

