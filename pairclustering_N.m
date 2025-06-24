function [ taxo_cell_out,taxo_single_node ] = pairclustering_N( tablem,N )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

taxo_single_node=[];
%N=max(max(tablem(:,1:2)));

group_line=zeros(N,1);



cell_num=0;
meta_aa=size(tablem);
for i=1:meta_aa(1)


  if group_line(tablem(i,1))~=0 & group_line(tablem(i,2))==0
   group_line(tablem(i,2))= group_line(tablem(i,1));

  elseif group_line(tablem(i,1))==0 & group_line(tablem(i,2))~=0
   group_line(tablem(i,1))= group_line(tablem(i,2));
   
  elseif group_line(tablem(i,1))==0 & group_line(tablem(i,2))==0
   cell_num = cell_num+1;
   group_line(tablem(i,1))= cell_num;
   group_line(tablem(i,2))= cell_num;    

    elseif group_line(tablem(i,1))~=group_line(tablem(i,2)) 
     meta=find(group_line == max(group_line(tablem(i,1)),group_line(tablem(i,2))) );    
     group_line(meta)= min(group_line(tablem(i,1)),group_line(tablem(i,2)));

  end

end


%clear ooxx
% for i=1:length(tablem)
% ooxx(i,1)=group_line(tablem(i,1));
% ooxx(i,2)=group_line(tablem(i,2));
% ooxx(i,3)= group_line(tablem(i,1)) - group_line(tablem(i,2));
% end


GN=unique(group_line);

if GN(1)==0
taxo_single_node= find(group_line==0);        

   for i=2:length(GN)
   taxo_cell_out{i-1}=find(group_line==GN(i));     
   end

else    
    
   for i=1:length(GN)
   taxo_cell_out{i}=find(group_line==GN(i));     
   end
end


end