function [ taxo_cell_sort, sorted_min_array_2,taxo_single_node ] = MSC_clustering( Dism, level,threshold_value )
%  ======= Usage ======== 
% 
%[ taxo_cell_sort, sorted_min_array_3 ] = MSC_clustering( Dism, level,threshold_value )
%
% @Input parameters 
% *Dism: The distance matrix 
%        ex:   0  2  3  1
%              2  0  5  3
%              3  5  0  2
%              1  3  2  0 
%
%  * level           : the level-th minimum distance are consideed here.
%  * threshold value : A  threshold value for sparification 
%
% @ output parameter
%  * taxo_cell_sort     : MSC taxonomy table, each cell contains a group of members. 
%  * sorted_min_array_3 : the pairwise table for clustering 
%
%  ======================

index_th=1;
if (nargin < 2)
level=1;    
end

if (nargin < 3)
threshold_value=2*max(max(Dism));    
index_th=0;
end



clear taxo_cell
N=max(size(Dism));

for i=1:N;Dism(i,i)=inf;end; 

if level ==1
sorted_min_array = find_min_of_matrix(Dism);
else
sorted_min_array = find_mins_of_matrix(Dism,level);    
end

   

if index_th ==1
 meta=size(sorted_min_array);
 index=1;
  for i=1:meta(1)
    if sorted_min_array(i,3) <threshold_value
    sorted_min_array_2(index,1)= sorted_min_array(i,1);
    sorted_min_array_2(index,2)= sorted_min_array(i,2);
    sorted_min_array_2(index,3)= sorted_min_array(i,3);
    sorted_min_array_3{index,1}= num2str(sorted_min_array(i,1));
    sorted_min_array_3{index,2}= num2str(sorted_min_array(i,2));
    sorted_min_array_3{index,3}= num2str(sorted_min_array(i,3));    
    index=index+1;
    end
  end
else
 sorted_min_array_2=sorted_min_array;    
end

if exist('sorted_min_array_2')==1   
 [taxo_cell,taxo_single_node]= pairclustering(sorted_min_array_2);


 clear ooxx
 for i=1:length(taxo_cell)
     ooxx(i)=length(taxo_cell{i});
 end

 [A,IXa]=sort(ooxx);

 for i=1:length(IXa)
     taxo_cell_sort{i}= taxo_cell{IXa(length(IXa) -i +1)};
 end
 
else
 'Note: The threshold value is too small, all nodes are single node in this run'   
 taxo_cell_sort=[];    
 taxo_single_node=1:length(Dism);
end



end




%% including functions


function sorted_min_array = find_min_of_matrix(dist_matrix)

m=length(dist_matrix);
min_array=[];
for i=1:m
    n=min(dist_matrix(i,:));
    minl=find(dist_matrix(i,:)==n);
    l=length(minl);
    for j=1:l
        min_array=[min_array;i,minl(j),dist_matrix(i,minl(j))];
    end
end
mm=length(min_array);

[a,index]=sort(min_array(:,3),'ascend');
for k=1:mm
    sorted_min_array(k,:)=min_array(index(k),:);
end

end

function sorted_min_array_out = find_mins_of_matrix(dist_matrix,level)

sorted_min_array_out =[];
max_v=max(max(dist_matrix));

for ii=1:level
m=length(dist_matrix);
min_array=[];

for i=1:m
    n=min(dist_matrix(i,:));
    minl=find(dist_matrix(i,:)==n);
    l=length(minl);
    for j=1:l
        min_array=[min_array;i,minl(j),dist_matrix(i,minl(j))];
        dist_matrix(i,minl(j)) = inf;
    end
end
mm=length(min_array);

[a,index]=sort(min_array(:,3),'ascend');
for k=1:mm
    sorted_min_array(k,:)=min_array(index(k),:);
end

sorted_min_array_out=[sorted_min_array_out',sorted_min_array']';


end


end

function [ taxo_cell_out,taxo_single_node ] = pairclustering( tablem )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

taxo_single_node=[];
N=max(max(tablem(:,1:2)));

group_line=zeros(N,1);



cell_num=0;
for i=1:length(tablem)


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