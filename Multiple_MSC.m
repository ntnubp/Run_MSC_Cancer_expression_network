function [ taxo_multiple_MSC,NE] = Multiple_MSC( Dis, Th_v, GGindex )
% 
%
%  * 20240822 set GGindex into each level (by Hu)
%
%Thid function is used to slove multiple level MSC clustering results with different threshold
% value (up to 9th level MSC, it would be enouth if the size of distance matrix
% < 70000 ).
% 
% == input ==
%  Dis  : the input Distance matrix
%  th_v : the threshold value (up to seven level, if the length of th_v < 7, the over length threshold value  would be setting as the last value of Th_v.)
%  GGindex: setting the mode of group-group relation
%        GGindex=0            ; min  - the group-group distance is the minimum of all inter-cluster distance (default) 
%        GGindex=1 (or others); mean - the group-group distance is the average of all inter-cluster distance
% 
% == output ==
% 
% taxo_multiple_MSC : Cell of MSC clusteing results;
% 
%                     taxo_multiple_MSC{1} :  leve l clustering result
%                     taxo_multiple_MSC{2} :  leve 2 clustering result
%                     taxo_multiple_MSC{3} :  leve 3 clustering result
%                                        ...
% 
%           NE     : The renormalization matrix of each level. 
%                             NE{1}        :  level 1 renormalization matrix 
%                             NE{2}        :  level 2 renormalization matrix 
%                             NE{3}        :  level 3 renormalization matrix 
%                                        ...
%  Note: Sometimes the size of taxo_multiple_MSC{i} > length of NE{i}, it
%  is caused by the single nodes effect (single nodes did not considered into NE{i}).
% 
% =============                                        ...
%
%
%                                        last modified in 2016 0523 by Geng-Ming Hu
%                                                 
%

if (nargin < 2)
Th_v =max(max(Dis))*100;   
end

if (nargin < 3)
GGindex=0;
end

% setting Th_v
meta_ll=length(Th_v);
if meta_ll <9
    for i=meta_ll+1:9
        Th_v(i)=Th_v(meta_ll);
    end    
end

% setting GGindex
meta_ll=length(GGindex);
if meta_ll <9
    for i=meta_ll+1:9
        GGindex(i)=GGindex(meta_ll);
    end    
end

aa=201605031200;
LTP=5;


%~~ level 1
if aa >2


taxo_multiple_MSC{1}= MSC_clustering(Dis,1,Th_v(1));

    % Renormalization group
[NE_1m,NE_1a]=clustering_distance_mm2( taxo_multiple_MSC{1},Dis );
if GGindex(1) ==0
NE_1=NE_1m;
else
NE_1=NE_1a;    
end
NE{1}=NE_1;

aa=length(NE_1);

     end

% ~~

%~~ level 2
if aa > LTP/2 
LL=2;    
 clear taxo_meta
 taxo_meta= MSC_clustering(NE_1,1,Th_v(2));
     % Renormalization group
     if length(taxo_meta)==0
    NE{2}=[];
    aa=0;
     else     
     
[NE_2m,NE_2a]=clustering_distance_mm2( taxo_meta, NE_1 );
if GGindex(2) ==0
NE_2=NE_2m;
else
NE_2=NE_2a;    
end
NE{2}=NE_2;
 
   meta_cell=[];
 for i=1:length(taxo_meta)
  meta_cell=[meta_cell',taxo_meta{i}']';    
 end
 singlegroup=setxor(1:aa,meta_cell);  % aa is the length of last NGD
 meta_LL=length(taxo_meta);
for i=1:length(singlegroup)
taxo_meta{meta_LL+i}= singlegroup(i);   
end
 
 aa=length(NE_2); 


 
 
 for i=1:length(taxo_meta)
    taxo_multiple_MSC{2}{i}=[];
    for j=1:length(taxo_meta{i})
    taxo_multiple_MSC{2}{i}=[taxo_multiple_MSC{2}{i}',taxo_multiple_MSC{1}{taxo_meta{i}(j)}']';    
    end
 end
     end

end
%~~


%~~ level 3
if aa > LTP/2 
LL=2;    
 clear taxo_meta
 taxo_meta= MSC_clustering(NE_2,1,Th_v(3));
     % Renormalization group
     if length(taxo_meta)==0
    NE{3}=[];
    aa=0;
     else     
[NE_3m,NE_3a]=clustering_distance_mm2( taxo_meta, NE_2 );
if GGindex(3) ==0
NE_3=NE_3m;
else
NE_3=NE_3a;    
end
NE{3}=NE_3;

 
 
   meta_cell=[];
 for i=1:length(taxo_meta)
  meta_cell=[meta_cell',taxo_meta{i}']';    
 end
 singlegroup=setxor(1:aa,meta_cell);  % aa is the length of last NGD
 meta_LL=length(taxo_meta);
for i=1:length(singlegroup)
taxo_meta{meta_LL+i}= singlegroup(i);   
end
 
 aa=length(NE_3); 


 
 
 for i=1:length(taxo_meta)
    taxo_multiple_MSC{3}{i}=[];
    for j=1:length(taxo_meta{i})
    taxo_multiple_MSC{3}{i}=[taxo_multiple_MSC{3}{i}',taxo_multiple_MSC{2}{taxo_meta{i}(j)}']';    
    end
 end

     end
end
%~~



%~~ level 4
if aa > LTP/2 
LL=2;    
 clear taxo_meta

 taxo_meta= MSC_clustering(NE_3,1,Th_v(4));
     % Renormalization group
     if length(taxo_meta)==0
    NE{4}=[];
    aa=0;
     else

[NE_4m,NE_4a]=clustering_distance_mm2( taxo_meta, NE_3 );
if GGindex(4) ==0
NE_4=NE_4m;
else
NE_4=NE_4a;    
end
NE{4}=NE_4;

 
   meta_cell=[];
 for i=1:length(taxo_meta)
  meta_cell=[meta_cell',taxo_meta{i}']';    
 end
 singlegroup=setxor(1:aa,meta_cell);  % aa is the length of last NGD
 meta_LL=length(taxo_meta);
for i=1:length(singlegroup)
taxo_meta{meta_LL+i}= singlegroup(i);   
end
 
 aa=length(NE_4); 


 
 
 for i=1:length(taxo_meta)
    taxo_multiple_MSC{4}{i}=[];
    for j=1:length(taxo_meta{i})
    taxo_multiple_MSC{4}{i}=[taxo_multiple_MSC{4}{i}',taxo_multiple_MSC{3}{taxo_meta{i}(j)}']';    
    end
 end

     end
end
%~~


%~~ level 5
if aa > LTP/2 
LL=2;    
 clear taxo_meta
 taxo_meta= MSC_clustering(NE_4,1,Th_v(5));
% NE_5=clustering_distance_mm2( taxo_meta,NE_4 );
     % Renormalization group
     if length(taxo_meta)==0
    NE{5}=[];
    aa=0;
     else     
     
[NE_5m,NE_5a]=clustering_distance_mm2( taxo_meta, NE_4 );
if GGindex(5) ==0
NE_5=NE_5m;
else
NE_5=NE_5a;    
end
NE{5}=NE_5;
 

   meta_cell=[];
 for i=1:length(taxo_meta)
  meta_cell=[meta_cell',taxo_meta{i}']';    
 end
 singlegroup=setxor(1:aa,meta_cell);  % aa is the length of last NGD
 meta_LL=length(taxo_meta);
for i=1:length(singlegroup)
taxo_meta{meta_LL+i}= singlegroup(i);   
end
 
 aa=length(NE_5); 


 
 
 for i=1:length(taxo_meta)
    taxo_multiple_MSC{5}{i}=[];
    for j=1:length(taxo_meta{i})
    taxo_multiple_MSC{5}{i}=[taxo_multiple_MSC{5}{i}',taxo_multiple_MSC{4}{taxo_meta{i}(j)}']';    
    end
 end
     end

end
%~~


%~~ level 6
if aa > LTP/2 
LL=2;    
 clear taxo_meta
 taxo_meta= MSC_clustering(NE_5,1,Th_v(6));
% NE_6=clustering_distance_mm2( taxo_meta,NE_5 );
      % Renormalization group
     if length(taxo_meta)==0
    NE{6}=[];
    aa=0;
     else      
      
[NE_6m,NE_6a]=clustering_distance_mm2( taxo_meta, NE_5 );
if GGindex(6) ==0
NE_6=NE_6m;
else
NE_6=NE_6a;    
end
NE{6}=NE_6;

   meta_cell=[];
 for i=1:length(taxo_meta)
  meta_cell=[meta_cell',taxo_meta{i}']';    
 end
 singlegroup=setxor(1:aa,meta_cell);  % aa is the length of last NGD
 meta_LL=length(taxo_meta);
for i=1:length(singlegroup)
taxo_meta{meta_LL+i}= singlegroup(i);   
end
 
 aa=length(NE_6); 


 
 
 for i=1:length(taxo_meta)
    taxo_multiple_MSC{6}{i}=[];
    for j=1:length(taxo_meta{i})
    taxo_multiple_MSC{6}{i}=[taxo_multiple_MSC{6}{i}',taxo_multiple_MSC{5}{taxo_meta{i}(j)}']';    
    end
 end

     end
end
%~~


%~~ level 7
if aa > LTP/2 
LL=2;    
 clear taxo_meta
 taxo_meta= MSC_clustering(NE_6,1,Th_v(7));
% NE_7=clustering_distance_mm2( taxo_meta,NE_6 );
     % Renormalization group
     if length(taxo_meta)==0
    NE{7}=[];
    aa=0;
     else     
[NE_7m,NE_7a]=clustering_distance_mm2( taxo_meta, NE_6 );
if GGindex(7) ==0
NE_7=NE_7m;
else
NE_7=NE_7a;    
end
NE{7}=NE_7;
 
   meta_cell=[];
 for i=1:length(taxo_meta)
  meta_cell=[meta_cell',taxo_meta{i}']';    
 end
 singlegroup=setxor(1:aa,meta_cell);  % aa is the length of last NGD
 meta_LL=length(taxo_meta);
for i=1:length(singlegroup)
taxo_meta{meta_LL+i}= singlegroup(i);   
end
 
 aa=length(NE_7); 


 
 
 for i=1:length(taxo_meta)
    taxo_multiple_MSC{7}{i}=[];
    for j=1:length(taxo_meta{i})
    taxo_multiple_MSC{7}{i}=[taxo_multiple_MSC{7}{i}',taxo_multiple_MSC{6}{taxo_meta{i}(j)}']';    
    end
 end

     end
end
%~~

%~~ level 8
if aa > LTP/2 
LL=2;    
 clear taxo_meta
 taxo_meta= MSC_clustering(NE_7,1,Th_v(8));
% NE_8=clustering_distance_mm2( taxo_meta,NE_7 );
     % Renormalization group
     if length(taxo_meta)==0
    NE{8}=[];
    aa=0;
     else     
[NE_8m,NE_8a]=clustering_distance_mm2( taxo_meta, NE_7 );
if GGindex(8) ==0
NE_8=NE_8m;
else
NE_8=NE_8a;    
end
NE{8}=NE_8;
 
   meta_cell=[];
 for i=1:length(taxo_meta)
  meta_cell=[meta_cell',taxo_meta{i}']';    
 end
 singlegroup=setxor(1:aa,meta_cell);  % aa is the length of last NGD
 meta_LL=length(taxo_meta);
for i=1:length(singlegroup)
taxo_meta{meta_LL+i}= singlegroup(i);   
end
 
 aa=length(NE_8); 


 
 
 for i=1:length(taxo_meta)
    taxo_multiple_MSC{8}{i}=[];
    for j=1:length(taxo_meta{i})
    taxo_multiple_MSC{8}{i}=[taxo_multiple_MSC{8}{i}',taxo_multiple_MSC{7}{taxo_meta{i}(j)}']';    
    end
 end

     end
end
%~~

%~~ level 9
if aa > LTP/2 
LL=2;    
 clear taxo_meta
 taxo_meta= MSC_clustering(NE_8,1,Th_v(9));
% NE_9=clustering_distance_mm2( taxo_meta,NE_8 );
     % Renormalization group
     if length(taxo_meta)==0
    NE{9}=[];
    aa=0;
     else     
[NE_9m,NE_9a]=clustering_distance_mm2( taxo_meta, NE_8 );
if GGindex(9) ==0
NE_9=NE_9m;
else
NE_9=NE_9a;    
end
NE{9}=NE_9;
 
   meta_cell=[];
 for i=1:length(taxo_meta)
  meta_cell=[meta_cell',taxo_meta{i}']';    
 end
 singlegroup=setxor(1:aa,meta_cell);  % aa is the length of last NGD
 meta_LL=length(taxo_meta);
for i=1:length(singlegroup)
taxo_meta{meta_LL+i}= singlegroup(i);   
end
 
 aa=length(NE_9); 


 
 
 for i=1:length(taxo_meta)
    taxo_multiple_MSC{9}{i}=[];
    for j=1:length(taxo_meta{i})
    taxo_multiple_MSC{9}{i}=[taxo_multiple_MSC{9}{i}',taxo_multiple_MSC{8}{taxo_meta{i}(j)}']';    
    end
 end

     end
end
%~~


end


%===== including function ======


function [ clm_min,clm_mean ] = clustering_distance_mm2( taxo_ref,Dism )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for ii=1:length(taxo_ref)
for jj=1:length(taxo_ref)
    clear meta_m;
    for i= 1:length(taxo_ref{ii})
    for j= 1:length(taxo_ref{jj})
        meta_m(i,j)=Dism(taxo_ref{ii}(i),taxo_ref{jj}(j));
    end        
    end

clm_mean(ii,jj)=mean(mean(meta_m));    
clm_min(ii,jj)=min(min(meta_m));

end
end


end


function [ taxo_cell_sort, sorted_min_array_3 ] = MSC_clustering( Dism, level,threshold_value )
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
  
sorted_min_array_2=[];
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

if length(sorted_min_array_2) ==0
    taxo_cell_sort={};
    taxo_single_node=1:N;
else
[taxo_cell,taxo_single_node]= pairclustering_N(sorted_min_array_2,N);


clear ooxx
for i=1:length(taxo_cell)
    ooxx(i)=length(taxo_cell{i});
end

[A,IXa]=sort(ooxx);

for i=1:length(IXa)
    taxo_cell_sort{i}= taxo_cell{IXa(length(IXa) -i +1)};
end

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

function [ taxo_cell_out,taxo_single_node ] = pairclustering_N( tablem,NN )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

taxo_single_node=[];
N=NN;

group_line=zeros(N,1);


meta=size(tablem);
cell_num=0;
for i=1:meta(1)


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