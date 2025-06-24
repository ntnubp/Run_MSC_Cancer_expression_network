function [ min_Dis_matrix_dp,connection_table_out ] = MSC_tree( Dis, taxo_multiple_MSC,NE )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N=length(Dis);
M=length(taxo_multiple_MSC);
max_v=max(max(Dis));

%get minimum connetcion table of first level
min_Dis_matrix=zeros(N);
connection_table=zeros(N,5);
index=1;
for i=1:N
    meta_vv=Dis(i,:);
    meta_vv(i)=100*max_v;
    
    min_index=find(meta_vv==min(meta_vv));
    
    for ii=1:length(min_index)
        DD=Dis(i,min_index(ii));
    min_Dis_matrix(i,min_index(ii))=1;    
    if DD==0
        min_Dis_matrix(i,min_index(ii))=0.3;        
    end
    
    connection_table(index,1)=i;
    connection_table(index,2)=min_index(ii);
    connection_table(index,3)=DD;
    connection_table(index,4)=0;
    connection_table(index,5)=0;
    
    if connection_table(index,3)==0;
    min_Dis_matrix(i,min_index(ii))=0.3;
    end
    
    index=index+1;
    
    end
    
end

for i=1:M
  if length(taxo_multiple_MSC{i})>1
      
  meta_mm=1;
  for ii=1:length(NE{i})
    meta_vv=NE{i}(ii,:);
    meta_vv(ii)=100*max_v;   
    min_index=find(meta_vv==min(meta_vv));
    
    for jj=1:length(min_index)
    meta_dis=Dis(taxo_multiple_MSC{i}{ii},taxo_multiple_MSC{i}{min_index(jj)});
    
    min_gv=min(min(meta_dis));

    for qq=1:length(taxo_multiple_MSC{i}{ii})
        for kk=1:length(taxo_multiple_MSC{i}{min_index(jj)})
    DD=Dis(taxo_multiple_MSC{i}{ii}(qq),taxo_multiple_MSC{i}{min_index(jj)}(kk));
    if DD==min_gv
    min_Dis_matrix(taxo_multiple_MSC{i}{ii}(qq),taxo_multiple_MSC{i}{min_index(jj)}(kk))=meta_mm;

    
    
    connection_table(index,1)=taxo_multiple_MSC{i}{ii}(qq);
    connection_table(index,2)=taxo_multiple_MSC{i}{min_index(jj)}(kk);
    connection_table(index,3)=DD;
    connection_table(index,4)=0;
    connection_table(index,5)=i;

    
    index=index+1;
    end
        end
    end
    end
    
  end
  
    
  end
end



min_Dis_matrix_dp=min_Dis_matrix+ min_Dis_matrix' ;



  for i=1:length(connection_table)
  connection_table(i,4)= min_Dis_matrix_dp(connection_table(i,1),connection_table(i,2));   
  end


index=1;  
for i=1:length(connection_table)
    if connection_table(i,4)~=1
      if connection_table(i,1)>connection_table(i,2)
         connection_table_out(index,:)=connection_table(i,:);
         index=index+1; 
      end        
    else
        for k=1:5
        connection_table_out(index,k)=connection_table(i,k);
        end
    index=index+1;
    end
end
  

  
end

