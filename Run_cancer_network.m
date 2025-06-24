% ===============================================================
% Contents of this Script:
%----------------------------------------------------------------
%   1. Data Loading
%   2. Key Gene Identification
%   3. Distance matrix calcuation
%   4. MSClustering (MSC)
%   5. MSC Tree Construction and Cytoscape File Generation
% ===============================================================

clear all
% 1. Data Loading
    % Load mRNA expression data
    Table_in = readtable('mmR_data_mRNA9935.csv', 'ReadVariableNames', true, 'ReadRowNames', true,'VariableNamingRule', 'preserve');
    sample_list_mmR_mRNA = Table_in.Properties.VariableNames;
    gene_list_mmR_mRNA = Table_in.Properties.RowNames;
    data_matrix_mmR_mRNA = table2array(Table_in);
    clear Table_in
        

    % Load miRNA expression data
    Table_in = readtable('mmR_data_miRNA.csv', 'ReadVariableNames', true, 'ReadRowNames', true,'VariableNamingRule', 'preserve');
    sample_list_mmR_miRNA = Table_in.Properties.VariableNames;
    gene_list_mmR_miRNA = Table_in.Properties.RowNames;
    data_matrix_mmR_miRNA = table2array(Table_in);
    clear Table_in
        

    % Load RPPA (Reverse Phase Protein Array) data
    Table_in = readtable('mmR_data_RPPA.csv', 'ReadVariableNames', true, 'ReadRowNames', true,'VariableNamingRule', 'preserve');
    sample_list_mmR_RPPA = Table_in.Properties.VariableNames;
    gene_list_mmR_RPPA = Table_in.Properties.RowNames;
    data_matrix_mmR_RPPA = table2array(Table_in);
    clear Table_in
    

    % Load sample metadata
    T = readtable('ID_info_mmR_2439_o.xlsx', 'VariableNamingRule', 'preserve');
    ID_info_index_list_mmR_2439=T.Properties.VariableNames;
    ID_info_mmR_2439=table2cell(T);

    % Extract disease list and related information
    disease_list_mmR=unique( ID_info_mmR_2439(:,3));

    for i=1:length(disease_list_mmR)
        meta_aa=find(strcmp(ID_info_mmR_2439(:,3),disease_list_mmR{i})==1);
    disease_info_mmR{i,1}=disease_list_mmR{i};
    disease_info_mmR{i,2}=length(meta_aa);
    disease_info_mmR{i,3}=meta_aa;
    end    



% 2. Key Gene Identification
    % Identify key genes with Ht < 0.4
Ht=0.4;

clear gene_std_dd_mmR_mRNA
for i=1:length(disease_list_mmR)
meta_aa=find(strcmp(ID_info_mmR_2439(:,3),disease_info_mmR{i})==1);
gene_std_dd_mmR_mRNA{i,1}=std(data_matrix_mmR_mRNA');
gene_std_dd_mmR_mRNA{i,2}=std(data_matrix_mmR_mRNA(:,meta_aa)');
gene_std_dd_mmR_mRNA{i,3}=gene_std_dd_mmR_mRNA{i,2}./gene_std_dd_mmR_mRNA{i,1};
end

clear key_gene_index
key_gene_index_all=[];
for i=1:length(disease_list_mmR)
    key_gene_index{i}=find(gene_std_dd_mmR_mRNA{i,3}<Ht);
    key_gene_index_all=[key_gene_index_all,key_gene_index{i}];
end
key_gene_index_all=unique(key_gene_index_all);

data_matrix_mmR_mRNA_kgr40=data_matrix_mmR_mRNA(key_gene_index_all,:);



% 3. Distance matrix calcuation
    % key mRNA
Sim_cor_mmR_mRNA_kgr40=get_correlation_matrix(data_matrix_mmR_mRNA_kgr40);
Sim_dis_mmR_mRNA_kgr40_tr=(2./(1+Sim_cor_mmR_mRNA_kgr40))-1;

    % miRNA in log10 scale
Sim_cor_mmR_miRNA_log=get_correlation_matrix(log10(data_matrix_mmR_miRNA+10^(-5)));
Sim_dis_mmR_miRNA_log_tr=(2./(1+Sim_cor_mmR_miRNA_log))-1;

    % RPPA
Sim_cor_mmR_RPPA=get_correlation_matrix(data_matrix_mmR_RPPA);
Sim_dis_mmR_RPPA_tr=(2./(1+Sim_cor_mmR_RPPA))-1;

    % Consensus distance (key mRMA+ miRNA(log) + RPPA)
Sim_cor_mmR_mix2_kgr40=(Sim_cor_mmR_mRNA_kgr40+Sim_cor_mmR_miRNA_log+Sim_cor_mmR_RPPA)/3;
Sim_dis_mmR_mix2_kgr40_tr= (2./(1+Sim_cor_mmR_mix2_kgr40))-1;



%   4. MSClustering (MSC)
    % Get multiple leve MSC result using the consensus distance
[taxo_cell_multiple_mmR_mix2_kgr40_tr,NE_mmR_mix2_kgr40_tr]=Multiple_MSC(Sim_dis_mmR_mix2_kgr40_tr);

    % Extract connection table at MSC Level 2
[taxo_meta sort_min_array_MSC2]=MSC_clustering(NE_mmR_mix2_kgr40_tr{2});

    % Compute average distance matrix between MSC Level 2 clusters
clear NE_mmR_mix2_kgr40_tr_MSC1mMSC2a
for i=1:length(taxo_cell_multiple_mmR_mix2_kgr40_tr{2})
    for j=1:length(taxo_cell_multiple_mmR_mix2_kgr40_tr{2})
    NE_mmR_mix2_kgr40_tr_MSC1mMSC2a(i,j)=mean(mean(Sim_dis_mmR_mix2_kgr40_tr(taxo_cell_multiple_mmR_mix2_kgr40_tr{2}{i},taxo_cell_multiple_mmR_mix2_kgr40_tr{2}{j})));
    end
end

    % Annotate average distances in the Level 2 connection table
 for i=1:length(sort_min_array_MSC2)
     sort_min_array_MSC2(i,4)=NE_mmR_mix2_kgr40_tr_MSC1mMSC2a(sort_min_array_MSC2(i,1),sort_min_array_MSC2(i,2));
 end
     
    % Remove links with average distance > 0.4
sort_min_array_MSC2_cut40=sort_min_array_MSC2(find(sort_min_array_MSC2(:,4)<0.4),:);
 [ taxo_cell_out_sort_min_array_MSC2_cut40,taxo_single_node ] = pairclustering_N( sort_min_array_MSC2_cut40,length(NE_mmR_mix2_kgr40_tr{2}) );
 
    % Reconstruct clustering results:
        % Adding outlier groups
 LL=length(taxo_cell_out_sort_min_array_MSC2_cut40);
 for i=1:length(taxo_single_node)
    taxo_cell_out_sort_min_array_MSC2_cut40{LL+i}=taxo_single_node(i); 
 end

        % Mapping to node level
 clear taxo_cell_mmR_mix2_kgr40_tr_MSC3
 for i=1:length(taxo_cell_out_sort_min_array_MSC2_cut40)
     taxo_cell_mmR_mix2_kgr40_tr_MSC3{i}=[];
     
     for j=1:length(taxo_cell_out_sort_min_array_MSC2_cut40{i})
         taxo_cell_mmR_mix2_kgr40_tr_MSC3{i}=[taxo_cell_mmR_mix2_kgr40_tr_MSC3{i}',taxo_cell_multiple_mmR_mix2_kgr40_tr{2}{taxo_cell_out_sort_min_array_MSC2_cut40{i}(j)}']';
     end
 end

       % Group line information
taxo_cell_mmR_mix2_kgr40_tr_MSC3_cut40=taxo_cell_mmR_mix2_kgr40_tr_MSC3;
 clear  GL_mmR_mix2_kgr40_tr_MSC3_cut40
 for i=1:length(taxo_cell_mmR_mix2_kgr40_tr_MSC3)
     for j=1:length(taxo_cell_mmR_mix2_kgr40_tr_MSC3{i})
     GL_mmR_mix2_kgr40_tr_MSC3_cut40(taxo_cell_mmR_mix2_kgr40_tr_MSC3{i}(j))=i;
     end
 end
 
       % Remove outliers nodes in MSC1 level
 GL_mmR_mix2_kgr40_tr_MSC3_cut40_MSC1c25=GL_mmR_mix2_kgr40_tr_MSC3_cut40;
 meta_dis=Sim_dis_mmR_mix2_kgr40_tr; 
 for i=1:length(meta_dis);meta_dis(i,i)=inf;end
 meta_vv=min(meta_dis);

 GL_mmR_mix2_kgr40_tr_MSC3_cut40_MSC1c25(find(meta_vv>0.25))=0;

       % Generate MSC results at Level 3
 for i=1:length(taxo_cell_mmR_mix2_kgr40_tr_MSC3_cut40)
  taxo_cell_mmR_mix2_kgr40_tr_MSC3_cut40_MSC1c25{i}=find(GL_mmR_mix2_kgr40_tr_MSC3_cut40_MSC1c25==i);   
 end



% 5. MSC Tree Construction and Cytoscape Export
    % Generate MSC tree connection table
[ min_Dis_matrix_dp_mmR_mix2_kgr40_tr,connection_table_mmR_mix2_kgr40_tr ] = MSC_tree(  Sim_dis_mmR_mix2_kgr40_tr, taxo_cell_multiple_mmR_mix2_kgr40_tr,NE_mmR_mix2_kgr40_tr );

connection_table_mmR_mix2_kgr40_tr_recut=connection_table_mmR_mix2_kgr40_tr;

GL_mmR_mix2_kgr40_tr=zeros(length(Sim_dis_mmR_mix2_kgr40_tr),1);
for i=1:length(taxo_cell_multiple_mmR_mix2_kgr40_tr{2})
    GL_mmR_mix2_kgr40_tr(taxo_cell_multiple_mmR_mix2_kgr40_tr{2}{i})=i;
end

    % Annotate outlier links
for i=1:length(connection_table_mmR_mix2_kgr40_tr_recut)
    if connection_table_mmR_mix2_kgr40_tr_recut(i,3)>0.25
        connection_table_mmR_mix2_kgr40_tr_recut(i,4)=-1;
    end
    
    if NE_mmR_mix2_kgr40_tr_MSC1mMSC2a(GL_mmR_mix2_kgr40_tr(connection_table_mmR_mix2_kgr40_tr_recut(i,1)),GL_mmR_mix2_kgr40_tr(connection_table_mmR_mix2_kgr40_tr_recut(i,2)))>0.4
        connection_table_mmR_mix2_kgr40_tr_recut(i,4)=-1;
    end
    
end

    %  Generate network link table
clear link_table_meta
for i=1:length(connection_table_mmR_mix2_kgr40_tr_recut)
    for j=1:2
    link_table_meta{i,j}=ID_info_mmR_2439 {connection_table_mmR_mix2_kgr40_tr_recut(i,j),1};
    end
    
    for j=3:5
        link_table_meta{i,j}=connection_table_mmR_mix2_kgr40_tr_recut(i,j);
    end
end

    % Assign color codes and corresponding color list
 [ disease_color_list_out_mmR ] = get_auto_color_list(disease_list_mmR);

 for i=1:length(ID_info_mmR_2439)
     color_list_out_mmR{i,1}=ID_info_mmR_2439{i,1};
     color_list_out_mmR{i,2}=disease_color_list_out_mmR{find(strcmp(disease_color_list_out_mmR(:,1),ID_info_mmR_2439{i,3})==1),2};
 end

    % Load positional information for network nodes
 mix2_msc3_op_pos=get_cyto_pos('kgr40_mmR_mix2_cut025_op.xgmml');
 

    % Export network file for visualization in Cytoscape
    %   *Note: The output file "MSCtree_out.xml" can be visualized using Cytoscape via:
    %          File==>Import==>Network from file
[ cyto_out ] = cytoscape_transform_MSCtree_mixsimDis_bkc( link_table_meta,[ID_info_mmR_2439(:,1),ID_info_mmR_2439(:,1)],'MSCtree_out',color_list_out_mmR,1000,mix2_msc3_op_pos );


