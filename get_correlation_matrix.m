function [Sim_cor_out] = get_correlation_matrix(data_matrix_in)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
meta_ss=size(data_matrix_in);
Sim_cor_out=zeros(meta_ss(2));

tic
for i=2:meta_ss(2)
    for j=1:i-1
        vv1=data_matrix_in(:,i);
        vv2=data_matrix_in(:,j);
        
        [AA BB]=corrcoef(vv1,vv2);
        Sim_cor_out(i,j)=AA(1,2);
    end
end

Sim_cor_out=Sim_cor_out+Sim_cor_out';

for i=1:length(Sim_cor_out)
    Sim_cor_out(i,i)=1;
end

end