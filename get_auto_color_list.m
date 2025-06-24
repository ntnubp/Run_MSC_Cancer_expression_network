function [ color_list_out ] = get_auto_color_list( ID_name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

input_color_v{1}='#FF0000';
input_color_v{2}='#00FF00';
input_color_v{3}='#0000FF';
input_color_v{4}='#CCCC00';
input_color_v{5}='#FF00FF';
input_color_v{6}='#669999';
input_color_v{7}='#6666FF';
input_color_v{8}='#66664C';
input_color_v{9}='#FF6666';
input_color_v{10}='#FF9933';
input_color_v{11}='#CC66CC';
input_color_v{12}='#990033';
input_color_v{13}='#333399';
input_color_v{14}='#339933';
input_color_v{15}='#993334';
input_color_v{16}='#999933';
input_color_v{17}='#993399';
input_color_v{18}='#339999';
input_color_v{19}='#003300';
input_color_v{20}='#333333';
input_color_v{21}='#999999';
input_color_v{22}='#CCCCCC';



if length( ID_name) <=length(input_color_v)
    for i=1:length(ID_name)
    color_list_out{i,1}=ID_name{i};
    color_list_out{i,2}=input_color_v{i};    
    end
    
else
color_list_out=[];

['ID name is too long, please set a ID name list < ',num2str(length(input_color_v))]    
end


end

