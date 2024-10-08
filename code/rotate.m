
clear;clc
addpath rotate_code/freesurfer_matalb/;
addpath rotate_code/rotate/;
addpath data;

load('data/Shaefer.mat');


mapdata = xlsread('data/atlas_308.xlsx', 3);%注意需要比较的数据的位置
x=mapdata(:,1);
y=mapdata(:,2);

p_perm = perm_sphere_p(x,y,perm_id)

% DK coordinate
% clc;clear
% addpath code/freesurfer_matalb/;
% addpath code/rotate/;
% addpath data;
% addpath results;
% 
% path_sphere_l = 'data/lh.sphere';
% path_annot_l  = 'data/lh.aparc.annot';
% 
% coord_l=centroid_extraction_sphere(path_sphere_l,path_annot_l);
% 
% path_sphere_r = 'data/rh.sphere';
% path_annot_r  = 'data/rh.aparc.annot';
% 
% coord_r=centroid_extraction_sphere(path_sphere_r,path_annot_r);
% 
% % 设置要保存的文件名和工作表名
% filename = 'dk_coord.xlsx';
% sheetname = 'Sheet1';
% 
% % 将coord_l和coord_r按列合并为一个矩阵
% coord = [coord_l; coord_r];
% 
% % 使用xlswrite函数将矩阵保存到Excel表格中
% xlswrite(filename, coord, sheetname);

