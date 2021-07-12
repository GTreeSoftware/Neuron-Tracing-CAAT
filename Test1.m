clc;clear;close all;

%% load the segmentation result into images
image_name='F:\AdaptiveTracing_CMP\L1\21168_37926_1200_prob2.tif';
% image_name='F:\AdaptiveTracing_CMP\S3\3_prob2.tif';
% image_name='F:\AdaptiveTracing_CMP\S4\4_prob2.tif';
% image_name='F:\AdaptiveTracing_CMP\S7\7_prob2.tif';
% image_name='F:\AdaptiveTracing_CMP\S8\8_prob2.tif';
% image_name='F:\AdaptiveTracing_CMP\For_Show\6734_33262_4420_prob0.tif';

[image,image_size] = Read_3D_Tiff( image_name );

%% 1. calculate the threshold for the image
[cc_skel,Out_Label, gth1,gth_low ] = GBK_Threshold_Region( image );

%% start tracking 
tic
neurite_num=0;
K_num=cc_skel.NumObjects;
for ii=1:K_num
    % find a seed point 
    seed_ind=find(Out_Label>0 ,1,'first'); % select a prediction point as the seed
    if(isempty(seed_ind))
        break;
    end

    [X,Y,Z]=ind2sub(image_size,seed_ind);
    start_point1=[X,Y,Z];

    % voxel scooping with adaptative range for calculation
    [ Out_Label,region_cell,region_cell_num,cluster_points,node_points,connect_nodes,connect_num,connect_links,radius_points] = Ada_Voxel_Scooping1( image,Out_Label, start_point1,K_num,gth_low);

    %  build the graph and get the connection 
    line_threshold=6;
    [ Paths,branch_inds,end_inds,new_con_num,new_con_nodes ] = Build_Tracing_Graph( connect_links,connect_num,connect_nodes,line_threshold );

    branch_points=node_points(branch_inds,:);
    end_points=node_points(end_inds,:);

    % plot every branch of the lines
    idx1=find(image>=gth1);
    [xx1,yy1,zz1]=ind2sub(image_size,idx1);
    plot3(xx1,yy1,zz1,'g.');hold on;
    
    % plot the lines of the images
    for i=1:length( Paths)
        conneced_points=node_points(Paths{i},:);  
        plot3(conneced_points(:,1),conneced_points(:,2),conneced_points(:,3),'o-','MarkerSize',4,'LineWidth',3,'color',[rand() rand() rand()]);hold on;
    end
    plot3(branch_points(:,1),branch_points(:,2),branch_points(:,3),'ro','MarkerSize',9,'MarkerFaceColor','r');
    plot3(end_points(:,1),end_points(:,2),end_points(:,3),'bo','MarkerSize',9,'MarkerFaceColor','b');
 
end
toc
 








