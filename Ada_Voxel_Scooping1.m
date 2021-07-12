function [label_image, region_cell,region_cell_num,cluster_points,node_points,connect_nodes,connect_num,connect_links,radius_points] = Ada_Voxel_Scooping1( prob_image,label_image, start_points,K_num,gth_low)
% get the region and subregion cell to include the the region after each expension
max_iter=4000;
max_branch=500;
% just contains the subregion and subregion size 
region_cell=cell(max_iter,max_branch);
region_cell_num=zeros(max_iter,1);

% record the cluster points of each subregion
cluster_points=zeros(max_iter*max_branch,3);
node_points=zeros(max_iter*max_branch,3);
radius_points=zeros(max_iter*max_branch,1);

% record the connection 
connect_nodes=zeros(max_iter*max_branch,max_branch);
connect_num=zeros(max_iter*max_branch,1);
connect_links=zeros(max_iter*max_branch,1);

% record each expaned area
expanded_area=zeros(max_iter*max_branch,1);

% record the cell number times
kk1=0;
dfault_gsize=[1,1,1];

% calculate every itertions
for iter_num=1:max_iter    
    % for the first time: setting the parameters
    if iter_num==1
        grown_num=1;
        grown_size=dfault_gsize;
        in_points=start_points;
        
        % begin region grown of first culster: not using the wait points
        [ out_points,~,label_image ] =Ada_Region_Grown( in_points,label_image,grown_size,grown_num );

        % judge the cluster number of every grown points based on connection
        [ sub_region,connected_size ] = Judge_Connected_Points(out_points);

        % for every clusters: find the center points of every cluster
        [ center_new ] = Calculate_SubRegion_Center( sub_region,connected_size );
             
        %%%%%%% : voxel scooping: calculate the area of the points
        [sub_region_area ] = Estimate_Points_Area(sub_region,0 );
        expanded_area((kk1+1):(kk1+connected_size),:)=sub_region_area;
        
        % record the result of the first time
        region_cell_num(1)=connected_size;
        
        for i=1:region_cell_num(iter_num)
            region_cell{iter_num,i}=sub_region{1,i};
        end
        
        cluster_points((kk1+1):(kk1+connected_size),:)=center_new;
        node_points((kk1+1):(kk1+connected_size),:)=center_new;
        radius_points((kk1+1):(kk1+connected_size),:)=1;      
        kk1=kk1+connected_size;                      
    else 
        % define the paramters
        grown_num=1;
        grown_size=dfault_gsize;
        
        % judge whether the input points were empty
        if region_cell_num(iter_num-1)==0
            break;
        end
        
        % find its last connection points index
        cc_count=sum(region_cell_num(1:iter_num-1));
        last_sub_idx=(cc_count-region_cell_num(iter_num-1)+1):cc_count;
        
        % find last subregion area
        last_sub_area=expanded_area(last_sub_idx);
            
        % for every iterations with n culsters   
        for c_iter=1:region_cell_num(iter_num-1)
            in_points=region_cell{iter_num-1,c_iter};
                    
            % begin region grown of first culster: not using the wait points
            [ out_points,~,label_image ] =Ada_Region_Grown( in_points,label_image,grown_size,grown_num );

            % whether out_points is empty or not, judge the adaptive threshold and adaptive range
            if isempty(out_points)
                % CAAT score calculation
                dt=5;
                [ cscore,npoint2]=CAAT_Score( in_points,prob_image,label_image,dt,K_num,gth_low );
                if cscore>=0.5
                    out_points=npoint2;
                    label_image(out_points(1),out_points(2),out_points(3))=-1;
                end
                
                % still could not find the points
                if isempty(out_points)
                    continue;
                end

            end
            
            % judge the cluster number of every grown points based on connection
            [ sub_region,connected_size ] = Judge_Connected_Points(out_points);
            
            % for every clusters: find the center points of every cluster
            [ center_new ] = Calculate_SubRegion_Center( sub_region,connected_size);

            % record the Ci,k and renew Ni,k
            cluster_points((kk1+1):(kk1+connected_size),:)=center_new;

            % calculate the node points of the image
            last_point=node_points(last_sub_idx(c_iter),:);          
            last_culster_area=last_sub_area(c_iter);
            
            % modify the points using voxel-scooping method
            [sub_region_area ] = Estimate_Points_Area(sub_region,1 );
            expanded_area((kk1+1):(kk1+connected_size),:)=sub_region_area;
               
            [ node_new ] = Renew_NewNodes( last_point,center_new,last_culster_area,sub_region_area,connected_size );       
            node_points((kk1+1):(kk1+connected_size),:)=node_new;     
            
            % voxel scooping to renew the subregion
            [ sub_region,label_image,radisu_record ] = Renew_SubRegion( node_new,sub_region,connected_size,label_image);

            % record the radius
            radius_points((kk1+1):(kk1+connected_size),:)=radisu_record;
 
            % record the result of the new iternum
            for i=1:connected_size          
                region_cell{iter_num,region_cell_num(iter_num)+i}=sub_region{1,i};
            end

            region_cell_num(iter_num)=region_cell_num(iter_num)+connected_size;  
            
            connect_links((kk1+1):(kk1+connected_size),:)=last_sub_idx(c_iter);

            % calculate the subregion  node links for all points
            for ss_iter=1:connected_size
                c_ind1=kk1+ss_iter;
                l_ind1=last_sub_idx(c_iter);

                % record current subregion collection
                connect_nodes(c_ind1,connect_num(c_ind1)+1)=l_ind1;
                connect_num(c_ind1)=connect_num(c_ind1)+1;

                % record last points collection
                connect_nodes(l_ind1,connect_num(l_ind1)+1)=c_ind1;
                connect_num(l_ind1)=connect_num(l_ind1)+1;

            end
            kk1=kk1+connected_size; 
        end   
    end 
end

% save the region area
region_cell_num=region_cell_num(1:iter_num-2);
region_cell=region_cell(1:iter_num-2,1:max(region_cell_num));

% save the cluster centers
cluster_points=cluster_points(1:kk1,:);
node_points=node_points(1:kk1,:);

% get the connection relation: using the connected points can find the
% longest lines of 2 end points
connect_num=connect_num(1:kk1,:);
connect_nodes=connect_nodes(1:kk1,1:max(region_cell_num));
connect_links=connect_links(1:kk1,:);
radius_points=radius_points(1:kk1,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% renew the skeleton nodes according to the paper of 3D neuron tracing
% by voxel scooping
function [ sarea ] = Estimate_Points_Area(sub_region,flag )
% judge whether the sub_region is cell or not
nsize=size(sub_region,2);
sarea=zeros(nsize,1);
for i=1:nsize
    in_points=sub_region{1,i};
    if flag==0
        [ xx_min,xx_max,yy_min,yy_max,zz_min,zz_max ] = Bounding_Box( in_points );
        new_size=[xx_max-xx_min+1,yy_max-yy_min+1,zz_max-zz_min+1];
        sarea(i)=sqrt(sum(new_size.^2));        
    else
        sarea(i)=size(in_points,1);
    end 
end
end

function [ node_new ] = Renew_NewNodes( last_point,center_new,last_area,cur_area,connected_size)
% calculte the new node points 
% refer to  "Three-dimensional neuron tracing by voxel scooping"
node_new=zeros(connected_size,3);
for i=1:connected_size
    c1=max(last_area,cur_area(i));
    c2=min(last_area,cur_area(i));  
    rate1=0.5^(c2/c1);
    node_new(i,:)=last_point+rate1*(center_new(i,:)-last_point);
end
end

function [ new_sub_region,label_image,radisu_record ] = Renew_SubRegion( node_new,sub_region,connected_size,label_image)
% renew the subregions based on voxel scooping
new_sub_region=cell(1,connected_size);
radisu_record=zeros(connected_size,1);

for i=1:connected_size
    node1= node_new(i,:);
    in_points=sub_region{1,i};
    
    % calculate the distance to other points: suppose sphere
    distance=pdist2(node1,in_points);
    con_radius=max(distance);
    con_radius=max(con_radius,1);
    radisu_record(i)=1.5*con_radius;
    
    % begin to set the range for the new region growning
    [ out_points,label_image ] = Region_Grown_Constrain( in_points,label_image,node1,con_radius );   
    % record the new points
    new_sub_region{1,i}=out_points;     
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ out_points,label_image ] = Region_Grown_Constrain( in_points,label_image,con_node,con_radius )
% define the region grown results of the images
% in_points: N*2 or N*3
% label_image: 2D or 3D ; logical: here for 3D
% out_points: all the expanded region of this time
% get the constrain of the images
image_size=size(label_image);
out_points=zeros(2000,3);
% used for different iterations
tmp=in_points;
con_radius2=con_radius^2;
grown_size=[1,1,1];

% in_points must be concluded in the new images and in label image its
% belongs to 0
pnum=size(in_points,1);
out_points(1:pnum,:)=in_points;
grown_num=3;

if length(image_size)==3
    if size(tmp,2)~=3
        tmp=tmp';
    end
    % grown times: not this one
    for i=1:grown_num
        % for all the input points
        % record the number of this iteration
        rnum=0;
        for point_num=1:size(tmp,1)
            % calculate the range of the voxels
            xx_range=max(1,tmp(point_num,1)-grown_size(1)):min(image_size(1),tmp(point_num,1)+grown_size(1));
            yy_range=max(1,tmp(point_num,2)-grown_size(2)):min(image_size(2),tmp(point_num,2)+grown_size(2));
            zz_range=max(1,tmp(point_num,3)-grown_size(3)):min(image_size(3),tmp(point_num,3)+grown_size(3));
            
            % begin to region grown
            for xx1=xx_range(1):xx_range(end)
                for yy1=yy_range(1):yy_range(end) 
                    for zz1=zz_range(1):zz_range(end)
                        if (label_image(xx1,yy1,zz1)>0)
                            dist_tmp=(xx1-con_node(1))^2+(yy1-con_node(2))^2+(zz1-con_node(3))^2;
                            if dist_tmp<=con_radius2
                                pnum=pnum+1;
                                rnum=rnum+1;
                                out_points(pnum,:)=[xx1,yy1,zz1];
                                % mark the label_image
                                label_image(xx1,yy1,zz1)=-label_image(xx1,yy1,zz1); 
                            end
                        end                
                    end
                end
            end   
        end
        
        % end of one iteration: judge whether has grown or not
        if rnum>0
            tmp=out_points((pnum-rnum+1):pnum,:);    
        else
            break;
        end 
    end   
    % end of all the iterations
    out_points=out_points(1:pnum,:);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ center_new ] = Calculate_SubRegion_Center( sub_region,connected_size )

center_new=zeros(connected_size,3);
for i=1:connected_size
    points1=sub_region{1,i};
    [ center1 ] = Calculate_Mean_Center( points1 );
    center_new(i,:)=center1;  
end
end

function [ cluster_center ] = Calculate_Mean_Center( points1 )
if size( points1,2)~=3
    points1=points1';
end

% judge whether one point or not
if size( points1,1)==1
    cluster_center=points1(1,:);
else
    % calculate the center just by x y z
    xx1=mean(points1(:,1));
    yy1=mean(points1(:,2));
    zz1=mean(points1(:,3));    
    cluster_center=[xx1,yy1,zz1];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ sub_region,connected_size ] = Judge_Connected_Points( in_points)
% input: in_points 
% output: subregion contains the size of every numbers and connected size
if size(in_points,2)~=3
    in_points=in_points';
end

% judge wheter empty or not
if isempty(in_points)
    connected_size=0;
    sub_region={};
else 
    % begin to find the bounding box of the in_points
    [ xx_min,xx_max,yy_min,yy_max,zz_min,zz_max ] = Bounding_Box( in_points );

    % crop to only contain small regions
    in_points(:,1)=in_points(:,1)-xx_min+1;
    in_points(:,2)=in_points(:,2)-yy_min+1;
    in_points(:,3)=in_points(:,3)-zz_min+1;

    new_size=[xx_max-xx_min+1,yy_max-yy_min+1,zz_max-zz_min+1];
    in_points=sub2ind(new_size,in_points(:,1),in_points(:,2),in_points(:,3));

    region_img=false(new_size);
    region_img(in_points)=1;

    % calculate the subregion based on the connection
    cc_skel1=bwconncomp(region_img);
    connected_size=cc_skel1.NumObjects;
    sub_region=cell(1,connected_size);

    % save the result
    for i=1:connected_size
        ind_tmp=cc_skel1.PixelIdxList{1,i};
        [xx1,yy1,zz1]=ind2sub(new_size,ind_tmp);
        xx1=xx1+xx_min-1;
        yy1=yy1+yy_min-1;
        zz1=zz1+zz_min-1;
        sub_region{1,i}=[xx1,yy1,zz1];

    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ out_points,wait_points,label_image ] = Ada_Region_Grown( in_points,label_image,grown_size,grown_num)
% define the region grown results of the images
% in_points: N*2 or N*3
% label_image: 2D or 3D ; logical: here for 3D
% out_points: all the expanded region of this time
% wait_points: un grown points of the images
image_size=size(label_image);
out_points=zeros(2000,3);
tmp=in_points;
pnum=0;
% for 3D image
if length(image_size)==3
    if size(tmp,2)~=3
        tmp=tmp';
    end
    % grown times
    for iter_num=1:grown_num
        % for all the input points
        % record the number of this iteration
        rnum=0;
        for point_num=1:size(tmp,1)
            % calculate the range of the voxels
            xx_range=max(1,tmp(point_num,1)-grown_size(1)):min(image_size(1),tmp(point_num,1)+grown_size(1));
            yy_range=max(1,tmp(point_num,2)-grown_size(2)):min(image_size(2),tmp(point_num,2)+grown_size(2));
            zz_range=max(1,tmp(point_num,3)-grown_size(3)):min(image_size(3),tmp(point_num,3)+grown_size(3));
            
            % begin to region grown
            for xx1=xx_range(1):xx_range(end)
                for yy1=yy_range(1):yy_range(end) 
                    for zz1=zz_range(1):zz_range(end)
                        if (label_image(xx1,yy1,zz1)>0)
                            pnum=pnum+1;
                            rnum=rnum+1;
                            out_points(pnum,:)=[xx1,yy1,zz1];                     
                            % mark the label_image: -label value as the
                            % tracked area
                            label_image(xx1,yy1,zz1)=-label_image(xx1,yy1,zz1);                        
                        end
                    end
                end
            end   
        end
        
        % end of one iteration: judge whether has grown or not
        if rnum>0
            tmp=out_points((pnum-rnum+1):pnum,:);    
        else
            break;
        end    
    end
    
    % end of all the iterations
    out_points=out_points(1:pnum,:);
    wait_points=out_points((pnum-rnum+1):pnum,:);    
end

end
