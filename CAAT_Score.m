function [ cscore,npoint2] = CAAT_Score( in_points,prob_img,label_image,dt,K_num,th_low )
% calculateconnected score using distance, region connectivity, and continuity term
% 1. calculate whether there are unvisited points in the searching range
RC_Value =  Region_Connectivity_Term(in_points,label_image,dt,K_num);
cscore=0;
npoint2=[0,0,0];

% 2. begin to calculate CAAT score
if RC_Value>0
    % 3. calculate the distance term
    [ DC_Value,npoint2 ] = Distance_Term(in_points,label_image,dt,K_num);
    
    % 4. calculate the direction continuity term
    th_low=min(0.1,th_low);
    DPC_Value= Direction_probablity_continuity_Term(in_points,prob_img,label_image,dt,K_num,th_low);
    cscore=DC_Value*RC_Value*DPC_Value;    
end

end


function [ RC_Value ] = Region_Connectivity_Term(in_points,label_image,dt,K_num)
% this is used to find whether there exist unvisited pixels that belong to
% another region.
RC_Value=0;
if size(in_points,2)~=3
    in_points=in_points';
end

% begin to find the bounding box of the in_points
image_size=size(label_image);
[ xx_min1,xx_max1,yy_min1,yy_max1,zz_min1,zz_max1 ] = Bounding_Box_Expand( in_points,image_size,dt+3);
new_label1=label_image(xx_min1:xx_max1,yy_min1:yy_max1,zz_min1:zz_max1);
% renew inputs
in_points(:,1)=in_points(:,1)-xx_min1+1;
in_points(:,2)=in_points(:,2)-yy_min1+1;
in_points(:,3)=in_points(:,3)-zz_min1+1;
npoint1=in_points(1,:);

% judge whether there are unvisted points in nearby place
new_pind=find(new_label1>0 & new_label1<=K_num);
if ~isempty(new_pind)
    new_value1=new_label1(new_pind);    
    in_value1=abs(new_label1(npoint1(1),npoint1(2),npoint1(3)));
    in_ind1=find(new_value1~=in_value1);
    
    % whether there are new points need to be analysed
    if ~isempty(in_ind1)
        RC_Value=1;
    end   
end
end


function [ DC_Value,npoint2 ] = Distance_Term(in_points,label_image,dt,K_num)
% this is used to find whether there exist unvisited pixels that belong to
% another region.
if size(in_points,2)~=3
    in_points=in_points';
end

% begin to find the bounding box of the in_points
image_size=size(label_image);
[ xx_min1,xx_max1,yy_min1,yy_max1,zz_min1,zz_max1 ] = Bounding_Box_Expand( in_points,image_size,dt+3);
new_label1=label_image(xx_min1:xx_max1,yy_min1:yy_max1,zz_min1:zz_max1);
% begin to decide whether link or not
% renew inputs
in_points(:,1)=in_points(:,1)-xx_min1+1;
in_points(:,2)=in_points(:,2)-yy_min1+1;
in_points(:,3)=in_points(:,3)-zz_min1+1;
points1=in_points(1,:);
new_size=size(new_label1);
index1=sub2ind(new_size,points1(:,1),points1(:,2),points1(:,3));

% find the index of new searched points
index2=find(new_label1>0 & new_label1<=K_num);
[npoint1,npoint2]=Nearest_3D_Points(index1,index2,new_size);

% calculate their distance value
distance1=max(abs(npoint1-npoint2));
if distance1<=dt
    DC_Value=1;
else
    DC_Value=exp(-(distance1-dt)/3);
end

npoint2(1)=npoint2(1)+xx_min1-1;
npoint2(2)=npoint2(2)+yy_min1-1;
npoint2(3)=npoint2(3)+zz_min1-1;

end


function [DPC_Value] = Direction_probablity_continuity_Term(in_points,prob_img,label_image,dt,K_num,th_low)
% this is used to find whether there exist unvisited pixels that belong to
% another region.
if size(in_points,2)~=3
    in_points=in_points';
end

% begin to find the bounding box of the in_points
image_size=size(label_image);
[ xx_min1,xx_max1,yy_min1,yy_max1,zz_min1,zz_max1 ] = Bounding_Box_Expand( in_points,image_size,dt+3);
new_label1=label_image(xx_min1:xx_max1,yy_min1:yy_max1,zz_min1:zz_max1);
new_prob1=prob_img(xx_min1:xx_max1,yy_min1:yy_max1,zz_min1:zz_max1);
% begin to decide whether link or not
% renew inputs
in_points(:,1)=in_points(:,1)-xx_min1+1;
in_points(:,2)=in_points(:,2)-yy_min1+1;
in_points(:,3)=in_points(:,3)-zz_min1+1;
points1=in_points(1,:);
new_size=size(new_label1);
index1=sub2ind(new_size,points1(:,1),points1(:,2),points1(:,3));

% find the index of new searched points
index2=find(new_label1>0 & new_label1<=K_num);
[npoint1,npoint2]=Nearest_3D_Points(index1,index2,new_size);

% calculate their probability along line
[ new_points] = Points_In_Line(npoint1,npoint2,new_size);
index3=sub2ind(new_size,new_points(:,1),new_points(:,2),new_points(:,3));
CP=new_prob1(index3);
CP(CP>=th_low)=1;
Pd=(length(CP)-sum(CP))/length(CP);
DPC_Value=exp(-Pd);

end

function [ new_points] = Points_In_Line(npoint1,npoint2,new_size)
[value1,direction1]=max(abs(npoint1-npoint2));
new_points=zeros(value1+1,3);

tt_min=min(npoint1(direction1),npoint2(direction1));
% tt_max=max(npoint1(direction1),npoint2(direction1));

for i=1:(value1+1)
    % first find the direction
    new_points(i,direction1)=tt_min+i-1;
    for j=1:3
        if j~=direction1
            aa1=(new_points(i,direction1)-npoint1(direction1))/(npoint2(direction1)-npoint1(direction1));
            aa2=npoint1(j)+aa1*(npoint2(j)-npoint1(j));
            new_points(i,j)=round(aa2);
            new_points(i,j)=max(1,new_points(i,j));
            new_points(i,j)=min(new_points(i,j),new_size(j));
        end       
    end 
end
end

function [ npoints1,npoints2 ] = Nearest_3D_Points( index1,index2,isize )
[xx1,yy1,zz1]=ind2sub(isize,index1);
[xx2,yy2,zz2]=ind2sub(isize,index2);

points1=[xx1,yy1,zz1];
points2=[xx2,yy2,zz2];

% find the nearest points in points2 for every points in points1 
[idx1, dist1] = knnsearch(points2,points1);
[min_dist,min_dind]=min(dist1);
npoints1=points1(min_dind,:);
npoints2=points2(idx1(min_dind),:);
end
