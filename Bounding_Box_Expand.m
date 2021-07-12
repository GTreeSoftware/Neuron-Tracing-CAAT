function [ xx_min1,xx_max1,yy_min1,yy_max1,zz_min1,zz_max1 ] = Bounding_Box_Expand( in_points,image_size,ds1)
% begin to find the bounding box of the in_points
[ xx_min,xx_max,yy_min,yy_max,zz_min,zz_max ] = Bounding_Box( in_points );

% calculate the largest extednding area
xx_min1=max(1,xx_min-ds1);
yy_min1=max(1,yy_min-ds1);
zz_min1=max(1,zz_min-ds1);

xx_max1=min(image_size(1),xx_max+ds1);
yy_max1=min(image_size(2),yy_max+ds1);
zz_max1=min(image_size(3),zz_max+ds1);


end

