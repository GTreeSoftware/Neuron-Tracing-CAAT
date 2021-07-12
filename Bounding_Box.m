function [ xx_min,xx_max,yy_min,yy_max,zz_min,zz_max ] = Bounding_Box( in_points )
% find the bounding box of the input points
if size(in_points,2)~=3
    in_points=in_points';
end

xx_min=min(in_points(:,1));
xx_max=max(in_points(:,1));

yy_min=min(in_points(:,2));
yy_max=max(in_points(:,2));

zz_min=min(in_points(:,3));
zz_max=max(in_points(:,3));


end

