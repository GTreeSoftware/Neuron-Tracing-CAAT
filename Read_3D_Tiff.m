function [ Image,image_size ] = Read_3D_Tiff( image_name )
% read the tiff file into images
Info=imfinfo(image_name);

Slice=size(Info,1);                                                    %%获取图片z向帧�?Width=Info.Width;
Height=Info.Height;
Width=Info.Width;
Image=zeros(Height,Width,Slice);
image_size=[Height,Width,Slice];
for i=1:Slice
    Image(:,:,i)=imread(image_name,i);                                         %%�?���?��的读入图�?end
end

end

