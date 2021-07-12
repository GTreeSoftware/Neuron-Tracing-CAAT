function [cc_skel,Out_Label, gth1,gth_low ] = GBK_Threshold_Region( image1 )
% this is used to estimate the global threshold for the whole images: [0 255]
%% 1. calculate the bounding box of the images
value1=image1(image1>0);
[ ~,count1 ] = Calculate_Histogram( value1,[1,128] );
cum_prob1=cumsum(count1(1:127)/sum(count1(1:127)));

gth1=find(cum_prob1>=0.7, 1 );
gth1=3*(gth1+1);
pred1=image1>=gth1;

% analysis the dataset: contain seed points, 
cc_skel=bwconncomp(pred1);
Out_Label=labelmatrix(cc_skel);
Out_Label=double(Out_Label);

% calculate the lowest value for min 
gth_low=find(cum_prob1>=0.5, 1 );

end

function [ center1,count1 ] = Calculate_Histogram( value1,range1 )
center1=range1(1):range1(2);
count1=zeros(1,range1(2)-range1(1)+1);

for i=1:length(value1)
    num1=value1(i)-range1(1)+1;
    if value1(i)<=range1(2)
        count1(num1)=count1(num1)+1;
    end
end

end