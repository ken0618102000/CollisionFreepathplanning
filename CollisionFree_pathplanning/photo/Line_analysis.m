clc;clear all;close all;
Line = imread('Image 2.png');
Line = rgb2gray(Line);

[j_width,i_hight] = size(Line);
image_position_save = [];
image_rol = 0;
image_col = 0;
image_line = [];
count = 0;
count2 = 0;

for j=1:1:j_width
    for i=1:1:i_hight
        image_value = Line(j,i);
        
        if(image_value < 100)
            count  = count +1;
            image_position_save(1,count) = i;
        end
        center_point = round(size(image_position_save,2)/2);
        
        
    end
    if(center_point)
        count2 = count2 + 1;
        image_line(:,count2) = [image_position_save(1,center_point),j];
    end
    center_point = 0;
    count = 0;
    image_position_save = [];
end

path_straight_lenth = sqrt(((image_line(1,1)-image_line(1,size(image_line,2)))^2)+((image_line(2,1)-image_line(2,size(image_line,2)))^2));
jump = 1;
path_total_lenth = 0;

for k=1:jump:size(image_line,2)-jump
    path_total_lenth = sqrt(((image_line(1,k)-image_line(1,k+jump))^2)+((image_line(2,k)-image_line(2,k+jump))^2)) + path_total_lenth; 
end

path_efficiency = path_straight_lenth/path_total_lenth
figure('Name','Robot path','NumberTitle','off');
plot(image_line(1,:),i_hight - image_line(2,:)); hold on;
plot([image_line(1,1),image_line(1,size(image_line,2))],[i_hight - image_line(2,1),i_hight - image_line(2,size(image_line,2))]);
% 
% figure('Name','All Robot','NumberTitle','off');
% imshow(Line);