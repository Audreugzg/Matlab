clc;
clear;
close all;

img = imread('la_la_land.jpg');
gray_pix = im2double(rgb2gray(img));
% change to gray scale
[row,column,k] = size(img);
I(:,:,1) = gray_pix;
I(:,:,2) = gray_pix;
I(:,:,3) = gray_pix;
%imshow(I);
I = im2double(I);
figure;
imshow(I);
title('Original Image  ');

sample_size = [201,201]; % for the sake of graidant filter calculation
sample_ij = [100,200];



%%%%%% sum of I %%%%%%%
% sum_row(1) is the first element of the row 
sum_row(1:row,1) = gray_pix(1:row,1);
sum_column(1,1:column) = gray_pix(1,1:column);

% get sum of each row 
for a = 1:row
    for b = 2:column
        sum_row(a,b) = sum_row(a,b-1) + gray_pix(a,b);
    end
end
% get the sum of each column
for a = 2:row
    for b = 1:column
        sum_column(a,b) = sum_column(a-1,b) + gray_pix(a,b);
    end
end

sum_table = zeros(row,column);
sum_table(1:row,1) = sum_column(1:row,1);
sum_table(1,1:column) = sum_row(1,1:column);

for a = 2:row
    for b = 2:column
        sum_table(a,b) =  sum_table(a-1,b-1) + sum_column(a-1,b) +sum_row(a,b-1)+gray_pix(a,b);
    end
end
   

%%%%%% sum of I^2 %%%%%%%
% sum_row(1) is the first element of the row 
sum_row(1:row,1) = gray_pix(1:row,1).*gray_pix(1:row,1);
sum_column(1,1:column) = gray_pix(1,1:column).*gray_pix(1,1:column);

% get the sum of each row
for a = 1:row
    for b = 2:column
        sum_row(a,b) = sum_row(a,b-1) + gray_pix(a,b).*gray_pix(a,b);
    end
end
% get the sum of each column
for a = 2:row
    for b = 1:column
        sum_column(a,b) = sum_column(a-1,b) + gray_pix(a,b).*gray_pix(a,b);
    end
end

sum_table_I_sqr = zeros(row,column);
sum_table_I_sqr(1:row,1) = sum_column(1:row,1);
sum_table_I_sqr(1,1:column) = sum_row(1,1:column);

for a = 2:row
    for b = 2:column
        sum_table_I_sqr(a,b) =  sum_table_I_sqr(a-1,b-1) + sum_column(a-1,b) +sum_row(a,b-1)+gray_pix(a,b).*gray_pix(a,b);
    end
end
   


%%%%%% sum of iI %%%%%%%


% sum_row(1) is the first element of the row 
sum_row(1:row,1) = gray_pix(1:row,1).*(1:row)';
sum_column(1,1:column) = gray_pix(1,1:column);

% get sum of each row 
for a = 1:row
    for b = 2:column
        sum_row(a,b) = sum_row(a,b-1) + a*gray_pix(a,b);
    end
end
% get the sum of each column
for a = 2:row
    for b = 1:column
        sum_column(a,b) = sum_column(a-1,b) + a*gray_pix(a,b);
    end
end

sum_table_i = zeros(row,column);
sum_table_i(1:row,1) = sum_column(1:row,1);
sum_table_i(1,1:column) = sum_row(1,1:column);

for a = 2:row
    for b = 2:column
        sum_table_i(a,b) =  sum_table_i(a-1,b-1) + sum_column(a-1,b) +sum_row(a,b-1)+a*gray_pix(a,b);
    end
end




%%%%%% sum of jI %%%%%%%


% sum_row(1) is the first element of the row 
sum_row(1:row,1) = gray_pix(1:row,1);
sum_column(1,1:column) = gray_pix(1,1:column).*(1:column);

% get sum of each row 
for a = 1:row
    for b = 2:column
        sum_row(a,b) = sum_row(a,b-1) + b*gray_pix(a,b);
    end
end
% get the sum of each column
for a = 2:row
    for b = 1:column
        sum_column(a,b) = sum_column(a-1,b) + b*gray_pix(a,b);
    end
end

sum_table_j = zeros(row,column);
sum_table_j(1:row,1) = sum_column(1:row,1);
sum_table_j(1,1:column) = sum_row(1,1:column);

for a = 2:row
    for b = 2:column
        sum_table_j(a,b) =  sum_table_j(a-1,b-1) + sum_column(a-1,b) +sum_row(a,b-1)+ b*gray_pix(a,b);
    end
end



%%%%% all 4 sum tables are ready 
    
% get sample variance 
sampleIsqr =  1/(sample_size(1)*sample_size(2)) * (sum_table_I_sqr(sample_ij(1)+sample_size(1), sample_ij(2)+sample_size(2)) + sum_table_I_sqr(sample_ij(1),sample_ij(2)) - sum_table_I_sqr(sample_ij(1),sample_ij(2)+sample_size(2)) - sum_table_I_sqr(sample_ij(1)+sample_size(1),sample_ij(2)));
sampleI = 1/(sample_size(1)*sample_size(2)) * (sum_table(sample_ij(1)+sample_size(1), sample_ij(2)+sample_size(2)) + sum_table(sample_ij(1),sample_ij(2)) - sum_table(sample_ij(1),sample_ij(2)+sample_size(2)) - sum_table(sample_ij(1)+sample_size(1),sample_ij(2)));
sample_var = sampleIsqr - sampleI^2;

% get sample gradient 
%sample_size = [201,201]; % for the sake of graidant filter calculation
%sample_ij = [200,500];

sampleGx = (sum_table_i(sample_ij(1)+sample_size(1), sample_ij(2)+sample_size(2)) + sum_table_i(sample_ij(1),sample_ij(2)) - sum_table_i(sample_ij(1),sample_ij(2)+sample_size(2)) - sum_table_i(sample_ij(1)+sample_size(1),sample_ij(2)))
sampleGx = 1/sampleGx * (sampleGx - (sample_ij(1)+(sample_size(1)+1)/2)*sampleI*(sample_size(1)*sample_size(2)))

sampleGy = (sum_table_j(sample_ij(1)+sample_size(1), sample_ij(2)+sample_size(2)) + sum_table_j(sample_ij(1),sample_ij(2)) - sum_table_j(sample_ij(1),sample_ij(2)+sample_size(2)) - sum_table_j(sample_ij(1)+sample_size(1),sample_ij(2)))
sampleGy = 1/sampleGy * (sampleGy - (sample_ij(2)+(sample_size(2)+1)/2)*sampleI*(sample_size(1)*sample_size(2)))

sampleGmag = (sampleGx^2 + sampleGy^2)^0.5
sampleGdir = atan2(sampleGy,sampleGx)

figure;
imshow(gray_pix(sample_ij(1):(sample_ij(1)+sample_size(1)),sample_ij(2):sample_ij(2)+sample_size(2)));
title('Sample Image ');
counter = 0;


% use variance to find patch
for i = 1:row-sample_size(1)
    for j = 1:column -sample_size(2)
         % var = sum(x^2)/N - miu^2
        Isqr = 1/(sample_size(1)*sample_size(2)) * (sum_table_I_sqr(i+sample_size(1), j+sample_size(2)) + sum_table_I_sqr(i,j) - sum_table_I_sqr(i,j+sample_size(2)) - sum_table_I_sqr(i+sample_size(1),j));
        I = 1/(sample_size(1)*sample_size(2)) * (sum_table(i+sample_size(1), j+sample_size(2)) + sum_table(i,j) - sum_table(i,j+sample_size(2)) - sum_table(i+sample_size(1),j));
        patch_var = Isqr - I^2; 
                if((abs(patch_var - sample_var))==0)
                    figure
                     display_img = gray_pix(i:i+sample_size(1),j:j+sample_size(2));
                     imshow(display_img);
                    
                     title('Matched Image using Variance ');
                end
    end
end


% use gradient to find patch 
for i = 1:row-sample_size(1)
    for j = 1:column -sample_size(2)
         % var = sum(x^2)/N - miu^2
        Gx = (sum_table_i(i+sample_size(1), j+sample_size(2)) + sum_table_i(i,j) - sum_table_i(i,j+sample_size(2)) - sum_table_i(i+sample_size(1),j));
        I = (sum_table(i+sample_size(1), j+sample_size(2)) + sum_table(i,j) - sum_table(i,j+sample_size(2)) - sum_table(i+sample_size(1),j));
        patch_Gx = 1/Gx *(Gx - (i +(sample_size(1)+1)/2)*I ); 
        
        Gy = (sum_table_j(i+sample_size(1), j+sample_size(2)) + sum_table_j(i,j) - sum_table_j(i,j+sample_size(2)) - sum_table_j(i+sample_size(1),j));
        patch_Gy = 1/Gy *(Gy - (j +(sample_size(2)+1)/2)*I ); 
                
        Gmag = (patch_Gx^2 + patch_Gy^2)^0.5;
        Gdir = atan2(patch_Gy,patch_Gx);

        
                if(abs(Gmag - sampleGmag)< 1e-4 && abs(Gdir - sampleGdir)<1e-4)
                     figure
                     display_img = gray_pix(i:i+sample_size(1),j:j+sample_size(2));
                     imshow(display_img);
                    
                     title('Matched Image using Gradient ');
                end
    end
end




