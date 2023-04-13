close all;

%need to find a good way to import all the images
n = 55;
A = ones(size of each image columned, n)

%use a for loop to transform and concatinate each image
for i = i:n
    img = imread(i);
    img = im2gray(img);
    imgCol = img(:);
    A = [A imgcol];
end
m = mean(A,z);
A_modified = A - m;

%the svd line
[U,S,V] = svd(A_modified,'econ');

