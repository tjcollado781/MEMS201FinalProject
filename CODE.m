close all;
% Max Cusick, Tomas Collado, Claudia Markel, Natalia Klim

% I found the method for reading in images here, lets be sure to cite it in the report: https://www.mathworks.com/matlabcentral/answers/396955-read-all-images-in-directory
images = dir('data1Cropped\*.jpg');     % this reads in everything from the data1Cropped image folder as long as that folder is in your working directory.    

%set initial values based off of the images
n = length(images);    % Number of files found
p = 1201*901;    % get number of pixels per image
A = ones(p, n);     % matrix the size of each image as a column by the number of images

%this is a loop to read in all of the images and put them in A
for i = 1:n
    current_image = images(i).name;     % get the name of the current image
    current_image = imread(['data1Cropped\' current_image]);    % retrieve the current image from the actual folder
    img = im2gray(current_image);   % convert each image to grayscale
    imgCol = img(:);    % turn each grayscale image into a column
    A(:,i) = imgCol;    % throw all the images together into A
end

m = mean(A,2);        %take the mean of A

%the svd line
[U,S,V] = svd(A - m,'econ');





