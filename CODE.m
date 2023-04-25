close all;
clear; 
clc; 
% Max Cusick, Tomas Collado, Claudia Markel, Natalia Klim

% I found the method for reading in images here, lets be sure to cite it in the report: https://www.mathworks.com/matlabcentral/answers/396955-read-all-images-in-directory
% Read in all images
images = dir('data1Cropped\*.jpg');         

%Initial values 
n = length(images);            % Number of files found
p = 1201*901;                  % Number of pixels per image
A = ones(p, n);                % Initial matrix the size of each image as a column by the number of images
numBasis = 10;                 % Number of basis used to reconstruct images

%Loop to read in all images to matrix A
for i = 1:n
    current_image = images(i).name;                             % Get name current image
    current_image = imread(['data1Cropped\' current_image]);    % Retrieve the current image from folder
    img = im2gray(current_image);                               % Convert image to grayscale
    imgCol = img(:);                                            % Convert image to column vector
    A(:,i) = imgCol;                                            % Put image column into matrix A
end

%Find mean of A
m = mean(A,2);       

%Calculate eigenvectors and eigenvalues
[U,S,V] = svd(A - m,'econ');

%Create plot of singular values
figure
plot(diag(S));


%Reconstruct original images
IMeigen = U(:, [1:numBasis])*S(1:numBasis, 1:numBasis)*V(:, [1:numBasis])'; % apply the selected basis to the images
IM = IMeigen + m;                                                           % add the mean back to the images
IM_2D = reshape(IM, 1201, 901, []);                                         % reshape the images

%Maybe make this a for loop to show all images (?)
figure
imshow(IM_2D (:, :, 5)/255)


