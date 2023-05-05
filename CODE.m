close all;
clear; 
clc; 
% Max Cusick, Tomas Collado, Claudia Markel, Natalia Klim

%-------------------------------------------------------------------------
% All images

% Read in all images
imagesR = dir('data1Cropped\*.jpg');       

%Initial values 
n = length(imagesR);           % Number of files found
p = 1201*901;                  % Number of pixels per image
A = ones(p, n);                % Initial matrix the size of each image as a column by the number of images
numBasis = 20;                  % Number of basis used to reconstruct images

%Loop to read in all images to matrix A
for i = 1:n
    current_image = imagesR(i).name;                            % Get name current image
    current_image = imread(['data1Cropped\' current_image]);     % Retrieve the current image from folder
    img = im2gray(current_image);                               % Convert image to grayscale
    imgCol = img(:);                                            % Convert image to column vector
    A(:,i) = imgCol;                                            % Put image column into matrix A
end

%Find mean of A
mR = mean(A,2);       

%Calculate eigenvectors and eigenvalues
[U,S,V] = svd(A - mR,'econ');

%Create plot of singular values
plot(diag(S));
hold on; 
title('Singular values of All Eigenfaces');
grid on; 
hold off;

%plot V for first image
plot(V(1,1:20))
grid on
title('V values for first reconstructed eigenface')
grid on
xlabel('Principal Component')
ylabel('Corresponding V-Value')

%Reconstruct original images
IMeigen = U(:, [1:numBasis])*S(1:numBasis, 1:numBasis)*V(:, [1:numBasis])';  % apply the selected basis to the images
IM = IMeigen + mR;                                                           % add the mean back to the images
IM_2D = reshape(IM, 1201, 901, []);                                          % reshape the images

%Show all the Eigenfaces
implay(IM_2D/255) 


%-------------------------------------------------------------------------
%Resting images

% Read in all images
imagesR = dir('dataResting\*.jpg');   % folder of all resting face images      

%Initial values 
n = length(imagesR);           % Number of files found
p = 1201*901;                  % Number of pixels per image
A = ones(p, n);                % Initial matrix the size of each image as a column by the number of images
numBasis = 8;                 % Number of basis used to reconstruct images

%Loop to read in all images to matrix A
for i = 1:n
    current_image = imagesR(i).name;                            % Get name current image
    current_image = imread(['dataResting\' current_image]);     % Retrieve the current image from folder
    img = im2gray(current_image);                               % Convert image to grayscale
    imgCol = img(:);                                            % Convert image to column vector
    A(:,i) = imgCol;                                            % Put image column into matrix A
end

%Find mean of A
mR = mean(A,2);       

%Calculate eigenvectors and eigenvalues
[U,S,V] = svd(A - mR,'econ');

%Create plot of singular values
figure
plot(diag(S));
hold on; 
title('Singular values of Straight Eigenfaces');
xlabel('Index of Singular Value');
ylabel('Singular Value');
grid on; 
hold off;

%Reconstruct original images
IMeigen = U(:, [1:numBasis])*S(1:numBasis, 1:numBasis)*V(:, [1:numBasis])';  % apply the selected basis to the images
IM = IMeigen + mR;                                                           % add the mean back to the images
IM_2D = reshape(IM, 1201, 901, []);                                          % reshape the images

%Show all the Eigenfaces
implay(IM_2D/255) 


%-------------------------------------------------------------------------
% Smiling images

% Read in all Smiling images
imagesS = dir('dataSmiling\*.jpg'); % folder of all smiling images

%Initial values 
n = length(imagesS);           % Number of files found
p = 1201*901;                  % Number of pixels per image
A = ones(p, n);                % Initial matrix the size of each image as a column by the number of images
numBasis = 8;                 % Number of basis used to reconstruct images

%Loop to read in all images to matrix A
for i = 1:n
    current_image = imagesS(i).name;                            % Get name current image
    current_image = imread(['dataSmiling\' current_image]);     % Retrieve the current image from folder
    img = im2gray(current_image);                               % Convert image to grayscale
    imgCol = img(:);                                            % Convert image to column vector
    A(:,i) = imgCol;                                            % Put image column into matrix A
end

%Find mean of A
mS = mean(A,2);       

%Calculate eigenvectors and eigenvalues
[U,S,V] = svd(A - mS,'econ');

%Plot V value for first Eigenface
figure
plot(V(1,:));
hold on;
title('V values of First Reconstructed Smiling Eigenface');
grid on; 
xlabel('Principal Component')
ylabel('Corresponding V-Value')
hold off; 

%Create plot of singular values
figure
plot(diag(S));
hold on; 
title('Singular values of Smiling Eigenfaces');
xlabel('Index of Singular Value');
ylabel('Singular Value');
grid on; 
hold off;

%Reconstruct original images
IMeigen = U(:, [1:numBasis])*S(1:numBasis, 1:numBasis)*V(:, [1:numBasis])';  % apply the selected basis to the images
IM = IMeigen + mS;                                                           % add the mean back to the images
IM_2D = reshape(IM, 1201, 901, []);                                          % reshape the images

%Show all the Eigenfaces
implay(IM_2D/255) 


%-------------------------------------------------------------------------
% Frowning images

% Read in all Frowning images
imagesF = dir('dataFrowning\*.jpg'); % folder of all frowning images

%Initial values 
n = length(imagesF);           % Number of files found
p = 1201*901;                  % Number of pixels per image
A = ones(p, n);                % Initial matrix the size of each image as a column by the number of images
numBasis = 8;                 % Number of basis used to reconstruct images

%Loop to read in all images to matrix A
for i = 1:n
    current_image = imagesF(i).name;                            % Get name current image
    current_image = imread(['dataFrowning\' current_image]);    % Retrieve the current image from folder
    img = im2gray(current_image);                               % Convert image to grayscale
    imgCol = img(:);                                            % Convert image to column vector
    A(:,i) = imgCol;                                            % Put image column into matrix A
end

%Find mean of A
mF = mean(A,2);       

%Calculate eigenvectors and eigenvalues
[U,S,V] = svd(A - mF,'econ');

%Plot V value for first Eigenface
figure
plot(V(1,:));
hold on;
title('V values of First Reconstructed Frowning Eigenface');
xlabel('Principal Component')
ylabel('Corresponding V-Value')
grid on; 
hold off; 

%Create plot of singular values
figure
plot(diag(S));
hold on; 
title('Singular values of Frowning Eigenfaces');
xlabel('Index of Singular Value');
ylabel('Singular Value');
grid on; 
hold off;

%Reconstruct original images
IMeigen = U(:, [1:numBasis])*S(1:numBasis, 1:numBasis)*V(:, [1:numBasis])';  % apply the selected basis to the images
IM = IMeigen + mF;                                                           % add the mean back to the images
IM_2D = reshape(IM, 1201, 901, []);                                          % reshape the images

%Show all the Eigenfaces
implay(IM_2D/255) 


