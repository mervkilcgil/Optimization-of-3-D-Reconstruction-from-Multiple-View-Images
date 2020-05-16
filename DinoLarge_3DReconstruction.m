clc
clear all
close all
% Use |imageDatastore| to get a list of all image file names in a
% directory.
myFolder = 'C:\Users\e1943\Documents\merve\iam566\566 project\matlab\dino';
imageDir = fullfile(myFolder, '*.png');
imds = imageDatastore(imageDir);

% Convert the images to grayscale.
images = cell(1, numel(imds.Files));
for i = 1:numel(imds.Files)
    I = readimage(imds, i);
    images{i} = rgb2gray(I);
end

camdatas = fullfile(myFolder,'dino_par.txt');
T = readtable(camdatas);
cams = cameras(T,[640,480]);
[vSet, prevFeatures] = firstview(images,cams);
[vSet, xyzPoints, camPoses, refPoses, reprojectionErrors] = computeDenseConstruction(images,vSet,cams);
dispDenseConstruction(camPoses,reprojectionErrors, xyzPoints,refPoses);
Y = double(xyzPoints);
threedpoints = array2table(Y);
writetable(threedpoints, '3DPoints.txt');