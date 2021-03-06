clc
clear all
close all
% Use |imageDatastore| to get a list of all image file names in a
% directory.
myFolder = 'C:\Users\e1943\Documents\merve\iam566\566 project\matlab\dinoNew';
imageDir = fullfile(myFolder, '*.png');
imds = imageDatastore(imageDir);

% Convert the images to grayscale.
images = cell(1, numel(imds.Files));
for i = 1:numel(imds.Files)
    I = readimage(imds, i);
    images{i} = rgb2gray(I);
end

camdatas = fullfile(myFolder,'dinoSR_par.txt');
T = readtable(camdatas);
cams = cameras(T,[640,480]);

[vSet, prevFeatures] = firstviewtest(images,cams);
[vSet, xyzPoints, camPoses, reprojectionErrors] = computeDenseConstructiontest(images,vSet,cams);

save('dinotest');
dispDenseConstruction(camPoses,reprojectionErrors, xyzPoints,camPoses);

Y = double(xyzPoints);
threedpoints = array2table(Y);
writetable(threedpoints, '3DPointsSRTestSmall.txt');

% Visualize the scene in 3D.
ptCloud = pointCloud(refinedPoints); 
