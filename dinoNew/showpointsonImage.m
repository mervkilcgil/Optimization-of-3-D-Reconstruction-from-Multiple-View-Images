clear all 
close all
clc
load('C:\Users\e1943\Documents\merve\iam566\566 project\matlab\dinoNew\bundle_workspace_2.mat');
I1 = rgb2gray(imread('dinoSR0001.png'));

myFolder = 'C:\Users\e1943\Documents\merve\iam566\566 project\matlab\dinoNew';
fileDir = fullfile(myFolder, 'uv_r2.txt');
T = readtable(fileDir);
rp = T{1:end,[1:end]};
ReprojectedPoints = double(reshape(rp,[length(rp)/2,2]));

fileDir = fullfile(myFolder, 'uv_g2.txt');
T = readtable(fileDir);
Points = T{1:end,[1:end]};
imgPts = double(reshape(Points,[length(Points)/2,2])); 

[m,n] = size(ReprojectedPoints);
imshow(I1);
axis on
hold on;
for i = 1:m
plot(imgPts(i,1),imgPts(i,2), 'r+', 'MarkerSize', 2, 'LineWidth', 2);
plot(ReprojectedPoints(i,1),ReprojectedPoints(i,2), 'go', 'MarkerSize', 2, 'LineWidth', 2);
end
hold off
