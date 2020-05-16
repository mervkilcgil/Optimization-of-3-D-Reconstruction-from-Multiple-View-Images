imOrig = imread(fullfile(matlabroot,'toolbox','vision','visiondata', ...
        'calibration','slr','image9.jpg'));
imshow(imOrig,'InitialMagnification',30);

projectedPoints = worldToImage(cameraParams,R,t,worldPoints);
hold on
plot(projectedPoints(:,1),projectedPoints(:,2),'g*-');
legend('Projected points');
hold off