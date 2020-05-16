function  dispDenseConstruction(camPoses,reprojectionErrors, xyzPoints,refinedPoses)
%D�SPDENSECONSTRUCT�ON Summary of this function goes here
%   Detailed explanation goes here
% Display the refined camera poses.
figure;
%plotCamera(camPoses, 'Size', 0.2);


% Exclude noisy 3-D world points.
goodIdx = (reprojectionErrors < 1000);
R = refinedPoses.Orientation{2};
t = refinedPoses.Location{2};
plotCamera('Location', t, 'Orientation', R, 'Size', 0.2, ...
'Color', 'b', 'Label', '2', 'Opacity', 0);
hold on
% Display the dense 3-D world points.
pcshow(xyzPoints(goodIdx, :), 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
    'MarkerSize', 45);
grid on
hold off

% Specify the viewing volume.
loc1 = camPoses.Location{1};
xlim([loc1(1)-5, loc1(1)+4]);
ylim([loc1(2)-5, loc1(2)+4]);
zlim([loc1(3)-1, loc1(3)+20]);
camorbit(0, -30);

title('Dense Reconstruction');
savefig('Reconstruction.fig');
end

