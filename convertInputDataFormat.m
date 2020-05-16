%==========================================================================
% Convert inputs to internal data structure.
%
% measurements: 2-by-M packed 2-D points
% visibility(i,j): true if point i is visible in view j
% cameraMatrices: 6-by-V, rotation vector+ translation vector
% quaternionBases: 4-by-V, quaternions for initial rotations
% cameraParams: a structure of camera parameters
%
% Note, all returned values are double
%==========================================================================
function [measurements, visibility, cameraMatrices, quaternionBases, ...
    cameraParamStruct, returnType] = convertInputDataFormat(pointTracks, ...
    cameraPoses, cameraParams, isUndistorted)

    numPoints = numel(pointTracks);
    numViews = height(cameraPoses);

    % visibility(i,j): true if point i is visible in view j
    visibility = zeros(numPoints, numViews);
    viewIds = cameraPoses.ViewId;
    x = zeros(numPoints, numViews);
    y = zeros(numPoints, numViews);
    for m = 1:numPoints
        trackViewIds = pointTracks(m).ViewIds;
        for n = 1:length(trackViewIds)
            viewIndex = find(viewIds == trackViewIds(n), 1, 'first');
            if isempty(viewIndex)
                error(message('vision:absolutePoses:missingViewId', trackViewIds(n)));
            end
            visibility(m, viewIndex) = 1;
            x(m, viewIndex) =  pointTracks(m).Points(n, 1);
            y(m, viewIndex) =  pointTracks(m).Points(n, 2);
        end
    end

    isVisible = find(visibility);
    x = x(isVisible);
    y = y(isVisible);

    visibility = sparse(visibility);
    % measurements stores 2-D points in 1st view first, then 2nd view, ...
    measurements = double([x, y])';

    % Convert camera poses to a compact form of camera projection matrices
    % Use quaternion for numerical stability
    quaternionBases = zeros(4, numViews);
    cameraMatrices = zeros(6, numViews);
    for j = 1:numViews
        t = cameraPoses.Location{j};
        R = cameraPoses.Orientation{j};
        cameraMatrices(4:6, j) = -t*R';
        quaternionBases(:, j) = vision.internal.quaternion.rotationToQuaternion(R);
    end

    returnType = class(cameraPoses.Location{1});

    if ~iscell(cameraParams)
        cameraParams = {cameraParams};
    end

    numCameras = length(cameraParams);
    cameraParamStruct(numCameras) = struct('focalLength',[], ...
                                        'principalPoint', [], ...
                                        'radialDistortion', [], ...
                                        'tangentialDistortion', [], ...
                                        'skew', []);

    for n = 1:numCameras
        % Convert the cameraParams object to a simple structure
        cameraParamStruct(n).focalLength = double(cameraParams{n}.FocalLength);
        cameraParamStruct(n).principalPoint = double(cameraParams{n}.PrincipalPoint);
        if ~isUndistorted
            % Skip if the distortion coefficients are all zeros
            if (any(cameraParams{n}.RadialDistortion) || any(cameraParams{n}.TangentialDistortion))
                cameraParamStruct(n).radialDistortion = double(cameraParams{n}.RadialDistortion);
                cameraParamStruct(n).tangentialDistortion = double(cameraParams{n}.TangentialDistortion);
            end
        end
        % Note, the internal reprojection function uses a different definition
        % of skew factor, i.e., s = S / fc(1)
        cameraParamStruct(n).skew = double(cameraParams{n}.Skew / cameraParams{n}.FocalLength(1));
    end
end