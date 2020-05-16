function [xyzRefinedPoints, refinedPoses, reprojectionErrors] = ...
    BundleAdjusmentImplementation(camParams, intrinsics, camPoses, objPts, tracks, imNumber)
    
    K = camParams.IntrinsicMatrix;
    [measurements, visibility, cameraMatrices, quaternionBases, cameraParams, returnPoseType] = ...
    convertInputDataFormat(tracks, camPoses, intrinsics, true);
    %poses = cameraMatrices;
    
    % Initialize LM solver
    [numPoints, numCameras] = size(visibility);
    % objPts is 3*numPoints
    
    numCamParams = 9; %
    
    numJRows = 2*numPoints*numCameras;
    numJCols = (numCamParams*numCameras);
    
    J = zeros(numJRows, numJCols);

    refinedPoses = camPoses;
    % fill B and A
    
    B = objPts;
    s = 1;
    for indx = 1:numCameras
        R = camParams.RotationMatrices;
        T = camParams.TranslationVectors;

        %[wx wy wz] = rotation_matrix_to_axis_angle(R);
        angle_axis = RotationMatrix2AngleAxis(R);
        wx = angle_axis(1);
        wy = angle_axis(2);
        wz = angle_axis(3);

        A(indx,:) = [wx, wy, wz, T(1), T(2), T(3), s*K(1,1), K(3,1), K(3,2)];
    end
    % imgPts is a numPoints*2 vector
    imgPts = measurements;
    % the "ground truth" for this optimization
    
    u_g = s*imgPts(1,:);
    v_g = s*imgPts(2,:);
    uv_g = [u_g; v_g];
    uv_g = uv_g(:);
    uvg = double(uv_g);
    uv_gtable = array2table(uvg);
    name = strcat('imagePoints_',,num2str(imNumber),'.txt');
    writetable(uv_gtable, name);
    
    [numIter, Anew, reprojectionError, errors] = ...
    LevenbergMarquardt(A,B,J,uv_g,numCameras,numPoints,numCamParams, imNumber);

    % poses is 3*4*numCameras, each pose matrix is of the form [R T]
    poses = zeros(3,4,numCameras);
    % wx, wy, wz, tx, ty, tz, f, u0, v0
    for i = 1:numCameras
        R = AngleAxis2RotationMatrix([Anew(i,1); Anew(i,2); Anew(i,3)]);
        T = [Anew(i,4) Anew(i,5) Anew(i,6)];
        poses(:, 1:3,i) = R;
        poses(:, 4, i) = T';    
        newCamParams(i, 1:3) = [Anew(i,7) Anew(i,8) Anew(i,9)];
        refinedPoses.Location{i} = cast(-T*R, returnPoseType);
        refinedPoses.Orientation{i} = cast(R, returnPoseType);
    end
    xyzRefinedPoints = poses;
    reprojectionErrors = computeReprojectionError(errors, visibility);
    name = strcat('bundle_workspace_', num2str(imNumber));
    save(name);
end
%==========================================================================
% Compute the reprojection error for each 3-D point
%==========================================================================
function reprojectionErrors = computeReprojectionError(errors, visibility)
    reprojectionErrors = zeros(size(visibility, 1), 1);
    k = 1;
    for n = 1:size(visibility, 1)
        nViews = sum(visibility(n, :));
        e = sqrt(errors(k : k + nViews - 1));
        reprojectionErrors(n) = sum(e) / numel(e);
        k = k + nViews;
    end
end
%------------------------------------------------------------------
%------------------------------------------------------------------