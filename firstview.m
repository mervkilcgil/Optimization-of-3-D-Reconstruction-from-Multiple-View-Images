function [vSet, prevFeatures] = firstview(images,cams)
    % Undistort the first image.
    I = images{1}; 

    % Detect features. Increasing 'NumOctaves' helps detect large-scale
    % features in high-resolution images. Use an ROI to eliminate spurious
    % features around the edges of the image.
    border = 50;
    roi = [border, border, size(I, 2)- 2*border, size(I, 1)- 2*border];
    prevPoints   = detectSURFFeatures(I, 'NumOctaves', 8, 'ROI', roi);

    % Extract features. Using 'Upright' features improves matching, as long as
    % the camera motion involves little or no in-plane rotation.
    prevFeatures = extractFeatures(I, prevPoints, 'Upright', true);

    % Create an empty viewSet object to manage the data associated with each
    % view.
    vSet = viewSet;

    % Add the first view. Place the camera associated with the first view
    % and the origin, oriented along the Z-axis.
    viewId = 1;
    vSet = addView(vSet, viewId, 'Points', prevPoints, 'Orientation', ...
        eye(3, 'like', prevPoints.Location), 'Location', ...
        zeros(1, 3, 'like', prevPoints.Location));
    for i = 2:numel(images)
        % Undistort the current image.
        I = images{i};
        
        % Detect, extract and match features.
        currPoints   = detectSURFFeatures(I, 'NumOctaves', 8, 'ROI', roi);
        currFeatures = extractFeatures(I, currPoints, 'Upright', true);    
        indexPairs = matchFeatures(prevFeatures, currFeatures, ...
            'MaxRatio', .7, 'Unique',  true);
        
        % Select matched points.
        matchedPoints1 = prevPoints(indexPairs(:, 1));
        matchedPoints2 = currPoints(indexPairs(:, 2));
        
        % Estimate the camera pose of current view relative to the previous view.
        % The pose is computed up to scale, meaning that the distance between
        % the cameras in the previous view and the current view is set to 1.
        % This will be corrected by the bundle adjustment.
        [relativeOrient, relativeLoc, inlierIdx] = helperEstimateRelativePose2(...
            matchedPoints1, matchedPoints2, cams(i).camsparams);
        
        % Add the current view to the view set.
        vSet = addView(vSet, i, 'Points', currPoints);
        
        % Store the point matches between the previous and the current views.
        vSet = addConnection(vSet, i-1, i, 'Matches', indexPairs(inlierIdx,:));
        
        % Get the table containing the previous camera pose.
        prevPose = poses(vSet, i-1);
        prevOrientation = prevPose.Orientation{1};
        prevLocation    = prevPose.Location{1};
            
        % Compute the current camera pose in the global coordinate system 
        % relative to the first view.
        orientation = relativeOrient * prevOrientation;
        location    = prevLocation + relativeLoc * prevOrientation;
        vSet = updateView(vSet, i, 'Orientation', orientation, ...
            'Location', location);
        
        % Find point tracks across all views.
        tracks = findTracks(vSet);

        % Get the table containing camera poses for all views.
        camPoses = poses(vSet);

        % Triangulate initial locations for the 3-D world points.
        xyzPoints = triangulateMultiview(tracks, camPoses, cams(i).intrinsics);

        pstruct = toStruct(cams(i).camsparams);
        cams(i).camsparams = cameraParameters(pstruct);

        % Refine the 3-D world points and camera poses.
        [xyzRefinedPoints,camPoses, reprojectionErrors] = ...
        BundleAdjusmentImplementation(cams(i).camsparams, cams(i).intrinsics, camPoses, xyzPoints,tracks,i);

        % Store the refined camera poses.
        vSet = updateView(vSet, camPoses);

        prevFeatures = currFeatures;
        prevPoints   = currPoints;  
    end
end

