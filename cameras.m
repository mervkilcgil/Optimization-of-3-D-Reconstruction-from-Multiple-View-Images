function cams = cameras(T,imageSize)
    %CAMERAS Summary of this function goes here
    %   Detailed explanation goes here

    k = T{1:end,[2:10]};
    r = T{1:end,[11:19]};
    t = T{1:end,[20:22]};
    [m,n] = size(t);

    R1 = zeros(3);
    K1 = zeros(3);
    K1(:,1) = k(1,1:3);
    K1(:,2) = k(1,4:6);
    K1(:,3) = k(1,7:end);
    R1(:,1) = r(1,1:3);
    R1(:,2) = r(1,4:6);
    R1(:,3) = r(1,7:end);
    t1 = t(1,:);
    cameraParams1 = cameraParameters('IntrinsicMatrix',K1,...
                                    'RadialDistortion',[0 0 0]); 

    cameraParams1.ImageSize = imageSize;
    paramStruct = toStruct(cameraParams1);
    paramStruct.TranslationVectors = t1;
    paramStruct.RotationVectors = rotationMatrixToVector(R1);
    cameraParams1 = cameraParameters(paramStruct);
    % Get intrinsic parameters of the camera
    focalLength1 = [cameraParams1.IntrinsicMatrix(1,1), cameraParams1.IntrinsicMatrix(2,2)];
    principalPoint1 = [cameraParams1.IntrinsicMatrix(3,1), cameraParams1.IntrinsicMatrix(3,2)];
    intrinsics1 = cameraIntrinsics(focalLength1,principalPoint1,imageSize);
    cams = struct('camsparams',cameraParams1,'intrinsics', intrinsics1,'rotations', R1,'translations', t1);
    for indx = 2:m
        R = zeros(3);
        K = zeros(3);
        K(:,1) = k(indx,1:3);
        K(:,2) = k(indx,4:6);
        K(:,3) = k(indx,7:end);
        R(:,1) = r(indx,1:3);
        R(:,2) = r(indx,4:6);
        R(:,3) = r(indx,7:end);
        cameraParams = cameraParameters('IntrinsicMatrix',K,...
                                        'RadialDistortion',[0 0 0]); 

        cameraParams.ImageSize = imageSize;
        paramStructi = toStruct(cameraParams);
        paramStructi.TranslationVectors = t(indx,:);
        paramStructi.RotationVectors = rotationMatrixToVector(R);
        cameraParams = cameraParameters(paramStructi);
        % Get intrinsic parameters of the camera
        focalLength = [cameraParams.IntrinsicMatrix(1,1), cameraParams.IntrinsicMatrix(2,2)];
        principalPoint = [cameraParams.IntrinsicMatrix(3,1), cameraParams.IntrinsicMatrix(3,2)];
        intrinsics = cameraIntrinsics(focalLength,principalPoint,imageSize);
        cams(indx).camsparams = cameraParams;
        cams(indx).intrinsics = intrinsics;
        cams(indx).rotations = R;
        cams(indx).translations = t(indx,:);
    end
end

