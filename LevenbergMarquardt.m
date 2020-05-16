
function [numIter, Anew, energyFunctions, errors] = ...
    LevenbergMarquardt(A,B,J,uv_g,numCameras,numPoints,numCamParams, imageNumber)
    import java.util.*;
    numIter        = 0;
    tol            = 1e-5;
    stol           = 1e-5;
    maxit          = 100;

    % Damping factors
    lambda          = 0.001;
    mu = 0.001;
    delta = 0.001;
    energyFunctions = zeros(maxit,1);
    energyFunction = 1000;
    e = zeros(maxit,1);;

    %the LM iteration 
    
    [numIter, Anew, energyFunctions, errors] = ...
    stepSolver(A,B,J,uv_g,numCameras,numPoints,numCamParams, mu, delta tol, stol, maxit);
%==========================================================================
function [numIter, Anew, energyFunctions, errors] = ...
    stepSolver(A,B,J,uv_g,numCameras,numPoints,numCamParams, mu, delta tol, stol, maxit)
    while numIter < maxit
        numIter = numIter + 1;
        % Fill in the Jacobian matrix
        J = fillJacobianMatrix(A,B',numCameras,numPoints,numCamParams,J);
        H = J'*J;
        [m,n] = size(J);

        uv_r = computeReprojectedPts(A, B, numCameras,numPoints);
        [uv_r, uv_g] = checkSize(uv_r, uv_g, m, n);

        d = [uv_r - uv_g];

        % Apply the damping factor to the Hessian matrix
        %H_lm = H + (lambda*eye(numJCols,numJCols));
        d0 = dp = -inv(H)*(J'*d(:));
        dL = mu*diag(diag(H));
        H_lm = H + dL;
        dp = -inv(H_lm)*(J'*d(:));
        norm_dp = norm(dp);
        if (norm_dp < tol)
            uvr = double(uv_r);
            uv_rtable = array2table(uvr);
            name = strcat('reprojectedPoints_' ,num2str(imageNumber),'.txt');
            writetable(uv_rtable, name);
            break;
        end

        % Apply updates and compute new error. The update is contingent on
        % actually reducing the reprojection error
        [ar, ac] = size(A);
        A_temp = zeros(size(A));
        for i = 1: ar
            for j = 1: ac
                A_temp(i,j) = A(i,j) + dp((i-1)*ac + j);
            end
        end
        uv_r = computeReprojectedPts(A_temp, B, numCameras,numPoints);
     
        d = [uv_r - uv_g];

        errors = d;
        e_k = 0.5*dp'*dp;
        if (norm(0.5*dp'*dp) < norm(0.5*d0'*d0))
            A = A_temp;
            mu = norm(0.5*dp'*dp);
            d0 = dp;
        else
            alpha0 = 1;
            while (alpha0 > stol)
                alpha0 = alpha0/2;
                [numIter, Anew, energyFunctions, errors] = ...
                stepSolver(A,B,J,uv_g,numCameras,numPoints,numCamParams, alpha0, delta tol, stol, maxit)
            end
        end
        energyFunctions(numIter) = (e_k);
        fprintf('Energy function is computed as %d, iteration number is %d\n',e_k, numIter);
        
        fprintf('Error is %d \n',norm_dp);
        fprintf('Lambda is %d \n',lambda);
        e(numIter)= (norm_dp);
       % drawContoursOfIteration(B, uv_r,numIter);
        
    end % end of Levenberg Marquardt
end
%==========================================================================
function [u, v] = proj3dto2d(pts3d, wx, wy, wz, tx, ty, tz, K)
    % Expression for the rotation matrix based on the Rodrigues formula
    R = AngleAxis2RotationMatrix([wx; wy; wz]);
    % Expression for the translation vector
    t=[tx;ty;tz];
    
    % perspective projection of the model point (X,Y)
    if isa(R, 'double') || isa(t, 'double')...
        || isa(pts3d, 'double')
        R = double(R);
        tTemp = double(t);                
        pts = double(pts3d);
        K = double(K);
    else
        R = single(R);
        tTemp = single(t);
        pts = single(pts3d);
        K = single(K);
    end
    t = tTemp(:)';
    cameraMatrix = [R; t] * K;
    projectedPoints = [pts ones(size(pts, 1), 1)] * cameraMatrix;
    %imagePointsTmp = bsxfun(@rdivide, projectedPoints(:, 1:2), projectedPoints(:, 3));

    uvs=projectedPoints;
    u=uvs(:,1)./uvs(:,3);
    v=uvs(:,2)./uvs(:,3);
end
%==========================================================================
% Compute the reprojected points
%==========================================================================
function uv_r = computeReprojectedPts(A, B, numCameras, numPoints)
    u_= []; v_ = [];
    for j = 1: numCameras,
        wx = A(j,1); 
        wy = A(j,2); 
        wz = A(j,3);
        tx = A(j,4); 
        ty = A(j,5); 
        tz = A(j,6);
        f = A(j,7); 
        u0 = A(j,8); 
        v0 = A(j,9);

        K = [f 0 u0; 0 f v0;  0 0 1];
        [u v] = proj3dto2d(B, wx, wy, wz, tx, ty, tz,K');
        u_ = [u_ u]; v_ = [v_ v];
    end
    u_r = reshape(u_, numPoints, numCameras)';
    u_r = u_r(:)';
    
    v_r = reshape(v_, numPoints, numCameras)';
    v_r = v_r(:)';
    % Interleave u and vs
    
    uv_r = [u_r; v_r];
    uv_r = uv_r(:);
end

%==========================================================================
% Size Check
%==========================================================================
function [uv_r, uv_g] = checkSize(uv_r, uv_g, m, n)
    [mr,nr] = size(uv_r);
    if (mr < m)
        uv_r(mr+1:m) = zeros(m-mr,1);
        disp("uv_r size increased");
    elseif (mr > m)
        uv_r = uv_r(1:m);
        disp("uv_r size decreased");
    end
    [mg,ng] = size(uv_g);
    if (mg < m)
        uv_g(mg+1:m) = zeros(m-mg,1);
        disp("uv_g size increased");
    elseif (mg > m)
        uv_g = uv_g(1:m);
        disp("uv_g size decreased");
    end
    if (mg > mr)
        uv_r(mr+1:mg) = zeros(mg-mr,1);
        disp("uv_r size increased");
    elseif (mr > mg)
        uv_g(mg+1:mr) = zeros(mr-mg,1);
        disp("uv_g size increased");
    end
end