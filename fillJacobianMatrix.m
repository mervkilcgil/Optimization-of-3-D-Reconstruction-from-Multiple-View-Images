function J = fillJacobianMatrix(A,B,numCameras,numPoints,numCamParams,J)
    syms u0s v0s fs real
    syms txs tys tzs wxs wys wzs real
    syms Xs Ys Zs real
    % Ihe intrinsic parameter matrix
    K=[fs, 0, u0s; 0, fs, v0s; 0, 0, 1];
    % Expression for the rotation matrix based on the Rodrigues formula
    theta=sqrt(wxs^2+wys^2+wzs^2);
    omega=  [0 -wzs wys; wzs 0 -wxs; -wys wxs 0;];
    R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);
    % Expression for the translation vector
    
    t=[txs;tys;tzs];
    % perspective projection of the model point (X,Y,Z)
    
    uvs=K*[R t]*[Xs; Ys; Zs; 1];
    u=uvs(1)/uvs(3);
    v=uvs(2)/uvs(3);
    % calculate the geometric distance in x and y direction
    % u,v = the x and y positions of the projection of the corresponding model point
    dx=u;
    dy=v;
    % Evaluate the symbolic expression of the Jacobian w.r.t. the estimated parameters
    % Jx=jacobian(dx,[au,av,u0,v0,wx wy wz tx ty tz]);
    % Jy=jacobian(dy,[au,av,u0,v0,wx wy wz tx ty tz]);
    
    Jx_f = jacobian(dx,fs);             
    Jx_u0 = jacobian(dx,[u0s]);        
    Jx_v0 = jacobian(dx,[v0s]);    
    Jx_wx = jacobian(dx,[wxs]);
    Jx_wy = jacobian(dx,[wys]);
    Jx_wz = jacobian(dx,[wzs]);
    Jx_tx = jacobian(dx,[txs]);
    Jx_ty = jacobian(dx,[tys]);
    Jx_tz = jacobian(dx,[tzs]);
    Jx_X = jacobian(dx,[Xs]);
    Jx_Y = jacobian(dx,[Ys]);
    Jx_Z = jacobian(dx,[Zs]);
    
    Jy_f = jacobian(dy,fs);
    Jy_u0 = jacobian(dy,[u0s]);
    Jy_v0 = jacobian(dy,[v0s]);
    Jy_wx = jacobian(dy,[wxs]);
    Jy_wy = jacobian(dy,[wys]);
    Jy_wz = jacobian(dy,[wzs]);
    Jy_tx = jacobian(dy,[txs]);
    Jy_ty = jacobian(dy,[tys]);
    Jy_tz = jacobian(dy,[tzs]);
    Jy_X = jacobian(dy,[Xs]);
    Jy_Y = jacobian(dy,[Ys]);
    Jy_Z = jacobian(dy,[Zs]);

    % Fill in the Jacobian matrix
    for indx_j = 1: numCameras,
            wx = A(indx_j,1); 
            wy = A(indx_j,2); 
            wz = A(indx_j,3);
            tx = A(indx_j,4); 
            ty = A(indx_j,5); 
            tz = A(indx_j,6);
            f = A(indx_j,7); 
            u0 = A(indx_j,8); 
            v0 = A(indx_j,9);
            for indx_i = 1: numPoints
                X = B(1,indx_i); Y = B(2,indx_i); Z = B(3,indx_i);
                [Jx_wx1, Jx_wy1, Jx_wz1, Jx_tx1, Jx_tz1, ...
                Jx_f1, Jx_u01, Jx_v01, Jy_wx1, Jy_wy1, Jy_wz1, ...
                Jy_ty1, Jy_tz1, Jy_f1, Jy_u01, Jy_v01] = ...
                getValuesofJacobian(X,Y,Z,f,u0,v0,wx,wy,wz,tx,ty,tz);
                % x or u coordinates
                i1 = (indx_i-1)*2*numCameras+2*indx_j-1;
                j1s = (indx_j-1)*numCamParams+1;
                j1e = indx_j*numCamParams;
        
                J(i1, j1s) = Jx_wx1;
                J(i1, j1s + 1) = Jx_wy1;
                J(i1, j1s + 2) = Jx_wz1;
                J(i1, j1s + 3) = Jx_tx1;
                J(i1, j1s + 4) = 0;
                J(i1, j1s + 5) = Jx_tz1;
                J(i1, j1s + 6) = Jx_f1;
                J(i1, j1s + 7) = Jx_u01;
                J(i1, j1e) = Jx_v01;

                % y or v coordinates
                i2 = (indx_i-1)*2*numCameras+2*indx_j;
                j2s = (indx_j-1)*numCamParams+1;
                j2e = indx_j*numCamParams;
                J(i2,j2s) = Jy_wx1;
                J(i2,j2s + 1) = Jy_wy1;
                J(i2,j2s + 2) = Jy_wz1;
                J(i2,j2s + 3) = 0;
                J(i2,j2s + 4) = Jy_ty1;
                J(i2,j2s + 5) = Jy_tz1;
                J(i2,j2s + 6) = Jy_f1;
                J(i2,j2s + 7) = Jy_u01;
                J(i2,j2e) = Jy_v01;
            end
    end
end
