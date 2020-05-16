function [Jx_wx1, Jx_wy1, Jx_wz1, Jx_tx1, Jx_tz1, ...
    Jx_f1, Jx_u01, Jx_v01, Jy_wx1, Jy_wy1, Jy_wz1, ...
    Jy_ty1, Jy_tz1, Jy_f1, Jy_u01, Jy_v01] = ...
    getValuesofJacobian(X,Y,Z,f,u0,v0,wx,wy,wz,tx,ty,tz)
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

    Xs = X;
    Ys = Y;
    Zs = Z;
    fs = f;
    u0s = u0;
    v0s = v0;
    wxs = wx;
    wys = wy;
    wzs = wz;
    txs = tx;
    tys = ty;
    tzs = tz;

    Jx_wx1 = eval(Jx_wx);
    Jx_wy1 = eval(Jx_wy);
    Jx_wz1 = eval(Jx_wz);

    Jx_tx1 = eval(Jx_tx);
    Jx_ty1 = eval(Jx_ty);
    Jx_tz1 = eval(Jx_tz);

    Jx_f1 = eval(Jx_f);
    Jx_u01 = eval(Jx_u0);
    Jx_v01 = eval(Jx_v0);

    Jy_wx1 = eval(Jy_wx);
    Jy_wy1 = eval(Jy_wy);
    Jy_wz1 = eval(Jy_wz);

    Jy_tx1 = eval(Jy_tx);
    Jy_ty1 = eval(Jy_ty);
    Jy_tz1 = eval(Jy_tz);
    
    Jy_f1 = eval(Jy_f);
    Jy_u01 = eval(Jy_u0);
    Jy_v01 = eval(Jy_v0);
end