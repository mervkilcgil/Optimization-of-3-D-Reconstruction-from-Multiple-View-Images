function drawContoursOfIteration(imagePoints, reprojectedPoints,iter)
    L1 = sqrt(diff(imagePoints(1)).^2 + diff(imagePoints(2)).^2);
    t1 = [0 cumsum(L1)];
    t1 = linspace(0,t1(end),50);
    x1 = spline(t1,imagePoints(1),t1);            % interpolate X data
    y1 = spline(t1,imagePoints(2),t1);            % interpolate Y data

    L2 = sqrt(diff(reprojectedPoints(1)).^2 + diff(reprojectedPoints(2)).^2);
    t2 = [0 cumsum(L2)];
    t2 = linspace(0,t2(end),50);
    x2 = spline(t2,reprojectedPoints(1),t2);            % interpolate X data
    y2 = spline(t2,reprojectedPoints(2),t2);            % interpolate Y data
    figure(iter)
    plot(x1,y1)
    hold on
    plot(x2,y2)
    hold off
end