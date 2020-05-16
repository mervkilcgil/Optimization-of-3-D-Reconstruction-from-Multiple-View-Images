function [relativeOrient, relativeLoc, inlierIdx] = helperEstimateRelativePose(...
                            matchedPoints1, matchedPoints2, cameraParams,R,t);
[E,inlierIdx] = modifestimateEssentialMatrix(matchedPoints1,matchedPoints2,...
 cameraParams);
 relativeOrient = R';
 relativeLoc = -t * relativeOrient;
end

