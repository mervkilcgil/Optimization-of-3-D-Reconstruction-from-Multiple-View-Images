# Optimization-of-3-D-Reconstruction-from-Multiple-View-Images
In this project, the poses of a calibrated camera from a sequence of views are estimated,
and the 3-D structure of the scene is reconstructed up to an unknown scale. We use
the pairwise point matches to estimate the camera pose of the current view relative to the
previous view. It then links the pairwise matches into longer point tracks spanning multiple
views using the findTracks method of the viewSet object. These tracks then serve as inputs
to multiview triangulation using the triangulateMultiview function and the refinement of
camera poses and the 3-D scene points using the bundle adjustment method.
