Superpixel segmentation using agglomerative clustering implemented as a VTK filter. VTK doesn't provide a good minheap and neither does the std namespace so Mx v1 (free of use) from CMU is included.

You've got a few options with this filter.
- SetOutputType: The normal output of this filter is 0, 1, 2, 3, ... labeled superpixels. Optionally you can specify vtkSuperpixelFilter::RANDRGB for random rgb values to be assigned to every cluster (with the output type being a 3 component float image). Or you can specify vtkSuperpixelFilter::AVGCOLOR to average the grayscale color of each cluster.
- SetNumberOfSuperpixels: This is the only required input. The number of superpixels you want.
- SetWeight: The intensities of image are scaled by this factor before use. Optionally you could just alter the brightness of the image before this filter. By increasing this you are saying you want the color to be more meaningful than the position. The merging is done via an energy definition given via spatial (position) and color components.

Note: The filter only works with 2d or 3d grayscale float images. Would be easy to alter to other data types and color spaces. ie: rgb, lab, etc.

![Alt text](https://andaharoo.files.wordpress.com/2017/11/superpixel-segmentation.png?w=1140)

![Alt text](https://andaharoo.files.wordpress.com/2017/12/output.png?w=354&h=236)