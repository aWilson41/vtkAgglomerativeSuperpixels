Superpixel segmentation using agglomerative clustering implemented as a VTK filter. VTK doesn't provide a good minheap and neither does the std namespace so Mx v1 (free of use) from CMU is included.

You've got a few options with this filter.
- SetExcludeZero: Similar to giving a mask as secondary input. Basically any 0 value will not be merged and remain zero in the resulting image. This is useful if you only want to merge specific pixels.
- SetOutputType: The normal output of this filter is 0, 1, 2, 3, ... labeled superpixels. Optionally you can specify vtkSuperpixelFilter::RANDRGB for random rgb values to be assigned to every cluster (with the output type being a 3 component float image). Or you can specify vtkSuperpixelFilter::AVGCOLOR to average the grayscale color.
- SetNumberOfSuperpixels: This one is the only required input. The number of superpixels you want.
- SetWeight: The merging of pixels is done via an energy definition. We merge the pairs with the minimal energy. The energy defintion is given via a spatial (position) and color component. This simply weights the color component giving it more of an effect in the resulting image.

Note: The filter only works with 2d or 3d grayscale float images. Would be easy to alter to other data types and color spaces. ie: rgb, lab, etc.

![Alt text](https://andaharoo.files.wordpress.com/2017/11/superpixel-segmentation.png?w=1140)

![Alt text](https://andaharoo.files.wordpress.com/2017/12/output.png?w=354&h=236)