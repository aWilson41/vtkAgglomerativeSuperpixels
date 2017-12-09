Superpixel segmentation using agglomerative clustering implemented as a VTK filter. VTK doesn't provide a good minheap so Mx v1 (free of use) from CMU is included. Std's minheap can't have nodes changed without rebuilding either so it wasn't very inefficient.

You've got a few options with this filter.
- SetExcludeZero: Similar to giving a mask as secondary input. Basically any 0 value will not be merged and remain zero in the resulting image. This is useful if you only want to merge specific pixels.
- SetOutputAvg: The normal output of this filter is 0, 1, 2, 3, ... labelled superpixels. Setting this to true will instead set the color of the superpixel to the pixels average color.
- SetNumberOfSuperpixels: This one is the only required input. The number of superpixels you want.
- SetWeight: The merging of pixels is done via an energy definition. We merge the pairs with the minimal energy. The energy defintion is given via a spatial (position) and color component. This simply weights the color component giving it more of an effect in the resulting image.

Note: Exclude Zero doesn't work yet. The filter only works with 2d or 3d grayscale float images. Would be easy to alter to other data types and color spaces. ie: rgb, lab, etc.