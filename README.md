Superpixel segmentation using agglomerative clustering implemented as a VTK filter. VTK doesn't provide a good minheap and neither does the std namespace so Mx v1 (free of use) from CMU is included.

Note: This filter only works with 2d or 3d grayscale images. Would be easy to alter to other data types and color spaces. ie: rgb, lab, etc. This filter takes any input image, internally uses float, and outputs a float image.

You've got a few options with this filter.
- SetOutputType: Specifies how the output image should be shown.
	- AVGCOLOR: Outputs the average of the cluster.
	- LABEL: Outputs each clusters pixels with a unique value. 0, 1, 2, 3, ...
	- RANDRGB: Random rgb values. In this case uchar is outputted from the filter.
	- MAXCOLOR: The maximum of the cluster.
	- MINCOLOR: The minimum of the cluster.
- SetNumberOfSuperpixels: This is the only required input. Simply the number of superpixels you want.
- SetColorWeight: The intensities of image are scaled by this factor before use. Optionally you could just alter the brightness of the image before this filter. By increasing this you are saying you want the color to be more meaningful than the position. The merging is done via an energy definition given via spatial (position) and color components.

![Alt text](https://andaharoo.files.wordpress.com/2018/03/superpixel-random-rgb.png)

![Alt text](https://andaharoo.files.wordpress.com/2018/03/superpixel-average.png)