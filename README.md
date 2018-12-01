Superpixel segmentation implementation of "Superpixel Generation by Agglomerative Clustering with Quadratic Error Minimization." Neither VTK or std's minheap are very good so the minheap fr om Mx v1 (free of use) from CMU is included.

https://www.utdallas.edu/~xxg061000/Superpixel_CGF2018.pdf

Note: This filter only works with 2d or 3d grayscale images. This filter takes any input image, internally uses float, and outputs a float image of the same dimension. The only big difference from the paper is that I don't use matrix representation. Without it we can easily compute the sum of square shortcut.

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

![Alt text](https://andaharoo.files.wordpress.com/2018/03/screenshot.png)
