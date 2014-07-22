Microdissection-Automation
==========================

Towards Computer Assisted Laser Microdissection &amp; Segmentation Using Matlab

Matlab version requirements: Image Processing Toolbox, Statistics Toolbox

1. set your working directory to destination folder where has all your images
2 a). for at least three-color picture: 
	cordcut=Run_ALL(imagename, nColors, filename, linTol, dilation_shape, dilation_size, auto-contrast);

sample input:
	imagename='adjusted.bmp'; //your color image name
	
	cordcut=Run_ALL(imagename,3,'ori_file.xml',10,'disk',5, 1); 
	
	% this number 10 can be changed. Smaller number means low reduce rate and lead to more accurate contours.
	
 	% choose your cluster of interest (the nuclei). The number may be 2 or 3.
 	
	% clean up steps enter 1 if you are satisfied and 0 if you are not.
	
	% minimal area: minimum number of pixels comprising a shape. If a shape has pixels less than this min, it will not be included in the output.
		e.g.: 10
	% dilation_shape is the shape you are used to dilate. 
	
In the ori_file.xml, the number of shapes and the coordinates can be seen after the program is done.

2 b). for only two distinct colors 

	imagename='circle_fill.jpg' or imagename='shapes2.jpg' //your black-and-white image name
	
	cordcut=Run_ALL(iamgename, 2,'file.txt', 10);
	
	% choose the minimal area:choose 10 or whatever you like.
	
The output x-y coordinates is in file.xml.
