# Cortical Histology Analysis Toolbox (CHAT)

## What is it?
CHAT is a set of Matlab functions for measuring neuron density, neuron size, and cortical layer thickness 
from histology images of cerebral cortex, as described in:

[**[1]**] Abbass, M., Trought, K., Long, D., Semechko, A., Wong, A.H.C., 2018. Automated immunohistochemical 
method to analyze large areas of the human cortex. Journal of Neuroscience Methods. Vol. 294, pp. 81-90.

[**[2]**] Trought, K., 2017. Automated immunohistochemical analysis of the orbitofrontal cortex in patients 
with schizophrenia, bipolar disorder and major depressive disorder. MSc Thesis. University of Toronto.

[**[3]**] Abbass, M., 2014. Cytoarchitecture of the anterior cingulate cortex in patients with schizophrenia,
 bipolar disorder and major depression. MSc Thesis. University of Toronto.

## Quick Demo and Description of CHAT Functionality
> **IMPORTANT: CHAT runs out-of-the-box and requires no compiling. However, it is also heavily dependent on the [Image Processing Toolbox] for various image processing tasks and will not run without this toolbox.**
1. Add CHAT repository to your Matlab path.
2. Navigate to `...\Cortical-Histology-Analysis-Toolbox\Sample Data Folder\Subfolder\` and inspect contents of the sample image data folders contained therein: `sample_01\`, `sample_02\`, and `sample_03\`. In each of these folders, you will see 5 **_pre-processed_** images:
    * *DAPI*
    * *CUX2*
    * *ZNF312*
    * *pia-wm*
    * *polylines*
    
    The first three images are cell segmentations corresponding to the DAPI, CUX, and ZNF stains. The fourth is the image of the "outer" and "inner" boundaries of the cortex (i.e., pia and white matter, resp.). The last image contains four (or five) open contours partitioning the cortex "longitudinally" into five (or six) cortical layers. Note that the images in a given folder have exactly the same dimensions; because they were derived from the same histology slice. **_Methods to pre-process raw histology images prior to data extraction are not (yet) part of CHAT; see [[1]] for details._**

3. Enter `CHAT_demo` into Matlab command prompt and wait until Matlab extracts relevant data from `sample_01\`, `sample_02\`, and `sample_03\` (5-6 min). When Matlab finishes processing the images in a given data folder, it will deposit into the said folder two files: `cortex_boundaries.tif` and `cortex_data.txt`. The former should be used for visual inspection and confirmation that the implemented algorithm correctly identified cortical boundaries. The latter contains records of:
    * average white matter to pia distance
    * average width of the cortex
    * average width of cortical layers
    * area of cortical layers
    * number of cells per cortical layer
    * relative distances (range 0 to 1) of cells to pia and corresponding cell sizes

4. Use `CHAT_ProcessImageData.m` function to process your own data. Enter `help CHAT_ProcessImageData` into Matlab command window for additional details regarding its usage.

## License
[MIT] © 2019 Anton Semechko 
a.semechko@gmail.com

[1]: https://doi.org/10.1016/j.jneumeth.2017.10.024
[2]: http://hdl.handle.net/1807/79362
[3]: http://hdl.handle.net/1807/71965
[Image Processing Toolbox]: https://www.mathworks.com/products/image.html
[MIT]: https://github.com/AntonSemechko/Cortical-Histology-Analysis-Toolbox/blob/master/LICENSE.md
