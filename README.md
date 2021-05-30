cuFSDAF
========
Version 1.0

Overview
========

Spatiotemporal data fusion is a cost-effective way to produce remote sensing images with high spatial and temporal resolutions using multi-source images. Using spectral unmixing analysis and spatial interpolation, the Flexible Spatiotemporal DAta Fusion (FSDAF) algorithm is suitable for heterogeneous landscapes and capable of capturing abrupt land-cover changes. However, the extensive computational complexity of FSDAF prevents its use in large-scale applications and mass production. Besides, the domain decomposition strategy of FSDAF causes accuracy loss at the edges of sub-domains due to the insufficient consideration of edge effects. In this study, an enhanced FSDAF (cuFSDAF) is proposed to address these problems, and includes three main improvements: (1) The TPS interpolator is replaced by an accelerated inverse distance weighted interpolator to reduce computational complexity. (2) The algorithm is parallelized based on the Compute Unified Device Architecture (CUDA), a widely used parallel computing framework for graphics processing units (GPUs). (3) An adaptive domain decomposition method is proposed to improve the fusion accuracy at the edges of sub-domains, and to enable GPUs with varying computing capacities to deal with datasets of any size. Experiments showed while obtaining similar accuracies to FSDAF and an up-to-date deep-learning-based method, cuFSDAF reduced the computing time significantly and achieved speed-ups of 140.3–182.2 over the original FSDAF program. cuFSDAF is capable of efficiently producing fused images with both high spatial and temporal resolutions to support applications for large-scale and long-term land surface dynamics.   

Key features of cuFSDAF:
========
+ Requires only one pair of Coarse-Fine images
+ Decomposes input images adaptively when the size of images exceeds the GPU’s memory
+ Proposes an accelerated IDW interpolator for cuFSDAF
+ Improves fusion accuracy at the edges of sub-domains
+ Supports a wide range of CUDA-enabled GPUs (https://developer.nvidia.com/cuda-gpus)
+ Supports a wide range of image formats (see http://gdal.org/formats_list.html)
+ Supports both Windows and Linux/Unix operating systems

References
========
+ Zhu, X. et al., 2016. A flexible spatiotemporal method for fusing satellite images with different resolutions. Remote Sensing of Environment, 172: 165-177.  
+ Gao, H., Zhu, X., Guan, Q., Yang, X., Yao, Y., Zeng, W., Peng, X., 2021. cuFSDAF: An Enhanced Flexible Spatiotemporal Data Fusion Algorithm Parallelized Using Graphics Processing Units. IEEE Transactions on Geoscience and Remote Sensing. https://doi.org/10.1109/TGRS.2021.3080384  


To Cite cuFSDAF in Publications
========
+ Please cite the following reference:  
  Gao, H., Zhu, X., Guan, Q., Yang, X., Yao, Y., Zeng, W., Peng, X., 2021. cuFSDAF: An Enhanced Flexible Spatiotemporal Data Fusion Algorithm Parallelized Using Graphics Processing Units. IEEE Transactions on Geoscience and Remote Sensing. https://doi.org/10.1109/TGRS.2021.3080384  
+ If you use the datasets we used here, please cite the ref. [80] in the paper we mentioned above.
+ You may contact the e-mail <ghcug14@cug.edu.cn> if you have further questions about the usage of codes and datasets.
+ For any possible research collaboration, please contact Prof. Qingfeng Guan (<guanqf@cug.edu.cn>).  

Compilation
========
+ Requirements:
  -	A computer with a CUDA-enabled GPU (https://developer.nvidia.com/cuda-gpus)
  -	A C/C++ compiler (e.g., Microsoft Visual Studio for Windows, and gcc/g++ for Linux/Unix) installed and tested
  -	Nvidia CUDA Toolkit (https://developer.nvidia.com/cuda-downloads) installed and tested
  -	Geospatial Data Abstraction Library (GDAL, https://gdal.org) installed and tested
  -	ALGLIB Library (https://www.alglib.net/) installed and tested

+ For the Windows operating system (using MS Visual Studio 2017 as an example)
  1. Create a project that uses the CUDA runtime and named "cuFSDAF".
  2. Copy source codes of cuFSDAF to the path of your project, and open all source codes in VS 2017. If you had not compiled ALGLIB, you could open source codes of ALGLIB together. We prepared source codes of ALGLIB 3.13.0 for convenience, but you may also add another ALGLIB version.
  3. Click menu Project -> Properties -> VC++ Directories -> Include Directories, and add the “include” directory of GDAL and ALGLIB (e.g., C:\GDAL\include\, .\alglib\\).
  4. Click menu Project -> Properties -> VC++ Directories -> Lib Directories, and add the “lib” directory of GDAL (e.g., C:\GDAL\lib\). If you compiled ALGLIB, you may add your "lib" directory of ALGLIB as well (e.g., .\alglib\\).
  5. Click menu Project -> Properties -> Link -> Input, and add ".lib" files of GDAL (e.g., gdal_i.lib) and ALGLIB (if you compiled it).
  6. Click menu Build -> Build Solution.
  Once successfully compiled, an executable file, cuFSDAF.exe, is created.

+ For the Linux/Unix operating system (using the CUDA compiler --- nvcc) 
  1. Compile a lib file using the source code of ALGLIB (named "libalglib.a"), and copy it to the directory of your source codes. If not, add name of all ".cpp" files of ALGLIB to the command line in step 2 and delete "-lalglib".  
  2. In a Linux/Unix terminal, type in: 
  ```
  - $ cd /the-directory-of-your-source-codes/
  - $ nvcc -std=c++11 -o cuFSDAF main.cpp FSDAF.cpp kernel.cu -lgdal -lalglib -L /path-of-gdal/lib -I /path-of-gdal/include -I /include-path-of-alglib
  ```
  Once successfully compiled, an executable file, cuFSDAF, is created.  

Debug
========
1. "dataanalysis.h": No such file or directory  
  "dataanalysis.h" is a head file of ALGLIB, remember to add the include path of ALGLIB when compiling. 


Usage 
========
+ Before running the program, make sure that all Landsat and MODIS images have been pre-processed and co-registered. They must have:
  - the same spatial resolution (i.e., Landsat resolution --- 30m)
  - the same image size (i.e., numbers of rows and columns)
  - the same map projection
+ A text file must be manually created to specify input images, and other parameters for the cuFSDAF model.  
Example (# for comments):

>cuFSDAF_PARAMETER_START  
>
>\# The input fine image at t1  
> IN_F1_NAME = D:\\Datasets\\Daxing\\Landsat\\L-2018-12-4.tif  
>
>\# The input coarse image at t1  
> IN_C1_NAME = D:\\Datasets\\Daxing\\MODIS\\M-2018-12-4.tif  
>
>\# The input coarse image at t2  
> IN_C2_NAME = D:\\Datasets\\Daxing\\MODIS\\M-2018-10-1.tif  
>
>\# The classified image for the fine image at t1  
> IN_F1_CLASS_NAME = D:\\Datasets\\Daxing\\class_2018-12-4  
>
>\# Window's size for searching similar pixels
> W = 20
>
>\# Number of similar pixels for mitigating errors using neighborhood  
> NUM_SIMILAR_PIXEL = 20  
>
>\# Number of purest coarse pixels in each class for unmixing analysis  
> NUM_PURE = 100  
>
>\# Minimum of Digital Number (DN) value, used for correcting exterme values  
> DN_MIN = 0.0
>
>\# Maximum of Digital Number (DN) value, used for correcting exterme values  
> DN_MAX = 10000.0  
>
>\# The scale factor, it is integer=coarse resolution/fine resolution, e.g., 480/30=16    
> SCALE_FACTOR = 16
>
>\# The value of background pixels  
> BACKGROUND = 0  
>
>\# Which band with value = BACKGROUND indicating background pixels  
> BACKGROUND_BAND = 1  
>
>\# Search radius for IDW interpolator, recommend 2*SCALE_FACTOR  
> IDW_SEARCH_RADIUS = 32  
>
>\# Power to calculate the weight of known points in IDW, if 2, Weight = 1/Distance^2  
> IDW_POWER = 2  
>
>cuFSDAF_PARAMETER_END  

+ The program runs as a command line. You may use the Command (i.e., cmd) in Windows, or a terminal in Linux/Unix. 
   - For the Windows version:    
   $ .\cuFSDAF.exe Parameters.txt 
   - For the Linux/Unix version:   
   $ ./cuFSDAF Parameters.txt 

+ The fused image will be saved in the path of your input coarse image at t2, named "xxx_cuFSDAF.tif".  
+ Note: The computational performance of cuFSDAF largely depends on the GPU. The more powerful is the GPU, the better performance. 


