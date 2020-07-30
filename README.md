cuFSDAF
========
Version 1.0

Overview
========
Spatiotemporal data fusion is a cost-effective way to produce remote sensing images with high spatial and temporal resolutions using multi-source images. Using spectral unmixing analysis and thin plate spline (TPS) interpolation, the Flexible Spatiotemporal DAta Fusion (FSDAF) algorithm is suitable for heterogeneous landscapes and capable of capturing abrupt landcover changes. However, the extensive computational complexity of FSDAF prevents its use in large-scale applications and mass production. In addition, the domain decomposition strategy of FSDAF causes accuracy loss at the edges of sub-domains due to the insufficient consideration of edge effects. In this study, an enhanced FSDAF (cuFSDAF) is proposed to address these problems, and includes three main improvements: (1) The TPS interpolator is replaced by a modified inverse distance weighted interpolator to reduce computational complexity. (2) The algorithm is parallelized based on Compute Unified Device Architecture (CUDA), a widely used parallel computing framework for graphics processing units (GPUs). (3) An adaptive domain decomposition method is proposed to improve the fusion accuracy at the edges of sub-domains, and to enable GPUs with varying computing capacities to deal with datasets of any size. Experiments showed that, while maintaining accuracy, cuFSDAF reduced computing time significantly and achieved speed-ups of 115.2–133.9 over the IDL-implemented FSDAF, and speed-ups of 75.5–81.8 over the C++-implemented FSDAF. cuFSDAF is capable of efficiently producing fused images with both high spatial and temporal resolutions to support applications for large-scale and long-term land surface dynamics. 


Key features of cuFSDAF:
========
+ Requires only one pair of Coarse-Fine images
+ Decomposes input images adaptively when the size of images exceeds the GPU’s memory
+ Proposes a modified IDW interpolator for cuFSDAF
+ Improves fusion accuracy at the edges of sub-domains
+ Supports a wide range of CUDA-enabled GPUs (https://developer.nvidia.com/cuda-gpus)
+ Supports a wide range of image formats (see http://gdal.org/formats_list.html)
+ Supports both Windows and Linux/Unix operating systems

References
========
+ Zhu, X. et al., 2016. A flexible spatiotemporal method for fusing satellite images with different resolutions. Remote Sensing of Environment, 172: 165-177.  

To Cite cuFSDAF in Publications
========
+ A paper describing cuFSDAF will be submitted to a scientific journal for publication soon
+ For now, you may just cite the URL of the source codes of cuFSDAF (https://github.com/HPSCIL/cuFSDAF) in your publications  
+ You may contact the e-mail <ghcug14@cug.edu.cn> if you have further questions.

Compilation
========
+ Requirements:
  -	A computer with a CUDA-enabled GPU (https://developer.nvidia.com/cuda-gpus)
  -	A C/C++ compiler (e.g., Microsoft Visual Studio for Windows, and gcc/g++ for Linux/Unix) installed and tested
  -	Nvidia CUDA Toolkit (https://developer.nvidia.com/cuda-downloads) installed and tested
  -	Geospatial Data Abstraction Library (GDAL, https://gdal.org) installed and tested
  -	ALGLIB Library (https://www.alglib.net/) installed and tested

+ For the Windows operating system (using MS Visual Studio 2017 as an example)
  1. Create a project that uses the CUDA runtime
  2. Open all the source codes in VS 2017
  3. Click menu Project -> Properties -> VC++ Directories -> Include Directories, and add the “include” directory of GDAL and alglib (e.g., C:\GDAL\include\)
  4. Click menu Project -> Properties -> VC++ Directories -> Lib Directories, and add the “lib” directory of GDAL and alglib (e.g., C:\GDAL\lib\)
  5. Click menu Project -> Properties -> Link -> Input, and add .lib files of GDAL and alglib (e.g., C:\GDAL\lib\)
  6. Click menu Build -> Build Solution  
  Once successfully compiled, an executable file, cuFSDAF.exe, is created.
+ For the Linux/Unix operating system (using the CUDA compiler --- nvcc)  
  In a Linux/Unix terminal, type in: 
  ```
  - $ cd /the-directory-of-your-source-codes/
  - $ nvcc -std=c++11 -o cuFSDAF main.cpp FSDAF.cpp kernel.cu -lgdal -lalglib
  ```
  If your compiler fail to find CUDA/GDAL/ALGLIB, you may add their include path and library path (e.g. /path-of-gdal).
  ```
  - $ nvcc -std=c++11 -o cuFSDAF main.cpp FSDAF.cpp kernel.cu -lgdal -lalglib -L /path-of-gdal/lib -I /path-of-gdal/include
  ```
  Once successfully compiled, an executable file, cuFSDAF, is created.  
  Note: You may compile alglib firstly if it has not been compiled. When compiling the alglib library, name the library file as "libalglib.a". If not, remember to modify "-lalglib" in the above command.


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
> IN_F1_NAME = F:\\Datasets\\Daxing_Test\\L-2019-6-14.tif  
>
>\# The input coarse image at t1  
> IN_C1_NAME = F:\\Datasets\\Daxing_Test\\M-2019-6-14.tif  
>
>\# The input coarse image at t2  
> IN_C2_NAME = F:\\Datasets\\Daxing_Test\\M-2019-8-17.tif  
>
>\# The classified image for the fine image at t1  
> IN_F1_CLASS_NAME = F:\\Datasets\\Daxing_Test\\class  
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
>\# Search radius for IDW interpolator, recommend at least 2*SCALE_FACTOR  
> IDW_SEARCH_RADIUS = 32  
>
>\# Power to calculate the weight of known points in IDW, if 2, Weight = 1/Distance^2  
> IDW_POWER = 2  
>
>cuFSDAF_PARAMETER_END  

+ The program runs as a command line. You may use the Command (i.e., cmd) in Windows, or a terminal in Linux/Unix. 
   - For the Windows version:    
   $ cuFSDAF.exe parameters.txt 
   - For the Linux/Unix version:   
   $ ./cuFSDAF parameters.txt 

+ Note: The computational performance of cuFSDAF largely depends on the GPU. The more powerful is the GPU, the better performance. 


