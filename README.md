cuFSDAF
========
Version 1.0

Overview
========
Historical remote sensed images have either lower spatial resolution or lower temporal resolution limited by hardware technology and atmospheric conditions. Spatiotemporal data fusion algorithms provide feasible ways to produce historical images with both high spatial resolution and high temporal resolution. The Flexible Spatiotemporal DAta Fusion (FSDAF) developed by Zhu et al. (2016) is a  spatiotemporal data fusion algorithm suitable for heterogeneous landscapes and captures land cover changes. However, the extensive computational complexity of FSDAF prevents it from being applied in mass production. Besides, the domain-decomposition strategy of FSDAF causes accuracy loss at the edges of sub-domains due to insufficient consideration of edge effects.

To overcome the computational barrier and modify the accuracy loss, we designed and implemented an enhanced FSDAF algorithm parallelized using GPUs, named cuFSDAF.  Compared with FSDAF, cuFSDAF greatly reduced the computing time while maintaining the accuracy. Experiments showed that cuFSDAF achieved a speedup of 27.8 using a Nvidia Geforce 1050 (Notebook) GPU compared with a sequential C++ FSDAF algorithm running on an Intel Core I5-7400 CPU.

Key features of cuFSDAF:
========
+ Requires only one pair of Coarse-Fine images
+ Decomposes input images adaptively when the size of images exceeds the GPU’s memory
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
  5. Click menu Build -> Build Solution  
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
