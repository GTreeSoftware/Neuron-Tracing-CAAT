# Automated neuron tracing using content-aware adaptive voxel scooping on CNN predicted probability map
 A method for accurate neuron tracing using content-aware adaptive voxel scooping on convolutional neural network (CNN) predicted probability map.
# Automated neuron tracing using content-aware adaptive voxel scooping on CNN predicted probability map

## Overview:

Neuron tracing plays an important role in the understanding of brain organization and function. Current methods often had trouble in tracing the complex tree-like  structures and broken parts of neurite from noisy background. We propose a method for accurate neuron tracing using content-aware adaptive voxel scooping on convolutional neural network (CNN) predicted probability map. First, a 3D residual CNN was applied as preprocessing to predict the object probability and suppress high noise. Then, an adaptive voxel scooping method was presented for successive neurite tracing on the probability map, based on the internal content properties (distance, connectivity, and probability continuity along direction) of neurite. Last, the neuron tree graph was built using the length first criterion. The proposed method was evaluated on the public BigNeuron datasets and fluorescence micro-optical sectioning tomography (fMOST) datasets, and outperformed current state-of-art methods on images with neurites had broken parts and complex structures.

## System Requirements

### Hardware Requirements:

The deep learning algorithm requires enough RAM and GPU to support the calculation. For optimal performace, we recommenda computer with the following specs:
**RAM**: 16+GB
**CPU**: Intel i7 or better
**GPU**:  1080Ti or better

### Software Requirements:

**pycharm**:  used for deep learning (DL) code editing. DL algorithm is used for  neurite probability prediction
**matlab**:  version 2017 or better; the neuron tracing algortihm (referred as CAAT) is built using matlab language . 

### Environment Requirements:

**CUDA**: cuda 9.0
**cudnn**: cudnn 7
**Python**: 3.6
**pytorch**:0.4.1 
**visdom**:0.1.8.5
**Numpy**: 1.14.5
**tifffile**: 0.15.1
**Scikit-image**:0.13.1

### CNN code for neurite prediction:

The code and model that used for neurite prediction can be accessed via: https://github.com/GTreeSoftware/DB-Enhance. Detailed explanations can be accessed via the "readme" file in the website. 

### CAAT code for neurite tracing:

**Test1.m**: The main function used for neurite tracing, which is applied on the predicted neurite probability map by CNN.  This function support 3D image data with the format of 'tif'. You can just change the 'image_name' of your test data, and run the code. You can get the traced result of the neurite and visualize the traced result in 3D.

**Read_3D_Tiff.m**: This function is used to read 3D image with tif format:---i.e. predicted neurite probability.

**CBK_Threshold_Region.m**: This function is used to estimate the initial threshold of the predicted neurite probability.

**Ada_Voxel_Scooping1.m**: This function is used to for neurite tracing.

**Build_Tracing_Graph.m**: This function is used to build the neuron tree graph using the length first criterion. You can prun the short branches based on the branch skeleton points number (parameter line_threshold). It is set 6 as default. You can change it as you like.  

### Test Dataset:

We also include 7 testing images for testing  in the 'image' file under the 'CAAT_Image' file. 

The original image (uint8  format) are put in the 'image' file, and the corresponding CNN predicted neurite probability are put in the 'prob' file.  

Image 7 is a larger image with size 1000×1000 ×300. It contains multiple neurites with inhomogeneous intensities.
The datasets can be accessed via: https://github.com/GTreeSoftware/DB-Enhance/releases/tag/testdata1.
