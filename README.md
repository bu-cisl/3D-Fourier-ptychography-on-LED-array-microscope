# 3D-Fourier-ptychography-on-LED-array-microscope

MATLAB implementation of 3D-Fourier-ptychography-on-LED-array-microscope. 

### Citation

If you find this project useful in your research, please consider citing our paper:

[**Tian, Lei, and Laura Waller. "3D intensity and phase imaging from light field measurements in an LED array microscope." optica 2.2 (2015): 104-111.**](https://www.osapublishing.org/optica/abstract.cfm?uri=optica-2-2-104)


### Abstract

Realizing high resolution across large volumes is challenging for 3D imaging techniques with high-speed acquisition. Here, we describe a new method for 3D intensity and phase recovery from 4D light field measurements, achieving enhanced resolution via Fourier ptychography. Starting from geometric optics light field refocusing, we incorporate phase retrieval and correct diffraction artifacts. Further, we incorporate dark-field images to achieve lateral resolution beyond the diffraction limit of the objective (5×
 larger NA) and axial resolution better than the depth of field, using a low-magnification objective with a large field of view. Our iterative reconstruction algorithm uses a multislice coherent model to estimate the 3D complex transmittance function of the sample at multiple depths, without any weak or single-scattering approximations. Data are captured by an LED array microscope with computational illumination, which enables rapid scanning of angles for fast acquisition. We demonstrate the method with thick biological samples in a modified commercial microscope, indicating the technique’s versatility for a wide range of applications.


### Running the code

MATLAB is required to run this code. Run MultiSlice_SuperRes.m

### Dataset

Dataset contains images captured from sequential LED illumination up to 0.41 NA, using a 4x objective (0.1 NA). In our paper, we achieved resolution corresponding to the sum of the two NAs~0.5. Our 3D FPM algorithms is based on multislice model that accounts for multiple scattering effects from 3D model. It combines 3D light field refocusing with FPM to achieve super resolution in both lateral (x,y) and axial (z) dimensions across a large field of view.

Dataset 1: Two stacked resolution targets separated by ~100um//
Dataset 2: Spiral Algea

### Result 

<p align="center">
  <img src="/images/result.png">
</p>
