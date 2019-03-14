simpleCT: introduction
------------

simpleCT is a Matlab toolbox that allows to read and apply basic image processing steps to Computed Tomography (CT) datasets.

simpleCT comes as a folder that you should add (with all sub folders) to your Matlab path. This contains one main class (sctdata) for handling a CT dataset plus functions that implement specific processes.

Some of the simpleCT methods require external dependencies that you should add to your path before you start working with the toolbox.

Below, the basic functionalities and some example calls of simpleCT are presented.

Contents
--------

This example shows how to:

*   load CT volume data using the sctdata class
*   display the data histogram
*   set a property, calibrate the data to BMD
*   display cross-sections through the volume and print them as PNG images
*   segment the bone tissue
*   crop the volume
*   crop the volume using a bounding box of a binary mask

Keyworsds: CT, image processing, bone structure

0.1 load .mhd image data with the sctdata class
-----------------------------------------------

sctdata is the main class of the simpleCT toolbox. it allows to store, display and perform basic operation on CT volume data. To start working with it, the first thing that we need to do is to initialize one instance of the class. We do it calling:

`HRdata = sctdata;`

If you omit the semicolon in the statement above you can immediately see a list of the properties of an sctdata object.

sctdata has a load method that allows you to open and work with several biomedical image data formats. Load sample CT volume with:

```matlab
HRdata.load
```
or:
```matlab
HRdata.load('mouse\_femur2.mhd');
```

Specific input/output methods are collected in the outside of the class inside the io folder of the simplect toolbox. For each supported data type, a class provides read and write methods. These are not specific for sctdata, meaning that they can be used to load several data formats to the workspace. The load method of sctdata uses them to read the data and all information belonging to it in an sctdata object.

Supported file formats are:

*   DICOM data (\*.dcm,\*.DCM)
*   MetaImage data (\*.mhd,\*.MHD) (\*.mha,\*.MHA)
*   (Scanco) ISQ data (\*.isq,\*.ISQ)
*   (Scanco) AIM data (\*.aim,\*.AIM)
*   VolFloat32 data (\*.vol)
*   MATLAB matrices (\*.mat,\*.MAT)
*   image stacks (\*.jpg;\*.gif;\*.png;\*.bmp;\*.tif)

0.2 display data histogram
--------------------------

```matlab
HRdata.showhist;
```
0.3 calibrate to BMD
--------------------

The HU-BMD calibration is assumed to follow a linear relation of the type:

*BMD = HU x (density\_slope / Mu\_scaling) + density\_offset*

The three parameters:

*   mu\_scaling
*   slope\_dens
*   offset\_dens

define the linear conversion from HU to BMD. In files such as DICOM and ISQ this information is usually stored in the file header. On the contrary, MHD files do not allow to store scanner settings. We need to retrieve these parameters from the log book of the scanner and set the specific properties with:
```matlab
HRdata.setSlope_dens(9.7);
HRdata.setOffset_dens(-30);
HRdata.setMu_scaling(1.4);
```
Apply the calibration and display the histogram after BMD conversion with:

```matlab
HRdata.BMDcalibrate;
HRdata.showhist;
```

0.4 display volume cross-sections
---------------------------------

```matlab
HRdata.showslices;
HRdata.showslices(28, 18, 30);
```
You can save them as PNGs with:
```matlab
HRdata.writeMidplanes('test_mouse_femur.png');
```
0.5 thresholding
----------------

The function threshbone included in the simplect toolbox contains several methods for the segmentation of CT images. The easiest way is to apply a fixed global threshold (for example 40% of the maximumg Grey Value). You can calculate and set a mask of the bone tissue with:

```matlab
GV040_mask = threshbone(HRdata.data, 'fixed', 0.4*max(HRdata.data(:)));
HRdata.setMask(GV040_mask);
HRdata.showslices;
```
Note that sctdata allows you to store only one binary mask associated with each dataset.

The Matlab function multithresh calculates one (or multiple) image threshold levels using the Otsu method. This is also implemented in the threshbone function of simplect:
```matlab
HRdata.setMask(threshbone(HRdata.data, 'Otsu'));
```
You can have a look at all methods implemented in threshbone with:
```matlab
help threshbone
```
0.6 cropping
------------

sctdata contains several methods for cropping a volume of interest. To manually select the crop region (very basic.. only crop of X- and Y-dimensions allowed..):
```matlab
HRdata_cropped1 = HRdata.crop;
```
Note that by calling crop with an additional output you can create a separate version of the data.

To crop specifying origin and size of the VOI (in pixels) use the crop\_pixel method
```matlab
HRdata_cropped2 = HRdata.crop_pixel(18, 6, 1, 39, [], 50);
```
Other crop methods are crop\_pixel\_fromcenter and crop\_around\_pixel

0.7 bounding box of the bone tissue
-----------------------------------

The function bbox allows to calculate the limits of a bounding box of abinary mask:
```matlab
[row0, rowd, col0, cold, slice0, sliced] = bbox(HRdata.mask, 1, [], true);
```
The bbox procedure is implemented in the crop\_pixel for automatic cropping using a bounding box of the mask that is stored in your sctdata instance:
```matlab
HRdata_cropped3 = HRdata.crop_pixel('bbox');
```
