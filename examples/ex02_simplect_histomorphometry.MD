simpleCT: functions for bone histomorphometry
-----------

This example illustrates the use of some of the functions included in simpleCT for the quantitative analysis of bone's density and 3D histomorphometry.

Contents
--------

This example shows how to:

*   load portion of High Resolution pQCT (HR-pQCT) data of a human femur neck sample
*   obtain masks of different bone compartments
*   calculate bone thickness and porosity

Suggested readings
------------------

Book chapters on bone imaging and bone histomorphometry

*   [Cowin, Stephen C. Bone mechanics handbook. CRC press, 2001. - CHAPTER 9 - Imaging of Bone, Peter Rüegsegger](https://www.crcpress.com/Bone-Mechanics-Handbook/Cowin/p/book/9780849391170)
*   Cowin, Stephen C. Bone mechanics handbook. CRC press, 2001. - CHAPTER 14 - Quantification of Cancellous Bone Architecture - Anders Odgaard

Nomenclature for bone histomorphometry:

*   [Dempster DW, Compston JE, Drezner MK, Glorieux FH, Kanis JA, Malluche H, Meunier PJ, Ott SM, Recker RR, Parfitt AM. Standardized nomenclature, symbols, and units for bone histomorphometry: a 2012 update of the report of the ASBMR Histomorphometry Nomenclature Committee. Journal of bone and mineral research. 2013 Jan 1;28(1):2-17.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3672237/)

Kaywords: CT, bone histomorphometry, BMD, HR-pQCT

0.1 load .mhd image data with sctdata
-------------------------------------
```matlab
HRdata = sctdata;
HRdata.load;
HRdata.showslices;
```
The following operations have already been performed on the data that you just loaded:

*   GV of air bubbles voxels set to the GV of water
*   3-D gaussian filtering (removes high frequency noise)
*   downsampling with factor 2 (reduces image size by a factor of 8
*   rotate (allign the Z-axis of the image with the axis of the femur neck)
*   crop (extracts a 7-mm thick region around the neck section with minimum cross-sectional area)

0.2 thresholding
----------------

The method for automatic separation of the cortical compartment on HR-pQCT images implemented in the Scanco IPL is described here:

*   [Burghardt, A. J., Buie, H. R., Laib, A., Majumdar, S., & Boyd, S. K. (2010). Reproducibility of Direct Quantitative Measures of Cortical Bone Micro-architecture of the Distal Radius and Tibia by HR-pQCT. Bone, 47(3), 519–528.](http://doi.org/10.1016/j.bone.2010.05.034)

To start, a fixed global threshold is calculated as the 40% of the maximum Gray Value within the scanned volume.
```matlab
GV040_mask = threshbone(HRdata.data, 'fixed', 0.4*max(HRdata.data(:)));
```
0.3 compartment masks
---------------------

0.3.1 periosteum mask
---------------------

The function periosteummask included in simplect allows to track the periosteum boundary (ando obtain a mask of logical) from a binary mask of the bone tissue:
```matlab
peri_mask = periosteummask(GV040_mask);
HRdata.setMask(peri\_mask);
HRdata.showslices;
% HRdata.writeMidplanes('test_neck_peri.png');
```
0.3.2 cortex mask
-----------------

The function corticalmask implements the procedure of (Burghardt et al., 2010).
```matlab
cort_mask = corticalmask(GV040_mask, HRdata.voxelsize(1));
HRdata.setMask(cort_mask);
HRdata.showslices;
```
The procedure can be further adapted to account for the different (generally thinner) structure of cortical bone in the femoral neck (the procedure was developed for tibial cortical bone).
```matlab
help corticalmask
[cort_mask] = corticalmask(GV040_mask, HRdata.voxelsize(1), 1, 1, 4, 0.25, 1.5, peri_mask, 1);
HRdata.setMask(cort_mask);
HRdata.showslices;
```
0.3.3 trabecular mask
---------------------

Can be obtained combining the two masks calculated above..
```matlab
trab_mask = peri_mask & ~cort_mask;
HRdata.setMask(trab_mask);
HRdata.showslices;
```
At the interface between cortical and trabecular compartments exists a transition zone that we might want to exclude for the histomorphometric characterization of the trabecular core. We therefore erode (and reduce the size) 1-mm of the trabecular mask. The radius of the structuring element is:
```matlab
erode_trab_mm = 1; % [mm]
erode_trab_voxels = round(erode_trab_mm / HRdata.voxelsize(1)); % [voxels]
```
Create structuring element for erosion
```matlab
se = strel('disk', erode_trab_voxels); % 2-D disk
se3d = strel3d(2*erode_trab_voxels); % 3-D ball
```
Perform erosion and display the result:
```matlab
trab_mask_eroded = imerode(trab_mask, se3d);
HRdata.setMask(trab_mask_eroded);
HRdata.showslices;
```
0.4 cortical bone thickness (Ct.Th)
-----------------------------------

Calculate point cloud of the periosteum contour:
```matlab
peri_con = periosteumcontour(peri_mask);
```
For the endosteum we need to avoid small structures non-connected to the main bone region:
```matlab
endo_con = endosteumcontour(removenonconnectedstruts(endo_mask, 'regionprops2'));
```
The function CtTh calculates the cortical bone thickness from the histogram of the miminum distances separating periosteum and endosteum contours:
```matlab
[ctth, ctth_array, ctth_bins] = CtTh(endo, peri, 'pdist2', 1);
```
0.5 cortical bone porosity (Ct.Po)
----------------------------------

Threshold-based:
```matlab
CtPo0 = 100 * (1-((GV040_mask & cort_mask) / cort_mask));
```
As in Burghardt et al., 2010 (also tthreshold-based).
```matlab
CtPo1 = burghardt_CtPo(GV040_mask, cort_mask, HRdata.voxelsize(1));
```
