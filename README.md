# River GCS Toolkit

This repository contains a set of tools for processing geospatial data. While some tools may have a variety of fluvial geomorphological applications, an emphasis is placed on analysis of [geomorphic covariance structures](http://pasternack.ucdavis.edu/research/projects/geomorphic-covariance-structures/) (GCS). To launch the entire toolkit, run `master.py`. Otherwise, there is an individual GUI included for each tool.

Currently, these tools can analyze the GCS's between cross-sectional flow widths (W), velocities (V), and detrended bed elevations (Z), though a GCS could be made for any pair of morphological parameters.
## Prerequisites

-ArcGIS + Python 2.7; i.e. the `arcpy` python package and ArcGIS license.

-LASTools: functionality used for LiDAR Data Processing tool. LASTools can be downloaded [here](https://rapidlasso.com/lastools/).

## Getting Started

After installing any missing prerequisites, download this repository and launch `master.py`. An explanation of each tool is provided below.


## Tools

### 1. LiDAR Data Processing

This tool uses the open-source software LASTools to reclassify and separate LiDAR points.

### 2. Create Centerline

This tool uses a least-cost path approach on a DEM input to produce a river channel centerline shapefile.

### 3. Create Station Lines

This tool produces a set of evenly spaced cross-sectional lines perpendicular to a centerline input.

### 4. Detrend DEM

This tool creates a detrended DEM via linear or quadratic regression on the river channel's longitudinal profile. Detrending can be split into multiple regressions by reach.

### 5. Extract GCS Data

This tool produces a table of GCS data for each discharge.

### 6. Classify Landforms

This tool uses the GCS data in conjunction with a classification scheme to assign a GCS-based landform classification to each cross-section.

### 7. GCS Analysis

This tool uses the classified GCS data to run several analyses and aggregate data into tables and plots.