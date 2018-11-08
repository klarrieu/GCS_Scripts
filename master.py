'''
This combines the scripts for every GCS data processing/analysis step into one compact GUI.

Uses a GCSHandler class to store input data as class attributes, store processing/analysis functions as class methods.

Use an excel notebook to store input data for things like:
-landform classification thresholds?
-reach breaks?
-

The central GUI employs a tab for each processing step.

tabs for each of the following utilities:

-LiDAR reprocessing
-create centerline
-DEM detrending
-extract channel geometry data
-classify landforms by GCS
-GCS analysis
-hypsograph analysis
-landform model analysis (landform stratified velocity)
-reach-scale aggregate metrics
'''