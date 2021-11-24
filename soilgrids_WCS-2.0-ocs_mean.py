#!/usr/bin/env python
# coding: utf-8

# Working with WCS 2.0
# ====================
# 
# The WCS standard underwent a major update in 2017 that introduced important improvements. Among these is a more abstract treatment of data axes, thinking especially of space-time datasets and data cubes. The friendly [OWSLib](https://geopython.github.io/OWSLib/) package has been updated accordingly and can be easily be used with the version 2.0 of the standard.
# 
# Connection set up and service inspection functions as before:

# In[1]:


from owslib.wcs import WebCoverageService

wcs = WebCoverageService('http://maps.isric.org/mapserv?map=/map/ocs.map',
                         version='2.0.1')


# With the 2.0 version, supported are announced in slightly different way, let us thus inspect these in first place:

# In[2]:


fileout = './data/workfiles/soilgrids_ocs_mean_iowa_approx.tif'
cov_id = 'ocs_0-30cm_mean'
ph_0_5 = wcs.contents[cov_id]
ph_0_5.supportedFormats


# Perhaps the most noteciable change is the way segmentation is expressed. The infamous bounding boxes have been dropped at last and replaced by axis specific segmentation. The new argument `subsets` takes in a list of tuples, one tuple per axis. Each tuple is a triplet: (*axis identifier*, *minimum*, *maximum*). Here is an example with the same extent for Senegal used in the previous notebook:

# In[3]:


subsets = [('X', -11000000, -10000000), ('Y', 4400000, 5000000)]


# CRSs are also expressed in a different way in version 2.0, now referring to the [opengis.net](http://www.opengis.net) registry. Again the surrogate EPSG code 152160 is used, as the Petroleum folk remain unimpressed by the Homolosine projection:

# In[4]:


crs = "http://www.opengis.net/def/crs/EPSG/0/152160"


# All information is now in place to issue a `GetCoverage` request and save the result to disk: 

# In[5]:


response = wcs.getCoverage(
    identifier=[cov_id], 
    crs=crs,
    subsets=subsets, 
    resx=250, resy=250, 
    format=ph_0_5.supportedFormats[0])


# In[6]:


with open(fileout, 'wb') as file:
    file.write(response.read())


# [Index](index.ipynb) | [Previous](02-WCS-getExtent.ipynb)
