#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from setuptools import setup

setup(name='pyEpoch',
      version='0.1',
      description='Reconstructing GRNS from scRNA-seq data',
      url='http://github.com/pcahan1/PyEpoch/',
      author='Patrick Cahan',
      author_email='patrick.cahan@gmail.com',
      license='MIT',
      packages=['pyEpoch'],
      install_requires=[
          'pandas',
          'numpy',
          'scanpy',
          'pygam',
          'scipy',
          'sklearn',
          'scipy',
          'math',
          'igraph',
          'networkx',
          'matplotlib',
          'statsmodels',
          'researchpy',
          'bioinfokit',
          'sys',
          'skfda',
          'seaborn',
          'pyitlib'
      ],
      zip_safe=False)

