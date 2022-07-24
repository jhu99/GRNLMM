GRNLMM
=========================================================
We development a new method, GRNLMM, which models single-cell expression data by a linear mixed model and uses the covariance matrix of random effect terms to characterize the correlation between genes. To overcome the influence of randomness of intercellular expression and improve the accuracy of the predicted GRNs, we use a known correlation matrix to reflect the relationship between cells and add a noise term to the model. Our results show that GRNLMM has advantages in accurately identifying the co-expressive relationships between genes and can explore genes and gene function modules that play an essential role in biological processes.

GRNLMM Installation
==================
.. toctree::
   :maxdepth: 2

   installation

GRNLMM Tutorial
==================
.. toctree::
   :maxdepth: 2
   
   tutorials/index.rst

GRNLMM API
==================  
.. toctree::
   :maxdepth: 2

   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
