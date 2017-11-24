# SVCM

主要实现了`Spatially Varying Coefficient Model for Neuroimaging Data with Jump Discontinuities`这篇文章里的算法，其中，不同之处在于估计变系数的时候采用了三个系数一起估计的方法，想法来自`MWPCR-Multiscale-Weighted-Principal-Component-Regression-for-High-dimensional-Prediction`这篇文章，`example.m`这个文件主要是实现了文章中的64×64×8的模拟，需要注意的是8G的内存的机器对于这个程序内存可能不够，如果想看文章效果，可以运行`example.m`中的20×20×10得例子，为了成功运行20×20×10的必须修改`GenerateData.m`中的`N`参数。
