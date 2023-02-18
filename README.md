# DIP-HDR
Creation of HDR image Project for the course of Digitial Image Processing in Aristotle University of Thessaloniki\
(Dept. of Electrical and Computer Engineering) based on the paper "Recovering high dynamic range radiance maps from photographs",\
by Paul E. Debevec and Jitendra Malik, written in 2008, which can be found [here][paper]

# Brief Background on HDR images

In order to achieve a numerical approximation of the [irradiance][irr] $E$ of a scene depicted by a set of digital images, taken with different exposure times
and based on the [paper][paper] of Debevec et al., it holds that: 
$$g(Z_{i,j}^k) = \ln{E_{i,j}} + \ln{\Delta t_k}$$
where: $Z_{i,j}^k \in [Z_{min},Z_{max}]$ is the value of the pixel $(i,j)$ of the $k$-th image,\
$E_{i,j}$, is the irradiance for this pixel (assumed constant over the multiple shots of the scene),\
$\Delta t_k$, is the exposure time for the $k$-th shot\
and $f=e^{g}$, is a nonlinear function, assumed to be smooth and monotonic.



[paper]: <https://doi.org/10.1145/1401132.1401174> "Recovering high dynamic range radiance maps from photographs, by Debevec et al., 2008"
[irr]: <https://en.wikipedia.org/wiki/Irradiance> "Wikipedia: Irradiance"
