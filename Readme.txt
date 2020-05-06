

***************************************************************************
* authors: Nelly Pustelnik, Barbara Pascal, Patrice Abry                  *
* institution: laboratoire de Physique de l'ENS de Lyon                   *
* date: May 2020                                                          *
* License CeCILL-B                                                        *
***************************************************************************
*********************************************************
* RECOMMENDATIONS:                                   	*
* This toolbox is designed to work with Matlab 2018.b   *
*********************************************************

------------------------------------------------------------------------------------------------------------------------
DESCRIPTION:
Numerous fields of nonlinear physics, very different in nature, produce signals and images, that share the common feature of being essentially constituted of piecewise homogeneous phases.
Analyzing signals and images from corresponding experiments to construct relevant physical interpretations thus often requires detecting such phases and estimating accurately  their characteristics (borders, feature differences,\ldots). 
However, situations of physical relevance often comes with low to very low signal to noise precluding the standard use of classical linear filtering for analysis and denoising and thus calling for the design of advanced nonlinear signal/image filtering techniques. 
Additionally, when dealing with experimental physics signals/images, a second limitation is the large amount of data that need to be analyzed to yield accurate and relevant conclusions (e.g., in producing a phase diagram or in analyzing  video frames of large size images) requiring the design of fast algorithms.
The present work proposes a unified signal/image nonlinear filtering procedure, with fast algorithms and an data-driven automated hyperparameter tuning, based on proximal operators and Stein unbiased estimator principles.
The interest and potential of these tools are illustrated at work on low-confinement solid friction signals and porous media multiphase flows. 


------------------------------------------------------------------------------------------------------------------------
SPECIFICATIONS for using TV_SURE toolbox:

The main function is "nonlinear_filtering.m".

Demo files are provided in:
- "main_piecewise_constant_multivariate.m" : multivariate piecewise constant denoising
- "main_piecewise_linear_univariate.m" : univariate linear constant denoising
- "main_TV2D.m" : 2D TV denoising
- "main_steinmethods.m" allowing to perform piecewise constant and piecewise linear denoising.

This toolbox makes use of '' GRANSO-master '' (required) to perform BFGS quasi-Newton minimization of SURE. The toolbox can be dowloaded from http://www.timmitchell.com/software/GRANSO/.
. For convenience purposes, copy of the latter toolbox is provided in the subfolders.

------------------------------------------------------------------------------------------------------------------------
RELATED PUBLICATION:

# B. Pascal, N. Pustelnik, P. Abry, J.-C. Geminard, V. Vidal. 
Parameter-free and fast nonlinear piecewise filtering,
Annals of telecom, 2020

------------------------------------------------------------------------------------------------------------------------