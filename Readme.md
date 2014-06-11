#Coelacanth
###Codes for EMCCD and Lucky-Imaging Analysis
![logo](http://upload.wikimedia.org/wikipedia/commons/c/ce/Coelacanth-bgiu.png)

*There's life in the old fish yet.*

#### What is this?
*Coelacanth* is an ongoing effort to refactor and document a collection of C++
codes I wrote for lucky-imaging data-reduction and analysis, as described in 
my [PhD thesis]. 
The hope is that it might be useful to others: either in whole, in part, or as
a cautionary lesson in how we used to do things.

#### What's the fish got to do with it?
<em>Co</em>des for <em>E</em>MCCD and <em>l</em>ucky-imaging <em>a</em>nalysis == *COELA* , 
but 'coelacanth' is more memorable. 
That, and [C++03] is a bit of a dinosaur.

#### So what's buried in there?
- A framework for performing complex multi-coordinate-system pixel arithmetic on 
  two-dimensional pixel arrays.
- A minimal, C++-style implementation of the [FITS] image-file format, 
  complete with bit-packing decompression.
- A basic but reasonably performant implementation of the [Keys81]  algorithm for 
  [bicubic interpolation].
- Routines for bright-speckle location via cross-correlation with a carefully 
  tuned PSF core model.
- A simple (translation-only) but reasonably performant implementation of the 
  [Drizzle] algorithm.
- High-performance Monte-Carlo simulation code for simulating EMCCD data from 
  light-intensity maps (built on [UNURAN]).
- Fitting routines for extracting [EMCCD] parameters from raw data.

- A full (multi-threaded) lucky-imaging pipeline, including: 
    - Asynchronous file-buffering and decompression
    - Calibration routines for constructing per-pixel histograms from dark-frames 
      data
    - Image debias / flat-fielding routines built upon the EMCCD calibration
    - Photon-thresholding for faint targets
    - Various routines for multi-EMCCD data-synchronisation and mosaic building.
    

- A set of decomposable, [DRY]-principle-adherent CMake build scripts.
- Four years of sweat and C++ related tears.

#### This sounds useful. How do I get started?
Unfortunately there's no documentation yet, and progress is slow as this 
is effectively a side-project.
However, if you're interested in digging in, do [get in touch] and I'll help
out if I can.



[C++03]: http://en.wikipedia.org/wiki/C%2B%2B03
[Drizzle]: http://adsabs.harvard.edu/abs/2002PASP..114..144F
[bicubic interpolation]: http://en.wikipedia.org/wiki/Bicubic_interpolation
[DRY]: http://en.wikipedia.org/wiki/Don%27t_repeat_yourself
[EMCCD]: http://en.wikipedia.org/wiki/Charge-coupled_device#Electron-multiplying_CCD
[FITS]: http://en.wikipedia.org/wiki/FITS
[get in touch]: http://timstaley.co.uk/
[Keys81]: http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=1163711
[UNURAN]: http://statmath.wu.ac.at/unuran/
[PhD thesis]: http://uk.arxiv.org/abs/1404.5907
