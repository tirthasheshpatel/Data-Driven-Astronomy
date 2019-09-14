"""
This script shows how to use astropy
to extract and use FITS datafiles into
python and plot them using matplotlib.

"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# Let's load our first FITS file and print its information
# HDU means Header/Data Unit. It consists of a header and the coresponding data
HDUlist = fits.open("fits_images_all/image0.fits")
HDUlist.info()

# We can extract the data part of HDU by using index operator.
raw_image = HDUlist[0].data

# Now, we will plot our image using a colormap.
plt.imshow(raw_image)
plt.colorbar()
plt.show()