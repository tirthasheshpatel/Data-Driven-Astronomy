"""
This python script shows how to perform mean stacking
on a list of images and plot our final results to
make a detection more clear.

Remember, our images have been croped and the probable
pulsar detections have been centered.

"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

def mean_stack(image_list):
    """
    This function calculates the mean of each
    pixel accross all the images. The noise cancles
    out and the detection becomes more clear.
    This is used for crearer view of a detection.
    """

    # First, we make a numpy array of all the images.
    data = np.array([fits.open(image_file)[0].data for image_file in image_list]) # (nb_images, pixels, pixels)
    print(f"The shape of list of images: {data.shape}")

    # Now, we calculate the mean across the number of images and not accross number of pixels.
    mean_img = np.mean(data, axis = 0) 
    print(f"The shape of final stacked image: {mean_img.shape}")

    return mean_img


# HDU means Header/Data Unit. It consists of a header and the coresponding data
HDUlist = fits.open("big_dataset_for_stacking/0001.fits")

# We can extract the data part of HDU by using index operator.
# The data part consists of our image.
raw_image = HDUlist[0].data

# Now, we will plot our image using a colormap.
plt.imshow(raw_image)
plt.colorbar()
plt.show()

# Let's stack a 1000 images to get a clearer view of the detection.
images_to_stack = [f"big_dataset_for_stacking/{i:04d}.fits" for i in range(1000)]
stacked_image = mean_stack(images_to_stack)

# Now, we will plot our image using a colormap.
plt.imshow(stacked_image)
plt.colorbar()
plt.show()