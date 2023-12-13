import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

# Get a list of file names in the imagesVal folder
image_dir = 'training_data_noise00004/preprocessedDataset/imagesVal/'
image_filenames = [f for f in os.listdir(image_dir) if f.endswith('.mat')]

# Get a list of file names in the labelsVal folder
label_dir = 'training_data_noise00004/preprocessedDataset/labelsVal/'
label_filenames = [f for f in os.listdir(label_dir) if f.endswith('.mat')]

def makeplots(i):

    # Load the first image and label data from the .mat files
    image_data = scipy.io.loadmat(os.path.join(image_dir, image_filenames[i]))
    label_data = scipy.io.loadmat(os.path.join(label_dir, label_filenames[i]))

    max_val = np.max(label_data['labelone'])
    max_ind = np.argmax(label_data['labelone'])
    max_x, max_y, max_z = np.unravel_index(max_ind, label_data['labelone'].shape)

    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot the image data on the left subplot
    dataslice = image_data['dataone'][max_x, :, :]
    im1 = ax1.imshow(dataslice)
    ax1.set_aspect('equal')
    ax1.set_title('Image')
    cbar1 = fig.colorbar(im1, ax=ax1)

    # Plot the label data on the right subplot
    labelslice = label_data['labelone'][max_x, :, :]
    im2 = ax2.imshow(labelslice)
    ax2.set_aspect('equal')
    ax2.set_title('Label')
    cbar2 = fig.colorbar(im2, ax=ax2)

    # Save the figure with the same name as the input file in the Val folder
    filename, _ = os.path.splitext(image_filenames[i])
    fig.savefig(os.path.join('training_data_noise00004/Val', f'{filename}_{max_x}.png'))

    plt.close(fig)

    if i%100 == 0:
        print(i)

def main():
    with multiprocessing.Pool(25) as pool:
        pool_outputs = pool.map(makeplots, range(0, 19000))

if __name__ == '__main__':
    main()
