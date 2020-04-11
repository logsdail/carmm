from ase.visualize.plot import animate


def write_animation(filename, images, writer=None,
                    interval=200, save_count=100,
                    save_parameters=None, ax=None, **kwargs):
    import matplotlib.pyplot as plt

    if writer is None and filename.endswith('.gif'):
        # Alternative is pillow (PIL).
        writer = 'imagemagick'

    if save_parameters is None:
        save_parameters = {}

    if ax is None:
        ax = plt.gca()

    animation = animate(images, ax=ax,
                        interval=interval, save_count=save_count,
                        **kwargs)
    animation.save(filename, writer=writer,
                   **save_parameters)


# Shortcuts for ase.io.formats (guessing file type from extension):
write_gif = write_animation
write_mp4 = write_animation
