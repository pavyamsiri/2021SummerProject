"""
The main function.
Reads in a csv file that stores star data and plots their light curve and periodogram.
"""
# Standard libraries
import sys
import os
from time import time
from datetime import datetime
import logging
import argparse
import re

# External libraries
from matplotlib import pyplot as plt
import pandas as pd

# Interal libraries
from utils import grab_lightcurve, create_plots
from star import StarData


# Default save folder
DEFAULT_FOLDER = "../Data/Plots/Tim_White_Targets"
# Default lightkurve cache
DEFAULT_CACHE = "../cache"
# Default file format
DEFAULT_EXTENSION = "png"

# Argument parser
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", nargs=1,
                    help="Source data file.")
parser.add_argument("-s", "--save", nargs="?",
                    const=DEFAULT_FOLDER, default=DEFAULT_FOLDER,
                    help="Plot save folder.")
parser.add_argument("-c", "--cache", nargs="?",
                    const=DEFAULT_CACHE, default=DEFAULT_CACHE,
                    help="Cache folder.")
parser.add_argument("-e", "--extension", nargs="?",
                    const=DEFAULT_EXTENSION, default=DEFAULT_EXTENSION,
                    help="Save file format.")


def setup_logger():
    """Setups the logger."""

    # Log filename
    now = datetime.now()
    log_filename = "Logs/" + now.strftime("%d_%m_%Y__%H_%M_%S.txt")

    # Log to both terminal and file
    log = logging.getLogger(__name__)
    stream_handler = logging.StreamHandler(sys.stdout)
    log.addHandler(stream_handler)
    file_handler = logging.FileHandler(log_filename)
    log.addHandler(file_handler)

    # Debug level
    log.setLevel(logging.DEBUG)

    # Log format
    formatter = logging.Formatter(
        "{asctime} {levelname:<8} {message}",
        style="{"
    )
    stream_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # Also redirect print messages and error messages to file output
    file_output = open(log_filename, "a")

    return log, file_output


def main(source_file: str, save_dir: str, cache_dir: str, file_extension: str, log: logging.Logger) -> None:
    """Reads in a csv file that stores star data and plots them.

    Parameters
    ----------
    source_file : str
        a valid file path to the plotting data
    save_dir : str
        a directory path to output saved plots to
    cache_dir : str
        a directory path to act as the cache for light curve data
    file_extension : str
        the file type of the saved plots i.e. .png or .svg
    log : logging.Logger
        the logger
    """

    # Create save folder if it doesn't exist
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    # Create cache folder if it doesn't exist
    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)

    # Read in csv data
    data_table = pd.read_csv(source_file)

    # Illegal file/folder name characters
    illegal_characters = r'[\\~#%&*{}\/:<>?|"-]'

    start_time = time()

    for index, star in enumerate(data_table["Star"]):
        # Store star data in named tuple
        current = StarData(
            star,
            data_table["Main ID"][index],
            data_table["TIC"][index],
            data_table["Sectors"][index],
            data_table["SpT"][index],
            data_table["Diam"][index],
            data_table["Err"][index],
            data_table["Instrument"][index]
        )

        plot_name = re.sub(r"\s+", "_", current.main_id)
        plot_name = re.sub(r"V\*", "Vstar", plot_name)
        spectral_prefix = re.sub(illegal_characters, '_', current.spectral_type)
        plot_filename = f"{save_dir}/{spectral_prefix}-{plot_name}"
        if data_table["Fast"][index]:
            plot_filename += "-fast"
        plot_filename += f".{file_extension}"

        log.info("\nIndex %d\n%s", index, current)
        # Skip if plot already exists
        if os.path.exists(plot_filename):
            log.info("%s's plot already exists. Continuing to next star...\n", current.star)
            continue

        # Download lightcurve
        light_curve = grab_lightcurve(
            f"TIC{current.tic}",
            cache_dir,
            data_table["Fast"][index]
            )

        # If no valid light curves exist skip but note down
        if light_curve is None:
            log.error("Invalid light curve for %s. Continuing to next star...\n", current.star)
            continue

        numax = {"ATL": None, "TIC": None}
        if not pd.isnull(data_table.loc[index, "ATL numax (muHz)"]):
            numax["ATL"] = float(data_table.loc[index, "ATL numax (muHz)"])
        if not pd.isnull(data_table.loc[index, "TIC 0.8M numax (muHz)"]):
            numax["TIC"] = []
            numax["TIC"].append(float(data_table.loc[index, "TIC 0.8M numax (muHz)"]))
            numax["TIC"].append(float(data_table.loc[index, "TIC 1.2M numax (muHz)"]))
            numax["TIC"].append(float(data_table.loc[index, "TIC 2.4M numax (muHz)"]))

        # Create plots and save as a file
        fig = create_plots(current, light_curve, numax)
        fig.savefig(plot_filename)
        plt.close(fig)

        log.info("Saving Star %s...\n", index)

    # Close all plot output streams
    plt.close("all")

    # Record how long it took to run
    elapsed = time() - start_time
    log.info("This took %d minutes and %.2f seconds to complete!", elapsed//60, elapsed % 60)

if __name__ == "__main__":
    # Parse arguments
    args = parser.parse_args()

    # Set up logger
    logger, log_file = setup_logger()

    # Check path is valid
    if args.data is None or not os.path.isfile(str(args.data[0])):
        logger.error(f"The given path `{args.data[0]}` is invalid!")
        sys.exit()

    # Process arguments
    DATA_SOURCE = args.data[0]
    SAVE_FOLDER = args.save
    DATA_CACHE = args.cache
    FILE_EXTENSION = args.extension if os.path.exists(args.extension) else DEFAULT_EXTENSION

    # Start main process
    main(DATA_SOURCE, SAVE_FOLDER, DATA_CACHE, FILE_EXTENSION, logger)

    # Close log file
    log_file.close()
