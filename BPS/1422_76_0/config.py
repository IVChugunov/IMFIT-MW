from pathlib import Path

# Directory containing this configuration file.
# All relative paths below are interpreted relative to it.
BASE_DIR = Path(__file__).resolve().parent


# Bands parameters
bands = ["g", "r", "z"]
wavelengths = [4686, 6166, 8932] # Units are irrelevant, but should be consistent.

# Paths and properties of data.
# Single string means common file for all bands. 

# Standard Imfit input files containing initial guesses.
imfit_paths = [f"{BASE_DIR}/{band}/input_nfree.imfit" for band in bands]

# Images
image_paths = [f"{BASE_DIR}/{band}/image.fits" for band in bands]

# Error/noise maps
noise_paths = [f"{BASE_DIR}/{band}/invvar.fits" for band in bands]

# Error/noise map type: either "sigma", "variance" or "weight"
noise_types = "weight"

# PSFs
# Leave None for no PSF.
psf_paths = [f"{BASE_DIR}/{band}/psf.fits" for band in bands]

# Masks
# Leave None for empty mask.
mask_paths = f"{BASE_DIR}/mask_wdust.fits"

# Result, model and residual paths to save.
fit_paths = [f"{BASE_DIR}/{band}/fit.imfit" for band in bands]
model_image_paths = [f"{BASE_DIR}/{band}/model.fits" for band in bands]
residual_paths = [f"{BASE_DIR}/{band}/residual.fits" for band in bands]

# File with IMFIT-like formatting containing wavelength dependence polynomial degree and control bands
coupling_path = f"{BASE_DIR}/coupling.imfit"

# kwargs for lmfit minimizer
minimizer_kwargs = {
	"method": "least_squares",
	"xtol": 1e-6,
	"ftol": 1e-6,
	"gtol": 1e-6,
	"max_nfev": 50000,
	"verbose": 2
}