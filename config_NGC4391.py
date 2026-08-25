from pathlib import Path

# Directory containing this configuration file.
# All relative paths below are interpreted relative to it.
BASE_DIR = Path(__file__).resolve().parent


# Bands parameters
bands = ["g", "r", "i", "z"]
wavelengths = [481, 617, 752, 866] # Units are irrelevant, but should be consistent.

# Paths and properties of data.
# Single string means common file for all bands. 

# Standard Imfit input files containing initial guesses.
imfit_paths = [f"NGC4391/{band}/input.imfit" for band in bands]

# Images
image_paths = [f"NGC4391/{band}/image.fits" for band in bands]

# Error/noise maps
noise_paths = [f"NGC4391/{band}/sigma.fits" for band in bands]

# Error/noise map type: either "sigma", "variance" or "weight"
noise_types = "sigma"

# PSFs
# Leave None for no PSF.
psf_paths = [f"NGC4391/{band}/psf.fits" for band in bands]

# Masks
# Leave None for empty mask.
#mask_paths = "NGC4391/mask.fits"
mask_paths = None

# Result, model and residual paths to save.
fit_paths = [f"NGC4391/{band}/fit.imfit" for band in bands]
model_image_paths = [f"NGC4391/{band}/model.fits" for band in bands]
residual_paths = [f"NGC4391/{band}/residual.fits" for band in bands]

# File with IMFIT-like formatting containing wavelength dependence polynomial degree and control bands
coupling_path = "coupling_NGC4391.imfit"

# kwargs for lmfit minimizer
minimizer_kwargs = {
	"method": "least_squares",
	"xtol": 1e-6,
	"ftol": 1e-6,
	"gtol": 1e-6,
	"max_nfev": 50000,
	"verbose": 2
}