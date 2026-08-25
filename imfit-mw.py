import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
from concurrent.futures import ThreadPoolExecutor

import sys
import importlib.util

import numpy as np
from astropy.io import fits
import pyimfit
from lmfit import Parameters, Minimizer, report_fit

# ----------------------------------------------------------------------
# Initial setup: load config, data and transform it for fitting.
# ----------------------------------------------------------------------

def load_config(config_path):
    """Load a Python configuration file as a module."""

    spec = importlib.util.spec_from_file_location("imfit_mw_config", config_path)

    if spec is None or spec.loader is None:
        raise ImportError(
            f"Cannot load configuration file: {config_path}"
        )

    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)

    return config

def expand_per_band(value, bands):
    """
    Convert either a single path or a list of paths into a list
    with one entry per band.

    A single string means the same file is used for every band.
    """

    if isinstance(value, str) or (value is None):
        return [value] * len(bands)

    if len(value) != len(bands):
        raise ValueError(
            f"Paths in config must be either one path or a list of {len(bands)} paths."
        )

    return list(value)

def load_data(config):
    """Load all observational data, the corresponding Imfit models and output paths into a single dict."""

    bands = list(config.bands)
    wavelengths = np.asarray(config.wavelengths, dtype=float)

    if len(wavelengths) != len(bands):
        raise ValueError(
            "Number of wavelengths must equal number of bands."
        )

    imfit_paths = expand_per_band(config.imfit_paths, bands)
    image_paths = expand_per_band(config.image_paths, bands)
    noise_paths = expand_per_band(config.noise_paths, bands)
    noise_types = expand_per_band(config.noise_types, bands)
    psf_paths = expand_per_band(config.psf_paths, bands)
    mask_paths = expand_per_band(config.mask_paths, bands)

    fit_paths = expand_per_band(config.fit_paths, bands)
    model_image_paths = expand_per_band(config.model_image_paths, bands)
    residual_paths = expand_per_band(config.residual_paths, bands)

    model_descs = []
    images = []
    sigmas = []
    psfs = []
    masks = []

    for band, imfit_path, image_path, noise_path, noise_type, psf_path, mask_path, fit_path, model_image_path, residual_path in zip(
        bands, imfit_paths, image_paths, noise_paths, noise_types, psf_paths, mask_paths, fit_paths, model_image_paths, residual_paths):

        # PyImfit parses the ordinary Imfit configuration file.
        model_desc = pyimfit.parse_config_file(imfit_path)

        image = fits.getdata(image_path)
        noise = fits.getdata(noise_path)

        if noise_type == "sigma":
            sigma = noise
        elif noise_type == "variance":
            sigma = np.sqrt(noise)
        elif noise_type == "weight":
            sigma = 1 / np.sqrt(noise)
        else:
            raise ValueError(
                "noise_types should contain either 'sigma', 'variance' or 'weight'."
            )

        if psf_path is not None:
            psf = fits.getdata(psf_path)
        else:
            psf = None

        if mask_path is not None:
            mask = fits.getdata(mask_path)
        else:
            mask = None

        # Create directories at the beginning, so if there is an error, it is thrown before fitting
        if fit_path is not None:
            os.makedirs(os.path.dirname(fit_path), exist_ok=True)

        if model_image_path is not None:
            os.makedirs(os.path.dirname(model_image_path), exist_ok=True)

        if residual_path is not None:
            os.makedirs(os.path.dirname(residual_path), exist_ok=True)

        model_descs.append(model_desc)
        images.append(image)
        sigmas.append(sigma)
        psfs.append(psf)
        masks.append(mask)

    return {
        "bands": bands,
        "wavelengths": wavelengths,
        "model_descs": model_descs,
        "images": images,
        "sigmas": sigmas,
        "psfs": psfs,
        "masks": masks,
        "fit_paths": fit_paths,
        "model_image_paths": model_image_paths,
        "residual_paths": residual_paths
    }

def get_model_parameter_info(model):
    """
    Extract all parameters from a PyImfit ModelDescription to use in lmfit.

    Returns
    -------
    parameters : dict

        Keys are names such as:

            X0__0
            Y0__0
            Sersic__0__n
            Sersic__0__I_e
            Exponential__1__h

        Values are dictionaries containing:

            value
            limits
            fixed
    """

    parameters = {}

    function_index = 0

    for fset in model._functionSets:
        parameters[f"X0__{function_index}"] = {"value": fset.x0.value, "limits": fset.x0.limits, "fixed": fset.x0.fixed}
        parameters[f"Y0__{function_index}"] = {"value": fset.y0.value, "limits": fset.y0.limits, "fixed": fset.y0.fixed}

        for function in fset.functionList():
            function_id = f"{function._funcName}__{function_index}"
            
            for parameter in function.parameterList():
                name = f"{function_id}__{parameter.name}"
                parameters[name] = {"value": parameter.value, "limits": parameter.limits, "fixed": parameter.fixed}

            function_index += 1

    return parameters

def validate_model_structure(data):
    """Check that all bands contain the same set of model parameters."""

    reference = get_model_parameter_info(data["model_descs"][0])
    reference_names = list(reference.keys())

    for band, model_desc in zip(data["bands"], data["model_descs"]):

        current = get_model_parameter_info(model_desc)

        if list(current.keys()) != reference_names:
            raise ValueError(
                f"Model structure in band {band} differs from the first band.\nReference: {reference_names}\nCurrent: {list(current.keys())}"
            )

    return reference_names

def select_control_bands(n_bands, degree):
    """
    If degree is provided by coupling file, control bands are auto-assigned.
    This function selects approximately uniformly spaced control-band indices.

    degree = 0 -> first band
    degree = 1 -> first and last bands
    degree = 2 -> first, middle, last
    ...
    degree = n_bands - 1 -> all bands
    """
    if degree < 0 or degree >= n_bands:
        raise ValueError(
            f"Polynomial degree {degree} is invalid for {n_bands} bands."
        )

    n_controls = degree + 1

    if n_controls == 1:
        return [0]

    indices = np.linspace(0, n_bands - 1, n_controls)
    return np.rint(indices).astype(int).tolist()

def parse_coupling_config(filename):
    """
    Parse the IMFIT-like coupling configuration.

    Each parameter line has the form

        PARAMETER  specification

    where specification can be:

        r                  # one control band
        g,z                # explicit control bands
        1                  # polynomial degree
        expr:Sersic__0__PA # lmfit expression
        function:step      # custom wavelength dependence

    X0/Y0 lines preceding FUNCTION are associated with that
    function set in the same way as in an ordinary IMFIT file.
    """

    coupling = {}

    current_function = None
    function_index = -1

    with open(filename, "r") as f:

        for line_number, line in enumerate(f, start=1):

            # Remove comments and whitespace.
            line = line.split("#", 1)[0].strip()

            if not line:
                continue

            fields = line.split()

            controls = None
            shared = None
            function_name = None

            # ----------------------------------------------------------
            # FUNCTION
            # ----------------------------------------------------------

            if fields[0] == "FUNCTION":
                if len(fields) != 2:
                    raise ValueError(
                        f"{filename}:{line_number}: FUNCTION must be followed by a function name"
                    )

                current_function = fields[1]
                function_index += 1
                continue

            # ----------------------------------------------------------
            # X0 / Y0
            # ----------------------------------------------------------

            if fields[0] in ("X0", "Y0"):
                if len(fields) < 2:
                    raise ValueError(
                        f"{filename}:{line_number}: invalid {fields[0]} specification"
                    )

                key = f"{fields[0]}__{function_index + 1}"

            else:
            # ----------------------------------------------------------
            # Ordinary model parameter
            # ----------------------------------------------------------
                if current_function is None:
                    raise ValueError(
                        f"{filename}:{line_number}: parameter encountered before FUNCTION"
                    )

                parameter_name = fields[0]
                key = f"{current_function}__{function_index}__{parameter_name}"

            for specification in fields[1:]:
                if specification.startswith("function:"):
                    function_name = specification[len("function:"):] #to remove beginning of the string
                elif specification.startswith("shared:"):
                    shared = specification[len("shared:"):]
                else:
                    controls = specification

                if function_name is None:
                    function_name = "polynomial"

            coupling[key] = (controls, shared, function_name)

    return coupling

def build_optimizer_parameters(data, config):
    """
    Construct the lmfit.Parameters object.

    The optimizer parameters are the values at the control bands,
    NOT polynomial coefficients.

    Returns
    -------
    params : lmfit.Parameters
        Parameters that will eventually be passed to lmfit.

    parameter_info : dict
        Metadata describing how each optimizer parameter maps onto
        the corresponding wavelength-dependent physical parameter.
    """

    coupling = parse_coupling_config(config.coupling_path)

    bands = data["bands"]
    wavelengths = data["wavelengths"]
    model_descs = data["model_descs"]

    reference_info = get_model_parameter_info(model_descs[0])

    params = Parameters()
    parameter_info = {}

    wavelength_functions = {"polynomial": polynomial}
    try: #if there are user-defined functions
        wavelength_functions.update(config.wavelength_functions)
    except AttributeError:
        pass

    for parameter_name in reference_info:

        controls, shared, function_name = coupling[parameter_name]

        if shared is not None:
            if shared not in coupling.keys():
                raise KeyError(
                    f"Parameter '{parameter_name}': key '{shared}' is not found."
                )
            controls, _, function_name = coupling[shared]

        if controls.isdigit(): #controls is either degree or comma-separated control bands
            degree = int(controls)
            control_indices = select_control_bands(len(bands), degree)
            control_bands = [bands[i] for i in control_indices]
        else:
            control_bands = controls.split(",")
            control_indices = [bands.index(band) for band in control_bands]

        if function_name not in wavelength_functions:
            raise ValueError(
                f"Unknown wavelength-dependence function '{function_name}' for parameter '{parameter_name}'. Available functions: {list(wavelength_functions)}"
            )

        # Check that all control bands actually exist.
        for band in control_bands:
            if band not in bands:
                raise ValueError(
                    f"Control band '{band}' for parameter '{parameter_name}' is not in bands={bands}."
                )

        # Wavelengths corresponding to control bands.
        control_wavelengths = wavelengths[control_indices]

        # Wavelengths must be distinct.
        if len(np.unique(control_wavelengths)) != len(control_wavelengths):
            raise ValueError(
                f"Control bands for '{parameter_name}' have duplicate wavelengths."
            )

        # Information for the corresponding parameter in every band.
        band_info = []

        for model_desc in model_descs:
            info = get_model_parameter_info(model_desc)[parameter_name]
            band_info.append(info)

        initial_values = []
        limits = []
        fixed_flags = []

        for band_index in control_indices:
            info = band_info[band_index]
            initial_values.append(info["value"])
            limits.append(info["limits"])
            fixed_flags.append(info["fixed"])

        # Create one lmfit parameter for each control-band value.
        optimizer_names = []

        for band, value, limit, fixed in zip(control_bands, initial_values, limits, fixed_flags):

            optimizer_name = (f"{parameter_name}__{band}")
            optimizer_names.append(optimizer_name)
            kwargs = {"value": value, "vary": not fixed}

            if limit is not None and not fixed:
                kwargs["min"] = limit[0]
                kwargs["max"] = limit[1]

            if shared is not None:
                kwargs["expr"] = f"{shared}__{band}"

            params.add(optimizer_name, **kwargs)

        parameter_info[parameter_name] = {
            "control_bands": control_bands,
            "control_indices": control_indices,
            "control_wavelengths": control_wavelengths,
            "optimizer_names": optimizer_names,
            "function_name": function_name,
            "function": wavelength_functions[function_name],
        }

    return params, parameter_info

# ----------------------------------------------------------------------
# For each function evaluation, fitting parameters are transformed
# into photometric parameters in each band.
# ----------------------------------------------------------------------

def polynomial(wavelengths, control_wavelengths, control_values):
    """
    Polynomial interpolation through the control-band values.

    Custom functions are possible, but should have same parameters and returns as this one.

    Parameters
    ----------
    wavelengths : array
        Wavelengths where the polynomial should be evaluated.

    control_wavelengths : array
        Wavelengths of the fitted control parameters.

    control_values : array
        Fitted parameter values at the control wavelengths.

    Returns
    -------
    values : ndarray
        Parameter values at all requested wavelengths.
    """
    degree = len(control_values) - 1
    coefficients = np.polynomial.polynomial.polyfit(control_wavelengths, control_values, degree)

    return np.polynomial.polynomial.polyval(wavelengths, coefficients)

def unpack_parameters(params, data, parameter_info):
    """
    Convert the fitted control-band values into actual model parameters
    for every band.

    This is the operation that will later be used inside residual().

    Parameters
    ----------
    params : lmfit.Parameters
        Current optimizer parameters.

    Returns
    -------
    band_parameters : list of dict
        One dictionary per band containing actual Imfit parameter values.
    """

    bands = data["bands"]
    wavelengths = data["wavelengths"]

    band_parameters = [{} for _ in bands]

    for parameter_name, info in parameter_info.items():

        control_values = np.array([params[name].value for name in info["optimizer_names"]], dtype=float)

        # info["function"] is callable defining wavelength dependence, default is polynomial().
        values = info["function"](wavelengths, info["control_wavelengths"], control_values)

        for band_index, value in enumerate(values):
            band_parameters[band_index][parameter_name] = value

    return band_parameters


# ----------------------------------------------------------------------
# Perform fitting with lmfit
# ----------------------------------------------------------------------

class MultiBandModel:
    """
    Persistent PyImfit models for all bands.

    One PyImfit.Imfit object is created for each band and reused for
    every residual evaluation.
    """

    def __init__(self, data):

        self.bands = data["bands"]
        self.model_descs = data["model_descs"]
        self.images = data["images"]
        self.sigmas = data["sigmas"]
        self.psfs = data["psfs"]
        self.masks = data["masks"]

        self.fitters = []

        # Mapping:
        #
        #   band_parameters name -> index in PyImfit parameter vector
        #
        self.parameter_indices = []

        for model_desc, image, sigma, psf, mask in zip(self.model_descs, self.images, self.sigmas, self.psfs, self.masks):

            fitter = pyimfit.Imfit(model_desc, psf=psf)

            if mask is None:
                fitter.loadData(image, error = sigma)
            else:
                fitter.loadData(image, error = sigma, mask = mask)

            self.fitters.append(fitter)

            # Construct names in exactly the same way as
            # get_model_parameter_info().
            names = []

            function_index = 0

            for fset in model_desc._functionSets:

                names.append(f"X0__{function_index}")
                names.append(f"Y0__{function_index}")

                for function in fset.functionList():

                    function_id = f"{function._funcName}__{function_index}"

                    for parameter in function.parameterList():
                        names.append(
                            f"{function_id}__{parameter.name}"
                        )

                    function_index += 1

            # PyImfit's parameter vector follows the same ordering
            # as numberedParameterNames.
            pyimfit_names = fitter.numberedParameterNames

            if len(names) != len(pyimfit_names):
                raise RuntimeError(
                    f"Number of parameters differs between our model description and PyImfit for band."
                )

            self.parameter_indices.append({name: i for i, name in enumerate(names)})

        self.pool = ThreadPoolExecutor(max_workers=len(self.bands))

    def make_model_image(self, band_index, band_parameters):
        """
        Generate a model image for one band.

        Parameters
        ----------
        band_index : int
            Index of the band.

        band_parameters : dict
            Physical model parameters for this band.

        Returns
        -------
        model_image : ndarray
            Model image.
        """

        fitter = self.fitters[band_index]
        indices = self.parameter_indices[band_index]

        n_parameters = len(indices)

        parameter_vector = np.empty(n_parameters, dtype=float)

        for name, value in band_parameters.items():
            parameter_vector[indices[name]] = value

        fitter.computeFitStatistic(parameter_vector)

        return fitter.getModelImage()

    def residual_one_band(self, band_index, band_parameters):
        """
        Compute the residual vector for one band.

        Parameters
        ----------
        band_index : int
            Band index.

        band_parameters : dict
            Physical model parameters for this band.

        Returns
        -------
        residual : ndarray
            Flattened residual vector for this band.
        """

        model_image = self.make_model_image(band_index, band_parameters)

        image = self.images[band_index]
        sigma = self.sigmas[band_index]
        mask = self.masks[band_index]

        # Guard against numerical problems in the model.
        if not np.all(np.isfinite(model_image)):
            n = image.size if mask is None else np.count_nonzero(~mask.astype(bool))
            return np.full(n, 1.0e100)

        res = (image - model_image) / sigma

        if mask is not None:
            res = res[~mask.astype(bool)]

        # Guard against numerical problems in the residual itself.
        if not np.all(np.isfinite(res)):
            return np.full(res.size, 1.0e100)

        return res.ravel()

    def residual(self, params, data, parameter_info):

        band_parameters = unpack_parameters(params, data, parameter_info)

        futures = [
            self.pool.submit(self.residual_one_band, band_index, parameters)
            for band_index, parameters in enumerate(band_parameters)
            ]

        residuals = [future.result() for future in futures]

        return np.concatenate(residuals)

    def make_string_current_params(self, band_index, band_parameters):
        fitter = self.fitters[band_index]
        indices = self.parameter_indices[band_index]
        n_parameters = len(indices)
        parameter_vector = np.empty(n_parameters, dtype=float)
        for name, value in band_parameters.items():
            parameter_vector[indices[name]] = value

        model_desc = fitter.getModelDescription()
        fit_desc_init = model_desc.getStringDescription(noLimits = True)
        fit_desc = []

        i = 0
        for line in fit_desc_init:
            try:
                rep_cand = line.split()[1]
                rep_float = float(rep_cand)
                line_mod = line.replace(rep_cand, f"{parameter_vector[i]:.6f}")
                fit_desc.append(line_mod)
                i += 1
            except IndexError:
                fit_desc.append(line)
            except ValueError:
                fit_desc.append(line)

        return fit_desc

    def close(self):
        self.pool.shutdown(wait=True)


def run_fit(data, params, parameter_info, minimizer_kwargs):
    """
    Run the lmfit optimization.
    """
    mb_model = MultiBandModel(data)
    minimizer = Minimizer(mb_model.residual, params, fcn_args=(data, parameter_info))
    result = minimizer.minimize(**minimizer_kwargs)
    return result, mb_model

# ----------------------------------------------------------------------
# Save fit results
# ----------------------------------------------------------------------

def save_results(data, result, parameter_info, mb_model):
    """
    Save fitted per-band .imfit files, model images, and residual images.
    """

    # Convert optimizer vector into parameters for each individual band.
    band_parameters = unpack_parameters(result.params, data, parameter_info)

    for b, band in enumerate(data["bands"]):

        image = data["images"][b]
        model_desc = data["model_descs"][b]
        fit_desc = mb_model.make_string_current_params(b, band_parameters[b])
        model_image = mb_model.make_model_image(b, band_parameters[b])
        residual_image = image - model_image

        fit_path = data["fit_paths"][b]
        model_image_path = data["model_image_paths"][b]
        residual_path = data["residual_paths"][b]

        if fit_path is not None:
            with open(fit_path, "w") as file:
                for line in fit_desc:
                    file.write(line)

        if model_image_path is not None:
            fits.writeto(model_image_path, model_image, overwrite=True)

        if residual_path is not None:
            fits.writeto(residual_path, residual_image, overwrite=True)

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main(config_path):
    config = load_config(config_path)
    data = load_data(config)
    validate_model_structure(data)
    params, parameter_info = build_optimizer_parameters(data, config)
    band_parameters = unpack_parameters(params, data, parameter_info)

    print("\nStarting multi-band fit...")

    result, mb_model = run_fit(data, params, parameter_info, config.minimizer_kwargs)
    report_fit(result)
    save_results(data, result, parameter_info, mb_model)
    mb_model.close()
    return (data, result, parameter_info, mb_model)

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python imfit_mw.py config.py")
        sys.exit(1)

    main(sys.argv[1])
