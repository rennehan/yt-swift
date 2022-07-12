import os

import numpy as np

from yt.frontends.gadget.data_structures import GadgetHDF5Dataset
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.funcs import only_on_root
from yt.utilities.logger import ytLogger as mylog
from .fields import PkdgravFieldInfo


class PkdgravDataset(GadgetHDF5Dataset):
    _field_info_class = PkdgravFieldInfo

    def _parse_parameter_file(self):

        hvals = self._get_hvals()

        # Rennehan
        pkdgrav = False
        try:
            units_test = self._get_info_attributes("Units")
            save_test = units_test["MsolUnit"]
            cosmology_hdr = self._get_info_attributes("Cosmology")
            pkdgrav = True
        except:
            pass

        # Rennehan
        only_on_root(mylog.info, "PkdgravDataset!")
        self.dimensionality = 3
        self.dimensionality = 3

        #self.refine_by = 2
        #self.parameters["HydroMethod"] = "sph"

        # This should be done more flexible!
        self.domain_right_edge = np.array([0.5, 0.5, 0.5])
        self.domain_left_edge = np.array([-0.5, -0.5, -0.5])

        self.domain_dimensions = np.ones(3, "int32")
        self._periodicity = (True, True, True)

        # Again, account for non-cosmological setups
        self.cosmological_simulation = 1

        try:
            self.current_redshift = float(hvals["Redshift"])
        except KeyError:
            self.current_redshift = 0.0
            only_on_root(mylog.info, "Redshift is not set in Header. Assuming z=0.")

        self.omega_lambda = float(cosmology_hdr["Omega_lambda"])
        self.omega_matter = float(cosmology_hdr["Omega_m"])
        self.hubble_constant = float(cosmology_hdr["HubbleParam"])

        self.current_time = float(hvals['Time'])
        only_on_root(mylog.info, "Current time: %e", self.current_time)
        self.parameters = hvals

        prefix = os.path.abspath(
            os.path.join(
                os.path.dirname(self.parameter_filename),
                os.path.basename(self.parameter_filename).split(".", 1)[0],
            )
        )

        if hvals["NumFiles"] > 1:
            # Probably broken
            for t in (
                f"{prefix}.%(num)s{self._suffix}",
                f"{prefix}.gad.%(num)s{self._suffix}",
            ):
                if os.path.isfile(t % {"num": 0}):
                    self.filename_template = t
                    break
            else:
                raise RuntimeError("Could not determine correct data file template.")
        else:
            self.filename_template = self.parameter_filename

        self.file_count = hvals["NumFiles"]

    def _set_code_unit_attributes(self):
        units = self._get_info_attributes("Units")

        if self.cosmological_simulation == 1:
            msg = "Assuming length units are in comoving kpc"
            only_on_root(mylog.info, msg)
            self.length_unit = self.quan(
                float(units["KpcUnit"]), "kpccm"
            )
        else:
            msg = "Assuming length units are in physical kpc"
            only_on_root(mylog.info, msg)
            self.length_unit = self.quan(float(units["KpcUnit"]), "kpc")

        self.mass_unit = self.quan(float(units["MsolUnit"]), "Msun")
        self.time_unit = self.quan(float(units["SecUnit"]), "s")
        self.velocity_unit = self.quan(float(units["KmPerSecUnit"]), "km * s**-1 * a**-1")
        self.temperature_unit = self.quan(1.0, "K")
        self.specific_energy_unit = self.quan(float(units["ErgPerGmUnit"]), "erg/g")


    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        """
        Check header for PKDGRAV version.
        """
        valid = True
        # Attempt to open the file, if it's not a hdf5 then this will fail:
        try:
            handle = h5py.File(filename, mode="r")
            pkdgrav_version = handle["Header"].attrs["PKDGRAV version"].decode("utf-8")
            handle.close()
        except (OSError, KeyError, ImportError):
            valid = False

        return valid
