import os

from yt.frontends.gadget.data_structures import GadgetHDF5Dataset
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import PkdgravFieldInfo


class PkdgravDataset(GadgetHDF5Dataset):
    _field_info_class = PkdgravFieldInfo

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

    def _set_code_unit_attributes(self):
        super()._set_code_unit_attributes()
