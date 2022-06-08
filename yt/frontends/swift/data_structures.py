import os

from yt.frontends.gadget.data_structures import GadgetHDF5Dataset
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import SwiftFieldInfo


class SwiftDataset(GadgetHDF5Dataset):
    _field_info_class = SwiftFieldInfo

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        """
        Checks to see if the file is a valid output from SWIFT.
        This requires the file to have the Code attribute set in the
        Header dataset to "SWIFT".
        """
        valid = True
        # Attempt to open the file, if it's not a hdf5 then this will fail:
        try:
            handle = h5py.File(filename, mode="r")
            valid = handle["Header"].attrs["Code"].decode("utf-8") == "SWIFT"
            handle.close()
        except (OSError, KeyError, ImportError):
            valid = False

        return valid

    def _set_code_unit_attributes(self):
        super()._set_code_unit_attributes()
