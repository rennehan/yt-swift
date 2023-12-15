import numpy as np
from numpy.testing import assert_almost_equal

from yt import load
from yt.frontends.swift.api import SwiftDataset
from yt.testing import ParticleSelectionComparison, requires_file, requires_module
from yt.utilities.on_demand_imports import _h5py as h5py

# This maps from field names to weight field names to use for projections
fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), ("gas", "density")),
        (("gas", "metallicity"), ("gas", "density")),
        (("gas", "O_metallicity"), ("gas", "density")),
        (("gas", "velocity_magnitude"), None),
    ]
)

# Combined the tests for loading a file and ensuring the units have been
# implemented correctly to save time on re-loading a dataset
@requires_module("h5py")
@requires_file(keplerian_ring)
def test_non_cosmo_dataset():
    ds = load(keplerian_ring)
    assert type(ds) is SwiftDataset


@requires_ds(g64, big_data=True)
def test_gizmo_64():
    ds = yt.load(g64)
    assert isinstance(ds, GizmoDataset)
    for test in sph_answer(ds, "snap_N64L16_135", 524288, fields):
        test_gizmo_64.__name__ = test.description
        yield test


@requires_file(gmhd)
def test_gizmo_mhd():
    """
    Magnetic fields should be loaded correctly when they are present.
    """
    ds = yt.load(gmhd, bounding_box=gmhd_bbox, unit_system="code")
    ad = ds.all_data()
    ptype = "PartType0"

    # Test vector magnetic field
    fmag = "particle_magnetic_field"
    f = ad[ptype, fmag]
    assert str(f.units) == "code_magnetic"
    assert f.shape == (409013, 3)

    # Test component magnetic fields
    for axis in "xyz":
        f = ad[ptype, fmag + "_" + axis]
        assert str(f.units) == "code_magnetic"
        assert f.shape == (409013,)

@requires_module("h5py")
@requires_file(keplerian_ring)
def test_non_cosmo_dataset_selection():
    ds = load(keplerian_ring)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()


@requires_module("h5py")
@requires_file(EAGLE_6)
def test_cosmo_dataset():
    ds = load(EAGLE_6)
    assert type(ds) == SwiftDataset

    field = ("gas", "density")
    ad = ds.all_data()
    yt_density = ad[field]
    yt_coords = ad[(field[0], "position")]

    # load some data the old fashioned way
    fh = h5py.File(ds.parameter_filename, mode="r")
    part_data = fh["PartType0"]

    # set up a conversion factor by loading the unit mas and unit length in cm,
    # and then converting to proper coordinates
    units = fh["Units"]
    units = dict(units.attrs)
    density_factor = float(units["Unit mass in cgs (U_M)"])
    density_factor /= float(units["Unit length in cgs (U_L)"]) ** 3

    # add the redshift factor
    header = fh["Header"]
    header = dict(header.attrs)
    density_factor *= (1.0 + float(header["Redshift"])) ** 3

    # now load the raw density and coordinates
    raw_density = part_data["Density"][:].astype("float64") * density_factor
    raw_coords = part_data["Coordinates"][:].astype("float64")
    fh.close()

    # sort by the positions - yt often loads in a different order
    ind_raw = np.lexsort((raw_coords[:, 2], raw_coords[:, 1], raw_coords[:, 0]))
    ind_yt = np.lexsort((yt_coords[:, 2], yt_coords[:, 1], yt_coords[:, 0]))
    raw_density = raw_density[ind_raw]
    yt_density = yt_density[ind_yt]

    # make sure we are comparing fair units
    assert str(yt_density.units) == "g/cm**3"

    # make sure the actual values are the same
    assert_almost_equal(yt_density.d, raw_density)


@requires_module("h5py")
@requires_file(EAGLE_6)
def test_cosmo_dataset_selection():
    ds = load(EAGLE_6)
    psc = ParticleSelectionComparison(ds)
    psc.run_defaults()
