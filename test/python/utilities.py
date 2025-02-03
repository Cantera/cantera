import numpy as np
from pathlib import Path, PurePath
from pytest import approx
import warnings
from ruamel import yaml


def load_yaml(yml_file):
    """
    Load YAML data from file using the "safe" loading option.
    """
    yaml_parser = yaml.YAML(typ="safe")
    with open(yml_file, "rt", encoding="utf-8") as stream:
        return yaml_parser.load(stream)

def compare(data, reference_file, rtol=1e-8, atol=1e-12):
    """
    Compare an array with a reference data file, or generate the reference
    file if it does not exist.
    """
    data = np.array(data)
    if Path(reference_file).is_file():
        # Compare with existing output file
        ref = np.genfromtxt(reference_file)
        assert data.shape == ref.shape
        for i in range(ref.shape[0]):
            assert ref[i] == approx(data[i], rel=rtol, abs=atol)
    else:
        # Generate the output file for the first time
        warnings.warn('Generating test data file:' + Path(reference_file).resolve())
        np.savetxt(reference_file, data, fmt='%.10e')

def compareProfiles(reference, sample, rtol=1e-5, atol=1e-12, xtol=1e-5):
    """
    Compare two 2D arrays of spatial or time profiles. Each data set should
    contain the time or space coordinate in the first column and data series
    to be compared in successive columns.

    The coordinates in each data set do not need to be the same: The data from
    the second data set will be interpolated onto the coordinates in the first
    data set before being compared. This means that the range of the "sample"
    data set should be at least as long as the "reference" data set.

    After interpolation, each data point must satisfy a combined relative and absolute
    error criterion specified by `rtol` and `atol`.

    If the comparison succeeds, this function returns `None`. If the comparison
    fails, a formatted report of the differing elements is returned.
    """
    if isinstance(reference, (str, PurePath)):
        reference = np.genfromtxt(reference, delimiter=',').T
    else:
        reference = np.asarray(reference).T

    if isinstance(sample, (str, PurePath)):
        sample = np.genfromtxt(sample, delimiter=',').T
    else:
        sample = np.asarray(sample).T

    assert reference.shape[0] == sample.shape[0]

    nVars = reference.shape[0]
    nTimes = reference.shape[1]

    bad = []
    template = '{0:9.4e}  {1: 3d}   {2:14.7e}  {3:14.7e}  {4:9.3e}  {5:9.3e}  {6:9.3e}'
    for i in range(1, nVars):
        scale = max(max(abs(reference[i])), np.ptp(reference[i]),
                    max(abs(sample[i])), np.ptp(sample[i]))
        slope = np.zeros(nTimes)
        slope[1:] = np.diff(reference[i]) / np.diff(reference[0]) * np.ptp(reference[0])

        comp = np.interp(reference[0], sample[0], sample[i])
        for j in range(nTimes):
            a = reference[i,j]
            b = comp[j]
            abserr = abs(a-b)
            relerr = abs(a-b) / (scale + atol)

            # error that can be accounted for by shifting the profile along
            # the time / spatial coordinate
            xerr = abserr / (abs(slope[j]) + atol)

            if abserr > atol and relerr > rtol and xerr > xtol:
                bad.append((reference[0][j], i, a, b, abserr, relerr, xerr))

    footer = []
    maxrows = 10
    if len(bad) > maxrows:
        bad.sort(key=lambda row: -row[5])
        footer += ['Plus {0} more points exceeding error thresholds.'.format(len(bad)-maxrows)]
        bad = bad[:maxrows]

    if bad:
        header = ['Failed series comparisons:',
                  'coordinate  comp.  reference val    test value    abs. err   rel. err    pos. err',
                  '----------  ---   --------------  --------------  ---------  ---------  ---------']
        return '\n'.join(header + [template.format(*row) for row in bad] + footer)
    else:
        return None
