"""EXCEL CSV file utilities."""

def write_CSV_data(fname, names, npts, nvar, append, data):
    """
    Write CSV data that can be imported into Excel

    fname    -- file name
    names    -- sequence of variable names
    npts     -- number of data points
    nvar     -- number of variables
    append   -- if > 0, append to plot file, otherwise overwrite
    data     -- object to generate plot data. This object must have a
                method 'value', defined so that data.value(j,n) returns
                the value of variable n at point j.
    """

    if append > 0:
        f = open(fname,'a')
    else:
        f = open(fname,'w')
    for nm in names:
        f.write(nm+',')
    f.write('\n')
    for j in range(npts):
        for n in range(nvar):
            f.write('%10.4e, ' % data.value(j,n))
        f.write('\n')
    f.close()
