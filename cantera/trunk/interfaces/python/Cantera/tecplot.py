"""TECPLOT utilities."""

def write_TECPLOT_zone(fname, title, zone, names, npts, nvar, append, data):
    """
    Write a TECPLOT zone specification to generate line plots of multiple
    variables.
    fname    -- file name
    title    -- plot title
    zone     -- zone name
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
    f.write('TITLE     = "' + title + '"\n')
    f.write('VARIABLES = \n')
    for nm in names:
        f.write('"'+nm+'"\n')
    f.write('ZONE T="'+zone+'"\n')
    f.write(' I='+`npts`+',J=1,K=1,F=POINT\n')
    f.write('DT=(')
    for n in range(nvar):
        f.write('SINGLE ')
    f.write(')\n')
    for j in range(npts):
        for n in range(nvar):
            f.write('%10.4e ' % data[j,n])
        f.write('\n')
    f.close()
