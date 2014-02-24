"""
Reaction path diagrams.

To create a simple reaction path diagram:

>>> import rxnpath
>>> element = 'C'
>>> rxnpath.write(gas, element, file)

Object 'gas' must an instance of a class derived from class 'Kinetics'
(for example class IdealGasMix). The diagram layout is written to
'file'.  The output must be postprocessed with program 'dot', which is
part of the GraphViz package. To create a Postscript plot:

dot -Tps rp.dot > rp.ps

Other output formats are also supported by dot, including gif, pcl,
jpg, png, and svg

For more control over the graph properties, create a PathDiagam object
and pass it to the 'write' procedure.

PathDiagram keyword options:

  diagram type:
  -- detailed       'true' or 'false'
  -- type           'both' or 'net' (forward and reverse arrows,
                                     or net arrow)
  -- dot_options    options passed through to 'dot'

  colors:
  -- normal_color   color for normal-weight lines
  -- bold_color     color for bold-weight lines
  -- dashed_color   color for dashed lines

  thresholds:
  -- threshold         min relative strength for a path to be shown
  -- normal_threshold  min relative strength for normal-weight path
                       Below this value, paths are dashed.
  -- bold_threshold    min relative strength for bold-weight path

"""

import _cantera

class PathDiagram:
    def __init__(self, **options):
        self.__rdiag_id = _cantera.rdiag_new()
        self.setOptions({"detailed":"true",
                    "dashed_color":"gray",
                    "bold_color":"red",
                    "normal_color":"steelblue",
                    "scale":-1,
                    "dot_options":'center=1;margin=0;size="5,6";page="5,6";ratio=compress;fontname=Arial;',
                    "title":"-",
                    "arrow_width":-5,
                    "threshold":0.001,
                    "bold_threshold":0.2,
                    "normal_threshold":0.01,
                    "label_threshold":0.001,
                    "flow_type":"net"}
                    )
        self.setOptions(options)

    def __del__(self):
        _cantera.rdiag_del(self.__rdiag_id)

    def id(self):
        return self.__rdiag_id

    def write(self, fmt, file):
        _cantera.rdiag_write(self.__rdiag_id, fmt, file)

    def add(self, other):
        _cantera.rdiag_add(self.__rdiag_id, other.id())

    def findMajorPaths(self, a, threshold = 0.0):
        _cantera.rdiag_findMajor(self.__rdiag_id, threshold, a)

    def displayOnly(self, node=-1):
        _cantera.rdiag_displayOnly(self.__rdiag_id, node)

    def setOptions(self, options):
        for o in options.keys():
            v = options[o]
            if o == "detailed":
                if v == "true":
                    _cantera.rdiag_detailed(self.__rdiag_id)
                elif v == "false":
                    _cantera.rdiag_brief(self.__rdiag_id)
            elif o == "dashed_color":
                _cantera.rdiag_setDashedColor(self.__rdiag_id, v)
            elif o == "bold_color":
                _cantera.rdiag_setBoldColor(self.__rdiag_id, v)
            elif o == "normal_color":
                _cantera.rdiag_setNormalColor(self.__rdiag_id, v)
            elif o == "scale":
                _cantera.rdiag_setScale(self.__rdiag_id, v)
            elif o == "dot_options":
                _cantera.rdiag_setDotOptions(self.__rdiag_id, v)
            elif o == "title":
                _cantera.rdiag_setTitle(self.__rdiag_id, v)
            elif o == "arrow_width":
                _cantera.rdiag_setArrowWidth(self.__rdiag_id, v)
            elif o == "threshold":
                _cantera.rdiag_setThreshold(self.__rdiag_id, v)
            elif o == "bold_threshold":
                _cantera.rdiag_setBoldThreshold(self.__rdiag_id, v)
            elif o == "normal_threshold":
                _cantera.rdiag_setNormalThreshold(self.__rdiag_id, v)
            elif o == "label_threshold":
                _cantera.rdiag_setLabelThreshold(self.__rdiag_id, v)
            elif o == "font":
                _cantera.rdiag_setFont(self.__rdiag_id, v)
            elif o == "flow_type":
                if v == "one_way":
                    _cantera.rdiag_setFlowType(self.__rdiag_id, 0)
                else:
                    _cantera.rdiag_setFlowType(self.__rdiag_id, 1)
            else:
                raise("unknown attribute "+o)

class PathBuilder:

    def __init__(self, kin, logfile=""):
        if logfile == "":
            logfile = "rxnpath.log"
        self.__rbuild_id = _cantera.rbuild_new()
        self.kin = kin
        _cantera.rbuild_init(self.__rbuild_id, logfile, kin.ckin)

    def __del__(self):
        _cantera.rbuild_del(self.__rbuild_id)

    def build(self, diagram = None, element = "C",
              dotfile = "rxnpaths.dot", format="dot"):
        if diagram == None:
            diagram = PathDiagram()
        _cantera.rbuild_build(self.__rbuild_id, self.kin.ckin, element,
                             "buildlog", diagram.id(), 1)
        if format == "dot":
            diagram.write(0, dotfile)
            diagram.write(1, "rp.txt")
        elif format == "plain":
            diagram.write(1, dotfile)


def write(g, el, file, d=None, format="dot"):
    b = PathBuilder(g)
    b.build(element = el, diagram = d, dotfile = file, format = format)


def view(url):
    import webbrowser
    webbrowser.open(url)


if __name__ == "__main__":
    from Cantera.gases import GRI30
    gas = GRI30()
    x = [1.0] * gas.nSpecies()
    gas.setState_TPX(1800.0, 1.01325e5, x)
    write(gas, 'C', 'c:/users/dgg/test.dot')
