from pathlib import Path
from sphinx_gallery.scrapers import figure_rst

# Provide options to examples that only generate plots if an option is specified
class ResetArgv:
    wants_plot = {
        "adiabatic.py",
        "premixed_counterflow_twin_flame.py",
        "piston.py",
        "reactor1.py",
        "reactor2.py",
        "sensitivity1.py",
    }
    def __repr__(self):
        return 'ResetArgv'

    def __call__(self, sphinx_gallery_conf, script_vars):
        if Path(script_vars['src_file']).name in self.wants_plot:
            return ['--plot']
        else:
            return []


class GraphvizScraper():
    """
    Capture Graphviz objects that are assigned to variables in the global namespace.
    """
    def __init__(self):
        # IDs of graphviz objects that have already been seen and processed
        self.processed = set()

    def __repr__(self):
        return 'GraphvizScraper'

    def __call__(self, block, block_vars, gallery_conf):
        import graphviz
        # We use a list to collect references to image names
        image_names = list()

        # The `image_path_iterator` is created by Sphinx-Gallery, it will yield
        # a path to a file name that adheres to Sphinx-Gallery naming convention.
        image_path_iterator = block_vars['image_path_iterator']

        graph_types = (graphviz.Source, graphviz.graphs.Digraph, graphviz.graphs.Graph)
        for obj in block_vars["example_globals"].values():
            if isinstance(obj, graph_types) and id(obj) not in self.processed:
                self.processed.add(id(obj))
                image_path = Path(next(image_path_iterator)).with_suffix(".svg")
                obj.format = "svg"
                obj.render(image_path.with_suffix(""))
                image_names.append(image_path)

        # Use the `figure_rst` helper function to generate the reST for this
        # code block's figures.
        return figure_rst(image_names, gallery_conf['src_dir'])
