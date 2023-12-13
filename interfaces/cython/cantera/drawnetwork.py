# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import importlib.metadata
from functools import wraps

_graphviz = None
def _import_graphviz():
    # defer import of graphviz
    global _graphviz
    if _graphviz is not None:
        return
    try:
        importlib.metadata.version("graphviz")
    except importlib.metadata.PackageNotFoundError:
        raise ImportError("This requires the graphviz package.")
    else:
        import graphviz as _graphviz


def needs_graphviz(func):
    # decorator function to load graphviz when needed
    @wraps(func)
    def inner(*args, **kwargs):
        if not _graphviz:
            _import_graphviz()
        return func(*args, **kwargs)

    return inner


@needs_graphviz
def draw_reactor(r, dot=None, print_state=False, species=None, **kwargs):
    """
    Draw `ReactorBase` object as ``graphviz`` ``dot`` node.
    The node is added to an existing ``dot`` graph if provided.
    Optionally include current reactor state in the node.

    :param r:
        `ReactorBase` object or subclass.
    :param dot:
        ``graphviz.graphs.BaseGraph`` object to which the reactor is added.
        If not provided, a new ``DiGraph`` is created. Defaults to ``None``.
    :param print_state:
        Whether state information of the reactor is printed into the node.
        Defaults to ``False``.
    :param species:
        If ``print_state`` is ``True``, define how species are to be printed.
        Options are ``'X'`` and ``'Y'`` for mole and mass fractions of all
        species, respectively, or an iterable that contains the desired species
        names as strings. Defaults to ``None``.
    :param **kwargs:
        Keyword options can contain ``graph_attr`` and general ``node_attr`` to
        be passed on to the ``graphviz`` functions to control the appearance of
        the graph and reactor node. ``node_attr`` defined in the reactor object
        itself have priority.
    :return:
        ``graphviz.graphs.BaseGraph`` object with reactor

    """
    if not dot:
        dot = _graphviz.Digraph(name=r.name,
                                graph_attr=kwargs.get("graph_attr"))

    # attributes defined in Reactor.node_attr overwrite default attributes
    node_attr = dict(kwargs.get("node_attr", {}), **r.node_attr)

    # include full reactor state in representation if desired
    if print_state:
        T_label = f"T (K)\\n{r.T:.2f}"
        P_label = f"P (bar)\\n{r.thermo.P*1e-5:.3f}"
        s_label = ""

        if species == "X":
            X = r.thermo.mole_fraction_dict()
            s_percents = "\\n".join([f"{s}: {v*100:.2f}" for s, v in X.items()])
            s_label += "X (%)\\n" + s_percents
        elif species == "Y":
            Y = r.thermo.mass_fraction_dict()
            s_percents = "\\n".join([f"{s}: {v*100:.2f}" for s, v in Y.items()])
            s_label += "Y (%)\\n" + s_percents
        else:
            # If individual species are provided as iterable of strings
            if species:
                X = r.thermo.mole_fraction_dict(threshold=-1)
                s_percents = "\\n".join([f"{s}: {X[s]*100:.2f}" for s in species])
                s_label += "X (%)\\n" + s_percents

        # For full state output, shape must be 'Mrecord'
        node_attr = {k:v for k,v in node_attr.items() if k != "shape"}
        dot.node(r.name, shape="Mrecord",
                         label="{"+ T_label +"|"+ P_label +"}"+"|"+ s_label,
                         xlabel=r.name,
                         **node_attr)

    else:
        dot.node(r.name, **node_attr)

    return dot


@needs_graphviz
def draw_reactor_net(n, **kwargs):
    """
    Draw `ReactorNet` object as ``graphviz.graphs.DiGraph``. Connecting flow
    controllers and walls are depicted as arrows.

    :param n:
        `ReactorNet` object
    :param **kwargs:
        Keyword options can contain ``graph_attr`` and general ``node_attr``,
        ``edge_attr``, ``heat_flow_attr``, and ``mass_flow_attr`` to be passed
        on to the ``graphviz`` functions to control the appearance of the
        graph, reactor nodes, and connection edges. ``node_attr`` and
        ``edge_attr`` defined in the objects themselves have priority.
    :return:
        ``graphviz.graphs.BaseGraph`` object with reactor net.

    """
    dot = _graphviz.Digraph(graph_attr=kwargs.get("graph_attr"))

    # collect elements as set to avoid duplicates
    reactors = set(n.reactors)
    connections = set()

    for r in reactors:
        draw_reactor(r, dot, **kwargs)
        connections.update(r.walls + r.inlets + r.outlets)

    # some Reactors or Reservoirs only exist as connecting nodes
    connected_reactors = get_connected_reactors(connections)

    # remove already drawn reactors and draw new reactors
    connected_reactors.difference_update(reactors)
    for r in connected_reactors:
        draw_reactor(r, dot, **kwargs)

    draw_connections(connections, dot, **kwargs)

    return dot


def get_connected_reactors(connections):
    """
    Collect and returned all connected reactors.

    :param connections:
        Iterable containing connections that are either subtypes of
        `FlowDevice` or `WallBase`.
    :return:
        A ``set`` containing the connected reactors and reservoir

    """

    connected_reactors = set()

    for c in connections:
        # this block detects whether c is a Wall or a FlowController
        try:
            r_in, r_out = c.upstream, c.downstream
        except AttributeError:
            r_in, r_out = c.left_reactor, c.right_reactor
        connected_reactors.update((r_in, r_out))

    return connected_reactors


@needs_graphviz
def draw_connections(connections, dot=None, **kwargs):
    """
    Draw connections between reactors and reservoirs. This includes flow
    controllers and walls.

    :param connections:
        Iterable containing connections that are either subtypes of
        `FlowDevice` or `WallBase`.
    :param dot:
        ``graphviz.graphs.BaseGraph`` object to which the connection is added.
        If not provided, a new ``DiGraph`` is created. Defaults to ``None``
    :param **kwargs:
        Keyword options can contain ``graph_attr`` and general ``node_attr``,
        ``edge_attr``, ``heat_flow_attr``, and ``mass_flow_attr`` to be passed
        on to the ``graphviz`` functions to control the appearance of the
        graph, reactor nodes, and connection edges. ``node_attr`` and
        ``edge_attr`` defined in the objects themselves have priority
    :return:
        A ``graphviz.graphs.BaseGraph`` object depicting the connections.

    """
    if not dot:
        dot = _graphviz.Digraph(graph_attr=kwargs.get("graph_attr"))
    # set default style for all connections and nodes if provided
    dot.edge_attr.update(kwargs.get("edge_attr", {}))
    dot.node_attr.update(kwargs.get("node_attr", {}))

    # retrieve default attributes for all mass flow and heat connections
    mass_flow_attr = kwargs.get("mass_flow_attr", {})
    heat_flow_attr = kwargs.get("heat_flow_attr", {})

    # using a while loop instead of iterating over all connections allows to
    # remove duplicate connections once they have been detected.
    connections = set(connections)
    while True:
        try:
            c = connections.pop()
        except KeyError:
            break

        # this block detects whether c is a Wall or a FlowController
        try:
            inflow, outflow = "upstream", "downstream"
            r_in, r_out = getattr(c, inflow), getattr(c, outflow)
            rate_attr = "mass_flow_rate"
            edge_attr = dict(mass_flow_attr, **c.edge_attr)
        except AttributeError:
            inflow, outflow = "left_reactor", "right_reactor"
            r_in, r_out = getattr(c, inflow), getattr(c, outflow)
            rate_attr = "heat_rate"
            edge_attr = {"color": "red", "style": "dashed",
                         **heat_flow_attr, **c.edge_attr}

        # find "duplicate" connections of the same kind that connect the same
        # two objects to eventually draw them as a single connection
        duplicates = set()
        inv_duplicates = set()
        for c_ in connections:
            try:
                if getattr(c_, inflow) is r_in and getattr(c_, outflow) is r_out:
                    duplicates.add(c_)
                elif getattr(c_, outflow) is r_in and getattr(c_, inflow) is r_out:
                    inv_duplicates.add(c_)
            except AttributeError:
                continue

        # remove duplicates from the set of the connections still to be drawn
        connections.difference_update(duplicates | inv_duplicates)

        # sum up heat rate/mass flow rate while considering the direction
        rate = (getattr(c, rate_attr)
                + sum(getattr(c, rate_attr) for c in duplicates)
                - sum(getattr(c, rate_attr) for c in inv_duplicates))

        # ensure arrow always indicates a positive flow
        if rate >= 0:
            inflow_name, outflow_name = r_in.name, r_out.name
        else:
            inflow_name, outflow_name = r_out.name, r_in.name
            rate *= -1

        if rate_attr == "mass_flow_rate":
            label = f"m = {rate:.2g} kg/s"
        elif rate_attr == "heat_rate":
            label = f"q = {rate:.2g} W"

        dot.edge(inflow_name, outflow_name, **{"label": label, **edge_attr})

    return dot
