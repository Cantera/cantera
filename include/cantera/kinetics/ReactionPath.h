/**
 *  @file ReactionPath.h
 *  Classes for reaction path analysis.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RXNPATH_H
#define CT_RXNPATH_H

#include "cantera/numerics/DenseMatrix.h"
#include "Group.h"
#include "Kinetics.h"

namespace Cantera
{
enum flow_t { NetFlow, OneWayFlow };

// forward references
class Path;

/**
 *  Nodes in reaction path graphs.
 */
class SpeciesNode
{
public:
    //! Default constructor
    SpeciesNode() = default;

    //! Destructor
    virtual ~SpeciesNode() = default;

    // public attributes
    size_t number = npos; //!< Species number
    string name; //!< Label on graph
    double value = 0.0; //!< May be used to set node appearance
    bool visible = false; //!< Visible on graph;

    //! @name References
    //!
    //! Return a reference to a path object connecting this node
    //! to another node.
    //! @{
    Path* path(int n) {
        return m_paths[n];
    }
    const Path* path(int n) const {
        return m_paths[n];
    }
    //! @}

    //! Total number of paths to or from this node
    int nPaths() const {
        return static_cast<int>(m_paths.size());
    }

    //! add a path to or from this node
    void addPath(Path* path);

    double outflow() {
        return m_out;
    }
    double inflow() {
        return m_in;
    }
    double netOutflow() {
        return m_out - m_in;
    }

    void printPaths();

protected:
    double m_in = 0.0;
    double m_out = 0.0;
    vector<Path*> m_paths;
};


class Path
{
public:
    typedef map<size_t, double> rxn_path_map;

    /**
     *  Constructor. Construct a one-way path from @c begin to @c end.
     */
    Path(SpeciesNode* begin, SpeciesNode* end);

    //! Destructor
    virtual ~Path() {}

    /**
     * Add a reaction to the path. Increment the flow from this reaction, the
     * total flow, and the flow associated with this label.
     */
    void addReaction(size_t rxnNumber, double value, const string& label = "");

    //! Upstream node.
    const SpeciesNode* begin() const {
        return m_a;
    }
    SpeciesNode* begin() {
        return m_a;
    }

    //! Downstream node.
    const SpeciesNode* end() const {
        return m_b;
    }
    SpeciesNode* end() {
        return m_b;
    }

    /**
     *  If @c n is one of the nodes this path connects, then
     *  the other node is returned. Otherwise zero is returned.
     */
    SpeciesNode* otherNode(SpeciesNode* n) {
        return (n == m_a ? m_b : (n == m_b ? m_a : 0));
    }

    //! The total flow in this path
    double flow() {
        return m_total;
    }
    void setFlow(double v) {
        m_total = v;
    }

    //!  Number of reactions contributing to this path
    int nReactions() {
        return static_cast<int>(m_rxn.size());
    }

    //!  Map from reaction number to flow from that reaction in this path.
    const rxn_path_map& reactionMap() {
        return m_rxn;
    }

    /**
     * Write the label for a path connecting two species, indicating
     * the percent of the total flow due to each reaction.
     */
    void writeLabel(std::ostream& s, double threshold = 0.005);

protected:
    map<string, double> m_label;
    SpeciesNode* m_a, *m_b;
    rxn_path_map m_rxn;
    double m_total = 0.0;
};


/**
 *  Reaction path diagrams (graphs).
 */
class ReactionPathDiagram
{
public:
    ReactionPathDiagram() = default;

    /**
     * Destructor. Deletes all nodes and paths in the diagram.
     */
    virtual ~ReactionPathDiagram();

    //! The largest one-way flow value in any path
    double maxFlow() {
        return m_flxmax;
    }

    //! The net flow from node @c k1 to node @c k2
    double netFlow(size_t k1, size_t k2) {
        return flow(k1, k2) - flow(k2, k1);
    }

    //! The one-way flow from node @c k1 to node @c k2
    double flow(size_t k1, size_t k2) {
        return (m_paths[k1][k2] ? m_paths[k1][k2]->flow() : 0.0);
    }

    //! True if a node for species k exists
    bool hasNode(size_t k) {
        return (m_nodes[k] != 0);
    }

    void writeData(std::ostream& s);

    /**
     *  Export the reaction path diagram. This method writes to stream
     *  @c s the commands for the 'dot' program in the @c GraphViz
     *  package from AT&T. (GraphViz may be downloaded from www.graphviz.org.)
     *
     *  To generate a postscript reaction path diagram from the output of this
     *  method saved in file paths.dot, for example, give the command:
     *  @code
     *  dot -Tps paths.dot > paths.ps
     *  @endcode
     *  To generate a GIF image, replace -Tps with -Tgif
     */
    void exportToDot(std::ostream& s);

    void add(ReactionPathDiagram& d);
    SpeciesNode* node(size_t k) {
        return m_nodes[k];
    }
    Path* path(size_t k1, size_t k2) {
        return m_paths[k1][k2];
    }
    Path* path(size_t n) {
        return m_pathlist[n];
    }
    size_t nPaths() {
        return m_pathlist.size();
    }
    size_t nNodes() {
        return m_nodes.size();
    }

    void addNode(size_t k, const string& nm, double x = 0.0);

    void displayOnly(size_t k=npos) {
        m_local = k;
    }

    void linkNodes(size_t k1, size_t k2, size_t rxn, double value, string legend = "");

    void include(const string& aaname) {
        m_include.push_back(aaname);
    }
    void exclude(const string& aaname) {
        m_exclude.push_back(aaname);
    }
    void include(vector<string>& names) {
        for (size_t i = 0; i < names.size(); i++) {
            m_include.push_back(names[i]);
        }
    }
    void exclude(vector<string>& names) {
        for (size_t i = 0; i < names.size(); i++) {
            m_exclude.push_back(names[i]);
        }
    }
    vector<string>& included() {
        return m_include;
    }
    vector<string>& excluded() {
        return m_exclude;
    }
    vector<size_t> species();
    vector<int> reactions();
    void findMajorPaths(double threshold, size_t lda, double* a);
    void setFont(const string& font) {
        m_font = font;
    }
    // public attributes

    string title;
    string bold_color = "blue";
    string normal_color = "steelblue";
    string dashed_color = "gray";
    string element;
    string m_font = "Helvetica";
    double threshold = 0.005;
    double bold_min = 0.2;
    double dashed_max = 0.0;
    double label_min = 0.0;
    double x_size = -1.0;
    double y_size = -1.0;
    string name = "reaction_paths";
    string dot_options = "center=1;";
    flow_t flow_type = NetFlow;
    double scale = -1;
    double arrow_width = -5.0;
    bool show_details = false;
    double arrow_hue = 0.6666;

protected:
    double m_flxmax = 0.0;
    map<size_t, map<size_t, Path*>> m_paths;

    //! map of species index to SpeciesNode
    map<size_t, SpeciesNode*> m_nodes;
    vector<Path*> m_pathlist;
    vector<string> m_include;
    vector<string> m_exclude;
    vector<size_t> m_speciesNumber;

    //! Indices of reactions that are included in the diagram
    set<size_t> m_rxns;
    size_t m_local = npos;
};


class ReactionPathBuilder
{
public:
    ReactionPathBuilder() = default;
    virtual ~ReactionPathBuilder() = default;

    int init(std::ostream& logfile, Kinetics& s);

    int build(Kinetics& s, const string& element, std::ostream& output,
              ReactionPathDiagram& r, bool quiet=false);

    //! Analyze a reaction to determine which reactants lead to which products.
    int findGroups(std::ostream& logfile, Kinetics& s);

protected:
    void findElements(Kinetics& kin);

    size_t m_nr;
    size_t m_ns;
    size_t m_nel;
    vector<double> m_ropf;
    vector<double> m_ropr;
    vector<double> m_x;
    vector<vector<size_t>> m_reac;
    vector<vector<size_t>> m_prod;
    DenseMatrix m_elatoms;
    vector<vector<int>> m_groups;
    vector<Group> m_sgroup;
    vector<string> m_elementSymbols;

    //! m_transfer[reaction][reactant number][product number] where "reactant
    //! number" means the number of the reactant in the reaction equation. For example,
    //! for "A+B -> C+D", "B" is reactant number 1 and "C" is product number 0.
    map<size_t, map<size_t, map<size_t, Group>>> m_transfer;

    vector<bool> m_determinate;
    Array2D m_atoms;
    map<string, size_t> m_enamemap;
};

}

#endif
