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
    /// Default constructor
    SpeciesNode() : number(npos), value(0.0),
        visible(false), m_in(0.0), m_out(0.0) {}

    /// Destructor
    virtual ~SpeciesNode() {}

    // public attributes
    size_t number; ///< Species number
    std::string name; ///< Label on graph
    doublereal value; ///< May be used to set node appearance
    bool visible; ///< Visible on graph;

    /**
     *  @name References.
     * Return a reference to a path object connecting this node
     *  to another node.
     */
    //@{
    Path* path(int n) {
        return m_paths[n];
    }
    const Path* path(int n) const {
        return m_paths[n];
    }
    //@}

    /// Total number of paths to or from this node
    int nPaths() const {
        return static_cast<int>(m_paths.size());
    }

    /// add a path to or from this node
    void addPath(Path* path);

    doublereal outflow() {
        return m_out;
    }
    doublereal inflow() {
        return m_in;
    }
    doublereal netOutflow() {
        return m_out - m_in;
    }

    void printPaths();

protected:
    doublereal m_in;
    doublereal m_out;
    std::vector<Path*> m_paths;
};


class Path
{
public:
    typedef std::map<size_t, doublereal> rxn_path_map;

    /**
     *  Constructor. Construct a one-way path from \c begin to \c end.
     */
    Path(SpeciesNode* begin, SpeciesNode* end);

    /// Destructor
    virtual ~Path() {}

    /**
     * Add a reaction to the path. Increment the flow from this reaction, the
     * total flow, and the flow associated with this label.
     */
    void addReaction(size_t rxnNumber, doublereal value,
                     const std::string& label = "");

    /// Upstream node.
    const SpeciesNode* begin() const {
        return m_a;
    }
    SpeciesNode* begin() {
        return m_a;
    }

    /// Downstream node.
    const SpeciesNode* end() const {
        return m_b;
    }
    SpeciesNode* end() {
        return m_b;
    }

    /**
     *  If \c n is one of the nodes this path connects, then
     *  the other node is returned. Otherwise zero is returned.
     */
    SpeciesNode* otherNode(SpeciesNode* n) {
        return (n == m_a ? m_b : (n == m_b ? m_a : 0));
    }

    /// The total flow in this path
    doublereal flow() {
        return m_total;
    }
    void setFlow(doublereal v) {
        m_total = v;
    }

    ///  Number of reactions contributing to this path
    int nReactions() {
        return static_cast<int>(m_rxn.size());
    }

    ///  Map from reaction number to flow from that reaction in this path.
    const rxn_path_map& reactionMap() {
        return m_rxn;
    }

    /**
     * Write the label for a path connecting two species, indicating
     * the percent of the total flow due to each reaction.
     */
    void writeLabel(std::ostream& s, doublereal threshold = 0.005);

protected:
    std::map<std::string, doublereal> m_label;
    SpeciesNode* m_a, *m_b;
    rxn_path_map m_rxn;
    doublereal m_total;
};


/**
 *  Reaction path diagrams (graphs).
 */
class ReactionPathDiagram
{
public:
    ReactionPathDiagram();

    /**
     * Destructor. Deletes all nodes and paths in the diagram.
     */
    virtual ~ReactionPathDiagram();

    /// The largest one-way flow value in any path
    doublereal maxFlow() {
        return m_flxmax;
    }

    /// The net flow from node \c k1 to node \c k2
    doublereal netFlow(size_t k1, size_t k2) {
        return flow(k1, k2) - flow(k2, k1);
    }

    /// The one-way flow from node \c k1 to node \c k2
    doublereal flow(size_t k1, size_t k2) {
        return (m_paths[k1][k2] ? m_paths[k1][k2]->flow() : 0.0);
    }

    /// True if a node for species k exists
    bool hasNode(size_t k) {
        return (m_nodes[k] != 0);
    }

    void writeData(std::ostream& s);

    /**
     *  Export the reaction path diagram. This method writes to stream
     *  \c s the commands for the 'dot' program in the \c GraphViz
     *  package from AT&T. (GraphViz may be downloaded from www.graphviz.org.)
     *
     *  To generate a postscript reaction path diagram from the output of this
     *  method saved in file paths.dot, for example, give the command:
     *  \code
     *  dot -Tps paths.dot > paths.ps
     *  \endcode
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

    void addNode(size_t k, const std::string& nm, doublereal x = 0.0);

    void displayOnly(size_t k=npos) {
        m_local = k;
    }

    void linkNodes(size_t k1, size_t k2, size_t rxn, doublereal value,
                   std::string legend = "");

    void include(const std::string& aaname) {
        m_include.push_back(aaname);
    }
    void exclude(const std::string& aaname) {
        m_exclude.push_back(aaname);
    }
    void include(std::vector<std::string>& names) {
        for (size_t i = 0; i < names.size(); i++) {
            m_include.push_back(names[i]);
        }
    }
    void exclude(std::vector<std::string>& names) {
        for (size_t i = 0; i < names.size(); i++) {
            m_exclude.push_back(names[i]);
        }
    }
    std::vector<std::string>& included() {
        return m_include;
    }
    std::vector<std::string>& excluded() {
        return m_exclude;
    }
    std::vector<size_t> species();
    vector_int reactions();
    void findMajorPaths(doublereal threshold, size_t lda, doublereal* a);
    void setFont(const std::string& font) {
        m_font = font;
    }
    // public attributes

    std::string title;
    std::string bold_color;
    std::string normal_color;
    std::string dashed_color;
    std::string element;
    std::string m_font;
    doublereal threshold, bold_min, dashed_max, label_min;
    doublereal x_size, y_size;
    std::string name, dot_options;
    flow_t flow_type;
    doublereal scale;
    doublereal arrow_width;
    bool show_details;
    doublereal arrow_hue;

protected:
    doublereal m_flxmax;
    std::map<size_t, std::map<size_t, Path*> > m_paths;
    std::map<size_t, SpeciesNode*> m_nodes;
    std::vector<Path*> m_pathlist;
    std::vector<std::string> m_include;
    std::vector<std::string> m_exclude;
    std::vector<size_t> m_speciesNumber;
    std::map<size_t, int> m_rxns;
    size_t m_local;
};


class ReactionPathBuilder
{
public:
    ReactionPathBuilder() {}
    virtual ~ReactionPathBuilder() {}

    int init(std::ostream& logfile, Kinetics& s);

    int build(Kinetics& s, const std::string& element, std::ostream& output,
              ReactionPathDiagram& r, bool quiet=false);

    //! Analyze a reaction to determine which reactants lead to which products.
    int findGroups(std::ostream& logfile, Kinetics& s);

protected:
    void findElements(Kinetics& kin);

    size_t m_nr;
    size_t m_ns;
    size_t m_nel;
    vector_fp m_ropf;
    vector_fp m_ropr;
    vector_fp m_x;
    std::vector<std::vector<size_t> > m_reac;
    std::vector<std::vector<size_t> > m_prod;
    DenseMatrix m_elatoms;
    std::vector<vector_int> m_groups;
    std::vector<Group> m_sgroup;
    std::vector<std::string> m_elementSymbols;

    //! m_transfer[reaction][reactant number][product number] where "reactant
    //! number" means the number of the reactant in the reaction equation, e.g.
    //! for "A+B -> C+D", "B" is reactant number 1 and "C" is product number 0.
    std::map<size_t, std::map<size_t, std::map<size_t, Group> > > m_transfer;

    std::vector<bool> m_determinate;
    Array2D m_atoms;
    std::map<std::string, size_t> m_enamemap;
};

}

#endif
