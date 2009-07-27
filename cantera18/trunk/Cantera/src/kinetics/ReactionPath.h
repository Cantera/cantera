/**
 *  @file ReactionPath.h
 *
 *  Classes for reaction path analysis.
 *
 * $Author: dggoodwin $
 * $Revision: 1.1 $
 * $Date: 2007/05/04 14:27:23 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_RXNPATH_H
#define CT_RXNPATH_H

// Cantera includes 
#include "ct_defs.h"
#include "DenseMatrix.h"
#include "Group.h"
#include "Kinetics.h"

namespace Cantera {

    enum flow_t   { NetFlow, OneWayFlow };

    Group parseGroupString(std::string str, std::vector<std::string>& esyms);

    // forward references
    class Path; 

    /**
     *  Nodes in reaction path graphs.
     */
    class SpeciesNode {
    public:

        typedef std::vector<Path*> path_list;

        /// Default constructor
        SpeciesNode() : number(-1), name(""), value(0.0), 
                        visible(false), m_in(0.0), m_out(0.0) {}

        /// Destructor
        virtual ~SpeciesNode() {}

        // public attributes
        int         number;           ///<  Species number
        std::string      name;             ///<  Label on graph
        doublereal  value;            ///<  May be used to set node appearance
        bool        visible;          ///<  Visible on graph;


        // public methods

        /** 
         *  @name References.
         * Return a reference to a path object connecting this node
         *  to another node.
         */
        //@{
        Path*        path(int n)       { return m_paths[n]; }
        const Path*  path(int n) const { return m_paths[n]; }
        //@}


        /// Total number of paths to or from this node 
        int nPaths() const { return static_cast<int>(m_paths.size()); }

        /// add a path to or from this node
        void addPath(Path* path);

        double outflow() {return m_out;}
        double inflow() {return m_in;}
        double netOutflow() {return m_out - m_in;}
        
        void printPaths();

        
    protected:
        double m_in, m_out;
        path_list m_paths;
    };


        
    class Path {

    public:
 
        typedef std::map<int, doublereal> rxn_path_map;

        /**
         *  Constructor. Construct a one-way path from 
         *  \c begin to \c end.
         */
        Path(SpeciesNode* begin, SpeciesNode* end);

        /// Destructor
        virtual ~Path() {}

        void addReaction(int rxnNumber, doublereal value, std::string label = "");

        /// Upstream node.
        const SpeciesNode* begin() const { return m_a; }
        SpeciesNode* begin() { return m_a; }

        /// Downstream node.
        const SpeciesNode* end() const { return m_b; }
        SpeciesNode* end() { return m_b; }

        /**
         *  If \c n is one of the nodes this path connects, then
         *  the other node is returned. Otherwise zero is returned.
         */
        SpeciesNode* otherNode(SpeciesNode* n) { 
            return (n == m_a ? m_b : (n == m_b ? m_a : 0));
        }

        /// The total flow in this path
        doublereal flow() { return m_total; }
        void setFlow(doublereal v) { m_total = v; }

        ///  Number of reactions contributing to this path
        int nReactions() { 
			return static_cast<int>(m_rxn.size()); 
		}

        ///  Map from reaction number to flow from that reaction in this path.
        const rxn_path_map& reactionMap() { return m_rxn; }

        void writeLabel(std::ostream& s, doublereal threshold = 0.005);
        
    protected:

        std::map<std::string, doublereal> m_label;
        SpeciesNode *m_a, *m_b;
        rxn_path_map m_rxn;
        doublereal m_total;
    };


    /**
     *  Reaction path diagrams (graphs).
     */
    class ReactionPathDiagram {

    public:

        ReactionPathDiagram();
        
        virtual ~ReactionPathDiagram();

        /// The largest one-way flow value in any path
        doublereal maxFlow() { return m_flxmax; }

        /// The net flow from node \c k1 to node \c k2
        doublereal netFlow(int k1, int k2) { 
            return flow(k1, k2) - flow(k2, k1);
        }

        /// The one-way flow from node \c k1 to node \c k2
        doublereal flow(int k1, int k2) {
            return (m_paths[k1][k2] ? m_paths[k1][k2]->flow() : 0.0);
        }

        /// True if a node for species k exists
        bool hasNode(int k) {
            return (m_nodes[k] != 0);
        }

        void writeData(std::ostream& s);
        void exportToDot(std::ostream& s);
        void add(ReactionPathDiagram& d);
        SpeciesNode* node(int k) { return m_nodes[k]; }
        Path* path(int k1, int k2) { return m_paths[k1][k2]; }
        Path* path(int n) { return m_pathlist[n]; }
        int nPaths() { return static_cast<int>(m_pathlist.size()); }
        int nNodes() { return static_cast<int>(m_nodes.size()); }

        void addNode(int k, std::string nm, doublereal x = 0.0);

        void displayOnly(int k=-1) { m_local = k; }

        void linkNodes(int k1, int k2, int rxn, doublereal value,
            std::string legend = "");

        void include(std::string aaname) { m_include.push_back(aaname); }
        void exclude(std::string aaname) { m_exclude.push_back(aaname); }
        void include(std::vector<std::string>& names) { 
            int n = static_cast<int>(names.size());
            for (int i = 0; i < n; i++) m_include.push_back(names[i]);
        }
        void exclude(std::vector<std::string>& names) { 
            int n = static_cast<int>(names.size());
            for (int i = 0; i < n; i++) m_exclude.push_back(names[i]);
        }
        std::vector<std::string>& included() { return m_include; }
        std::vector<std::string>& excluded() { return m_exclude; }
        vector_int species();
        vector_int reactions();
        void findMajorPaths(doublereal threshold, int lda, doublereal* a);
        void setFont(std::string font) {
            m_font = font;
        }
        // public attributes

        std::string title;
        std::string bold_color;
        std::string normal_color;
        std::string dashed_color;
        std::string element;
        std::string m_font;
        doublereal threshold, 
            bold_min, dashed_max, label_min;
        doublereal x_size, y_size;
        std::string name, dot_options;
        flow_t flow_type;
        double scale;
        double arrow_width; 
        bool show_details;
        double arrow_hue;

    protected:

        doublereal                    m_flxmax;
        std::map<int, std::map<int, Path*> >    m_paths;
        std::map<int, SpeciesNode*>        m_nodes;
        std::vector<Path*>                 m_pathlist;
        std::vector<std::string>                m_include;
        std::vector<std::string>                m_exclude;
        vector_int                   m_speciesNumber;
        std::map<int, int>                 m_rxns;
        int                           m_local;
    };



    class ReactionPathBuilder {

    public:
        ReactionPathBuilder() {}
        virtual ~ReactionPathBuilder() {}
    
        int init(std::ostream& logfile, Kinetics& s);

        int build(Kinetics& s, std::string element, std::ostream& output, 
            ReactionPathDiagram& r, bool quiet=false);

        int findGroups(std::ostream& logfile, Kinetics& s);
 
        void writeGroup(std::ostream& out, const Group& g);

    protected:
        void findElements(Kinetics& kin);

        int m_nr;
        int m_ns;
        int m_nel;
        vector_fp m_ropf;
        vector_fp m_ropr;
        array_fp m_x;
        std::vector<vector_int> m_reac;
        std::vector<vector_int> m_prod;
        DenseMatrix m_elatoms;
        std::vector<std::vector<int> > m_groups;
        std::vector<Group> m_sgroup;
        std::vector<std::string> m_elementSymbols;
        //        std::map<int, int> m_warn;
        std::map<int, std::map<int, std::map<int, Group> > >  m_transfer;
        std::vector<bool> m_determinate;
        Array2D m_atoms;
        std::map<std::string,int> m_enamemap;
    };

}

#endif
