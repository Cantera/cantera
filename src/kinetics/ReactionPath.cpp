/**
 *  @file ReactionPath.cpp
 *  Implementation file for classes used in reaction path analysis.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionPath.h"
#include "cantera/kinetics/reaction_defs.h"

using namespace std;

namespace Cantera
{

void SpeciesNode::addPath(Path* path)
{
    m_paths.push_back(path);
    if (path->begin() == this) {
        m_out += path->flow();
    } else if (path->end() == this) {
        m_in += path->flow();
    } else {
        throw CanteraError("SpeciesNode::addPath", "path added to wrong node");
    }
}

void SpeciesNode::printPaths()
{
    for (size_t i = 0; i < m_paths.size(); i++) {
        cout << m_paths[i]->begin()->name << " -->  "
             << m_paths[i]->end()->name << ":   "
             << m_paths[i]->flow() << endl;
    }
}

Path::Path(SpeciesNode* begin, SpeciesNode* end)
    : m_a(begin), m_b(end), m_total(0.0)
{
    begin->addPath(this);
    end->addPath(this);
}

void Path::addReaction(size_t rxnNumber, doublereal value,
                       const string& label)
{
    m_rxn[rxnNumber] += value;
    m_total += value;
    if (label != "") {
        m_label[label] += value;
    }
}

void Path::writeLabel(ostream& s, doublereal threshold)
{
    if (m_label.size() == 0) {
        return;
    }
    doublereal v;
    for (const auto& label : m_label) {
        v = label.second/m_total;
        if (m_label.size() == 1) {
            s << label.first << "\\l";
        } else if (v > threshold) {
            s << label.first;
            int percent = int(100*v + 0.5);
            if (percent < 100) {
                s << " (" << percent << "%)\\l";
            } else {
                s << "\\l";
            }
        }
    }
}

ReactionPathDiagram::ReactionPathDiagram()
{
    name = "reaction_paths";
    m_flxmax = 0.0;
    bold_color = "blue";
    normal_color = "steelblue";
    dashed_color = "gray";
    dot_options = "center=1;";
    m_font = "Helvetica";
    bold_min = 0.2;
    dashed_max = 0.0;
    label_min = 0.0;
    threshold = 0.005;
    flow_type = NetFlow;
    scale = -1;
    x_size = -1.0;
    y_size = -1.0;
    arrow_width = -5.0;
    show_details = false;
    arrow_hue = 0.6666;
    title = "";
    m_local = npos;
}

ReactionPathDiagram::~ReactionPathDiagram()
{
    // delete the nodes
    for (const auto& node : m_nodes) {
        delete node.second;
    }

    // delete the paths
    size_t nn = nPaths();
    for (size_t n = 0; n < nn; n++) {
        delete m_pathlist[n];
    }
}

vector_int ReactionPathDiagram::reactions()
{
    doublereal flmax = 0.0;
    for (size_t i = 0; i < nPaths(); i++) {
        Path* p = path(i);
        flmax = std::max(p->flow(), flmax);
    }
    m_rxns.clear();
    for (size_t i = 0; i < nPaths(); i++) {
        for (const auto& rxn : path(i)->reactionMap()) {
            double flxratio = rxn.second/flmax;
            if (flxratio > threshold) {
                m_rxns[rxn.first] = 1;
            }
        }
    }
    vector_int r;
    for (const auto& rxn : m_rxns) {
        r.push_back(int(rxn.first));
    }
    return r;
}

void ReactionPathDiagram::add(ReactionPathDiagram& d)
{
    for (size_t n = 0; n < nPaths(); n++) {
        Path* p = path(n);
        size_t k1 = p->begin()->number;
        size_t k2 = p->end()->number;
        p->setFlow(p->flow() + d.flow(k1,k2));
    }
}

void ReactionPathDiagram::findMajorPaths(doublereal athreshold, size_t lda,
        doublereal* a)
{
    double netmax = 0.0;
    for (size_t n = 0; n < nNodes(); n++) {
        for (size_t m = n+1; m < nNodes(); m++) {
            size_t k1 = m_speciesNumber[n];
            size_t k2 = m_speciesNumber[m];
            double fl = fabs(netFlow(k1,k2));
            netmax = std::max(fl, netmax);
        }
    }
    for (size_t n = 0; n < nNodes(); n++) {
        for (size_t m = n+1; m < nNodes(); m++) {
            size_t k1 = m_speciesNumber[n];
            size_t k2 = m_speciesNumber[m];
            double fl = fabs(netFlow(k1,k2));
            if (fl > athreshold*netmax) {
                a[lda*k1 + k2] = 1;
            }
        }
    }
}

void ReactionPathDiagram::writeData(ostream& s)
{
    s << title << endl;
    for (size_t i1 = 0; i1 < nNodes(); i1++) {
        size_t k1 = m_speciesNumber[i1];
        s << m_nodes[k1]->name << " ";
    }
    s << endl;
    for (size_t i1 = 0; i1 < nNodes(); i1++) {
        size_t k1 = m_speciesNumber[i1];
        for (size_t i2 = i1+1; i2 < nNodes(); i2++) {
            size_t k2 = m_speciesNumber[i2];
            double f1 = flow(k1, k2);
            double f2 = flow(k2, k1);
            s << m_nodes[k1]->name << " " << m_nodes[k2]->name
              << " " << f1 << " " << -f2 << endl;
        }
    }
}

void ReactionPathDiagram::exportToDot(ostream& s)
{
    doublereal flmax = 0.0;
    s.precision(3);

    // a directed graph
    s << "digraph " << name << " {" << endl;

    // the graph will be no larger than x_size, y_size
    if (x_size > 0.0) {
        if (y_size < 0.0) {
            y_size = x_size;
        }
        s << "size = \""
          << x_size << ","
          << y_size << "\";"
          << endl;
    }

    if (dot_options != "") {
        s << dot_options << endl;
    }

    // draw paths representing net flows
    if (flow_type == NetFlow) {
        // if no scale was specified, normalize net flows by the maximum net
        // flow
        if (scale <= 0.0) {
            for (size_t i1 = 0; i1 < nNodes(); i1++) {
                size_t k1 = m_speciesNumber[i1];
                node(k1)->visible = false;
                for (size_t i2 = i1+1; i2 < nNodes(); i2++) {
                    size_t k2 = m_speciesNumber[i2];
                    double flx = netFlow(k1, k2);
                    if (flx < 0.0) {
                        flx = -flx;
                    }
                    flmax = std::max(flx, flmax);
                }
            }
        } else {
            flmax = scale;
        }
        flmax = std::max(flmax, 1e-10);

        // loop over all unique pairs of nodes
        for (size_t i1 = 0; i1 < nNodes(); i1++) {
            size_t k1 = m_speciesNumber[i1];
            for (size_t i2 = i1+1; i2 < nNodes(); i2++) {
                size_t k2 = m_speciesNumber[i2];
                double flx = netFlow(k1, k2);
                if (m_local != npos && k1 != m_local && k2 != m_local) {
                    flx = 0.0;
                }
                if (flx != 0.0) {
                    double flxratio;
                    size_t kbegin, kend;
                    // set beginning and end of the path based on the sign of
                    // the net flow
                    if (flx > 0.0) {
                        kbegin = k1;
                        kend = k2;
                        flxratio = flx/flmax;
                    } else {
                        kbegin = k2;
                        kend = k1;
                        flxratio = -flx/flmax;
                    }

                    // write out path specification if the net flow is greater
                    // than the threshold
                    if (flxratio >= threshold) {
                        // make nodes visible
                        node(kbegin)->visible = true;
                        node(kend)->visible = true;

                        s << "s" << kbegin << " -> s" << kend;
                        s <<  "[fontname=\""+m_font+"\", penwidth=";

                        if (arrow_width < 0) {
                            double lwidth = 1.0 - 4.0
                                     * log10(flxratio/threshold)/log10(threshold) + 1.0;
                            s << lwidth;
                            s << ", arrowsize="
                              <<  std::min(6.0, 0.5*lwidth);
                        } else {
                            s << arrow_width;
                            s << ", arrowsize=" << flxratio + 1;
                        }

                        doublereal hue = 0.7;
                        doublereal bright = 0.9;
                        s << ", color=" << "\"" << hue << ", "
                          << flxratio + 0.5
                          << ", " << bright << "\"" << endl;

                        if (flxratio > label_min) {
                            s << ", label=\" " << flxratio;
                            if (show_details) {
                                if (flow(kbegin, kend) > 0.0) {
                                    s << "\\l fwd: "
                                      << flow(kbegin, kend)/flmax << "\\l";
                                    path(kbegin, kend)->writeLabel(s);
                                }
                                if (flow(kend, kbegin) > 0.0) {
                                    s << " \\l rev: "
                                      << flow(kend,kbegin)/flmax << "\\l";
                                    path(kend, kbegin)->writeLabel(s);
                                }
                            }
                            s << "\"";
                        }
                        s << "];" << endl;
                    }
                }
            }
        }
    } else {
        if (scale < 0) {
            for (size_t i = 0; i < nPaths(); i++) {
                flmax = std::max(path(i)->flow(), flmax);
            }
        } else {
            flmax = scale;
        }

        for (size_t i = 0; i < nPaths(); i++) {
            Path* p = path(i);
            double flxratio = p->flow()/flmax;
            if (m_local != npos) {
                if (p->begin()->number != m_local
                        && p->end()->number != m_local) {
                    flxratio = 0.0;
                }
            }
            if (flxratio > threshold) {
                p->begin()->visible = true;
                p->end()->visible = true;
                s << "s" << p->begin()->number
                  << " -> s" << p->end()->number;

                if (arrow_width < 0) {
                    double lwidth = 1.0 - 4.0 * log10(flxratio/threshold)/log10(threshold)
                             + 1.0;
                    s <<  "[fontname=\""+m_font+"\", penwidth="
                      << lwidth;
                    s << ", arrowsize="
                      <<  std::min(6.0, 0.5*lwidth);
                } else {
                    s <<  ", penwidth="
                      <<  arrow_width;
                    s << ", arrowsize=" << flxratio + 1;
                }
                doublereal hue = 0.7;
                doublereal bright = 0.9;
                s << ", color=" << "\"" << hue << ", " << flxratio + 0.5
                  << ", " << bright << "\"" << endl;

                if (flxratio > label_min) {
                    s << ", label = \" " << flxratio;
                    if (show_details) {
                        s << "\\l";
                        p->writeLabel(s);
                    }
                    s << "\"";
                }
                s << "];" << endl;
            }
        }
    }
    s.precision(2);
    for (const auto& node : m_nodes) {
        if (node.second->visible) {
            s << "s" << node.first << " [ fontname=\""+m_font+"\", label=\"" << node.second->name
              << "\"];" << endl;
        }
    }
    s << " label = " << "\"" << "Scale = "
      << flmax << "\\l " << title << "\";" << endl;
    s  << " fontname = \""+m_font+"\";" << endl << "}" << endl;
}


void ReactionPathDiagram::addNode(size_t k, const string& nm, doublereal x)
{
    if (!m_nodes[k]) {
        m_nodes[k] = new SpeciesNode;
        m_nodes[k]->number = k;
        m_nodes[k]->name = nm;
        m_nodes[k]->value = x;
        m_speciesNumber.push_back(k);
    }
}

void ReactionPathDiagram::linkNodes(size_t k1, size_t k2, size_t rxn,
                                    doublereal value, string legend)
{
    Path* ff = m_paths[k1][k2];
    if (!ff) {
        ff= new Path(m_nodes[k1], m_nodes[k2]);
        m_paths[k1][k2] = ff;
        m_pathlist.push_back(ff);
    }
    ff->addReaction(rxn, value, legend);
    m_rxns[rxn] = 1;
    m_flxmax = std::max(ff->flow(), m_flxmax);
}

std::vector<size_t> ReactionPathDiagram::species()
{
    return m_speciesNumber;
}

int ReactionPathBuilder::findGroups(ostream& logfile, Kinetics& s)
{
    m_groups.resize(m_nr);
    for (size_t i = 0; i < m_nr; i++) { // loop over reactions
        logfile << endl << "Reaction " << i+1 << ": "
                << s.reactionString(i);

        if (m_determinate[i]) {
            logfile << " ... OK." << endl;
        } else if (m_reac[i].size() == 2 && m_prod[i].size() == 2) {
            // indices for the two reactants
            size_t kr0 = m_reac[i][0];
            size_t kr1 = m_reac[i][1];

            // indices for the two products
            size_t kp0 = m_prod[i][0];
            size_t kp1 = m_prod[i][1];

            // references to the Group objects representing the reactants
            const Group& r0 = m_sgroup[kr0];
            const Group& r1 = m_sgroup[kr1];
            const Group& p0 = m_sgroup[kp0];
            const Group& p1 = m_sgroup[kp1];

            const Group* group_a0=0, *group_b0=0, *group_c0=0,
                         *group_a1=0, *group_b1=0, *group_c1=0;
            Group b0 = p0 - r0;
            Group b1 = p1 - r0;
            if (b0.valid() && b1.valid()) {
                logfile << " ... ambiguous." << endl;
            } else if (!b0.valid() && !b1.valid()) {
                logfile << " ... cannot express as A + BC = AB + C" << endl;
            } else {
                logfile << endl;
            }

            if (b0.valid()) {
                if (b0.sign() > 0) {
                    group_a0 = &r0;
                    group_b0 = &b0;
                    group_c0 = &p1;
                    m_transfer[i][0][0] = r0;
                    m_transfer[i][1][0] = b0;
                    m_transfer[i][1][1] = p1;
                } else {
                    group_a0 = &r1;
                    group_c0 = &p0;
                    b0 *= -1;
                    group_b0 = &b0;
                    m_transfer[i][1][1] = r1;
                    m_transfer[i][0][1] = b0;
                    m_transfer[i][0][0] = p0;
                }
                logfile << "     ";
                group_a0->fmt(logfile, m_elementSymbols);
                logfile << " + ";
                group_b0->fmt(logfile,m_elementSymbols);
                group_c0->fmt(logfile, m_elementSymbols);
                logfile << " = ";
                group_a0->fmt(logfile, m_elementSymbols);
                group_b0->fmt(logfile, m_elementSymbols);
                logfile << " + ";
                group_c0->fmt(logfile, m_elementSymbols);
                if (b1.valid()) {
                    logfile << "   [<= default] " << endl;
                } else {
                    logfile << endl;
                }
            }

            if (b1.valid()) {
                if (b1.sign() > 0) {
                    group_a1 = &r0;
                    group_b1 = &b1;
                    group_c1 = &p0;
                    if (!b0.valid()) {
                        m_transfer[i][0][1] = r0;
                        m_transfer[i][1][1] = b0;
                        m_transfer[i][1][0] = p0;
                    }
                } else {
                    group_a1 = &r1;
                    group_c1 = &p1;
                    b1 *= -1;
                    group_b1 = &b1;
                    if (!b0.valid()) {
                        m_transfer[i][1][0] = r1;
                        m_transfer[i][0][0] = b0;
                        m_transfer[i][0][1] = p1;
                    }
                }
                logfile << "     ";
                group_a1->fmt(logfile, m_elementSymbols);
                logfile << " + ";
                group_b1->fmt(logfile, m_elementSymbols);
                group_c1->fmt(logfile, m_elementSymbols);
                logfile << " = ";
                group_a1->fmt(logfile, m_elementSymbols);
                group_b1->fmt(logfile, m_elementSymbols);
                logfile << " + ";
                group_c1->fmt(logfile, m_elementSymbols);
                logfile << endl;
            }
        } else {
            logfile << "... cannot parse. [ignored]" << endl;
        }
    }
    return 1;
}

void ReactionPathBuilder::findElements(Kinetics& kin)
{
    m_enamemap.clear();
    m_nel = 0;
    for (size_t i = 0; i < kin.nPhases(); i++) {
        ThermoPhase* p = &kin.thermo(i);
        // iterate over the elements in this phase
        for (size_t m = 0; m < p->nElements(); m++) {
            string ename = p->elementName(m);

            // if no entry is found for this element name, then it is a new
            // element. In this case, add the name to the list of names,
            // increment the element count, and add an entry to the
            // name->(index+1) map.
            if (m_enamemap.find(ename) == m_enamemap.end()) {
                m_enamemap[ename] = m_nel + 1;
                m_elementSymbols.push_back(ename);
                m_nel++;
            }
        }
    }
    m_atoms.resize(kin.nTotalSpecies(), m_nel, 0.0);
    // iterate over the elements
    for (size_t m = 0; m < m_nel; m++) {
        size_t k = 0;
        // iterate over the phases
        for (size_t ip = 0; ip < kin.nPhases(); ip++) {
            ThermoPhase* p = &kin.thermo(ip);
            size_t mlocal = p->elementIndex(m_elementSymbols[m]);
            for (size_t kp = 0; kp < p->nSpecies(); kp++) {
                if (mlocal != npos) {
                    m_atoms(k, m) = p->nAtoms(kp, mlocal);
                }
                k++;
            }
        }
    }
}

int ReactionPathBuilder::init(ostream& logfile, Kinetics& kin)
{
    m_transfer.clear();
    m_elementSymbols.clear();
    findElements(kin);
    m_ns = kin.nTotalSpecies();
    m_nr = kin.nReactions();

    // all reactants / products, even ones appearing on both sides of the
    // reaction
    vector<vector<size_t> > allProducts(m_nr);
    vector<vector<size_t> > allReactants(m_nr);
    for (size_t i = 0; i < m_nr; i++) {
        for (size_t k = 0; k < m_ns; k++) {
            for (int n = 0; n < kin.reactantStoichCoeff(k, i); n++) {
                allReactants[i].push_back(k);
            }
            for (int n = 0; n < kin.productStoichCoeff(k, i); n++) {
                allProducts[i].push_back(k);
            }
        }
    }

    // m_reac and m_prod exclude indices for species that appear on
    // both sides of the reaction, so that the diagram contains no loops.
    m_reac.resize(m_nr);
    m_prod.resize(m_nr);
    m_ropf.resize(m_nr);
    m_ropr.resize(m_nr);
    m_determinate.resize(m_nr);
    m_x.resize(m_ns); // not currently used ?
    m_elatoms.resize(m_nel, m_nr);

    for (size_t i = 0; i < m_nr; i++) {
        // construct the lists of reactant and product indices, not including
        // molecules that appear on both sides.
        m_reac[i].clear();
        m_prod[i].clear();
        map<size_t, int> net;
        size_t nr = allReactants[i].size();
        size_t np = allProducts[i].size();
        for (size_t ir = 0; ir < nr; ir++) {
            net[allReactants[i][ir]]--;
        }
        for (size_t ip = 0; ip < np; ip++) {
            net[allProducts[i][ip]]++;
        }

        for (size_t k = 0; k < m_ns; k++) {
            if (net[k] < 0) {
                size_t nmol = -net[k];
                for (size_t jr = 0; jr < nmol; jr++) {
                    m_reac[i].push_back(k);
                }
            } else if (net[k] > 0) {
                size_t nmol = net[k];
                for (size_t jp = 0; jp < nmol; jp++) {
                    m_prod[i].push_back(k);
                }
            }
        }

        size_t nrnet = m_reac[i].size();

        // compute number of atoms of each element in each reaction, excluding
        // molecules that appear on both sides of the reaction. We only need to
        // compute this for the reactants, since the elements are conserved.
        for (size_t n = 0; n < nrnet; n++) {
            size_t k = m_reac[i][n];
            for (size_t m = 0; m < m_nel; m++) {
                m_elatoms(m,i) += m_atoms(k,m);
            }
        }
    }

    // build species groups
    vector_int comp(m_nel);
    m_sgroup.resize(m_ns);
    for (size_t j = 0; j < m_ns; j++) {
        for (size_t m = 0; m < m_nel; m++) {
            comp[m] = int(m_atoms(j,m));
        }
        m_sgroup[j] = Group(comp);
    }

    // determine whether or not the reaction is "determinate", meaning that
    // there is no ambiguity about which reactant is the source for any element
    // in any product. This is false if more than one reactant contains a given
    // element, *and* more than one product contains the element. In this case,
    // additional information is needed to determine the partitioning of the
    // reactant atoms of that element among the products.
    for (size_t i = 0; i < m_nr; i++) {
        size_t nr = m_reac[i].size();
        size_t np = m_prod[i].size();
        m_determinate[i] = true;
        for (size_t m = 0; m < m_nel; m++) {
            int nar = 0;
            int nap = 0;
            for (size_t j = 0; j < nr; j++) {
                if (m_atoms(m_reac[i][j],m) > 0) {
                    nar++;
                }
            }
            for (size_t j = 0; j < np; j++) {
                if (m_atoms(m_prod[i][j],m) > 0) {
                    nap++;
                }
            }
            if (nar > 1 && nap > 1) {
                m_determinate[i] = false;
                break;
            }
        }
    }

    findGroups(logfile, kin);
    return 1;
}

string reactionLabel(size_t i, size_t kr, size_t nr,
                     const std::vector<size_t>& slist, const Kinetics& s)
{
    string label = "";
    for (size_t j = 0; j < nr; j++) {
        if (j != kr) {
            label += " + "+ s.kineticsSpeciesName(slist[j]);
        }
    }
    if (s.reactionType(i) == THREE_BODY_RXN) {
        label += " + M ";
    } else if (s.reactionType(i) == FALLOFF_RXN) {
        label += " (+ M)";
    }
    return label;
}

int ReactionPathBuilder::build(Kinetics& s, const string& element,
                               ostream& output, ReactionPathDiagram& r, bool quiet)
{
    map<size_t, int> warn;
    doublereal threshold = 0.0;
    size_t m = m_enamemap[element]-1;
    r.element = element;
    if (m == npos) {
        return -1;
    }

    s.getFwdRatesOfProgress(m_ropf.data());
    s.getRevRatesOfProgress(m_ropr.data());

    // species explicitly included or excluded
    vector<string>& in_nodes = r.included();
    vector<string>& out_nodes = r.excluded();

    vector_int status(s.nTotalSpecies(), 0);
    for (size_t ni = 0; ni < in_nodes.size(); ni++) {
        status[s.kineticsSpeciesIndex(in_nodes[ni])] = 1;
    }
    for (size_t ne = 0; ne < out_nodes.size(); ne++) {
        status[s.kineticsSpeciesIndex(out_nodes[ne])] = -1;
    }

    for (size_t i = 0; i < m_nr; i++) {
        double ropf = m_ropf[i];
        double ropr = m_ropr[i];

        // loop over reactions involving element m
        if (m_elatoms(m, i) > 0) {
            size_t nr = m_reac[i].size();
            size_t np = m_prod[i].size();

            for (size_t kr = 0; kr < nr; kr++) {
                size_t kkr = m_reac[i][kr];
                string fwdlabel = reactionLabel(i, kr, nr, m_reac[i], s);

                for (size_t kp = 0; kp < np; kp++) {
                    size_t kkp = m_prod[i][kp];
                    string revlabel = "";
                    for (size_t j = 0; j < np; j++) {
                        if (j != kp) {
                            revlabel += " + "+ s.kineticsSpeciesName(m_prod[i][j]);
                        }
                    }
                    if (s.reactionType(i) == THREE_BODY_RXN) {
                        revlabel += " + M ";
                    } else if (s.reactionType(i) == FALLOFF_RXN) {
                        revlabel += " (+ M)";
                    }

                    // calculate the flow only for pairs that are not the same
                    // species, both contain atoms of element m, and both are
                    // allowed to appear in the diagram
                    if ((kkr != kkp) && (m_atoms(kkr,m) > 0
                                         && m_atoms(kkp,m) > 0)
                            && status[kkr] >= 0 && status[kkp] >= 0) {
                        // if neither species contains the full number of atoms
                        // of element m in the reaction, then we must consider
                        // the type of reaction to determine which reactant
                        // species was the source of a given m-atom in the
                        // product
                        double f;
                        if ((m_atoms(kkp,m) < m_elatoms(m, i)) &&
                                (m_atoms(kkr,m) < m_elatoms(m, i))) {
                            map<size_t, map<size_t, Group> >& g = m_transfer[i];
                            if (g.empty()) {
                                if (!warn[i] && !quiet) {
                                    output << endl;
                                    output << "*************** REACTION IGNORED ***************" << endl;
                                    output << "Warning: no rule to determine partitioning of " << element
                                           << endl << " in reaction " << s.reactionString(i) << "." << endl
                                           << "*************** REACTION IGNORED **************" << endl;
                                    output << endl;
                                    warn[i] = 1;
                                }
                                f = 0.0;
                            } else {
                                if (!g[kr][kp]) {
                                    f = 0.0;
                                } else {
                                    f = g[kr][kp].nAtoms(m);
                                }
                            }
                        } else {
                            // no ambiguity about where the m-atoms come from or
                            // go to. Either all reactant m atoms end up in one
                            // product, or only one reactant contains all the
                            // m-atoms. In either case, the number of atoms
                            // transferred is given by the same expression.
                            f = m_atoms(kkp,m) * m_atoms(kkr,m) / m_elatoms(m, i);
                        }

                        double fwd = ropf*f;
                        double rev = ropr*f;
                        bool force_incl = ((status[kkr] == 1) || (status[kkp] == 1));

                        bool fwd_incl = ((fwd > threshold) ||
                                         (fwd > 0.0 && force_incl));
                        bool rev_incl = ((rev > threshold) ||
                                         (rev > 0.0 && force_incl));
                        if (fwd_incl || rev_incl) {
                            if (!r.hasNode(kkr)) {
                                r.addNode(kkr, s.kineticsSpeciesName(kkr), m_x[kkr]);
                            }
                            if (!r.hasNode(kkp)) {
                                r.addNode(kkp, s.kineticsSpeciesName(kkp), m_x[kkp]);
                            }
                        }
                        if (fwd_incl) {
                            r.linkNodes(kkr, kkp, int(i), fwd, fwdlabel);
                        }
                        if (rev_incl) {
                            r.linkNodes(kkp, kkr, -int(i), rev, revlabel);
                        }
                    }
                }
            }
        }
    }
    return 1;
}

}
