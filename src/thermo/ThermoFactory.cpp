/**
 *  @file ThermoFactory.cpp
 *     Definitions for the factory class that can create known ThermoPhase objects
 *     (see \ref thermoprops and class \link Cantera::ThermoFactory ThermoFactory\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ThermoFactory.h"

#include "cantera/thermo/Species.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/PDSSFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/PlasmaPhase.h"

#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/MaskellSolidSolnPhase.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/thermo/RedlichKisterVPSSTP.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/thermo/PureFluidPhase.h"
#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/PengRobinson.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/CoverageDependentSurfPhase.h"
#include "cantera/thermo/EdgePhase.h"
#include "cantera/thermo/MetalPhase.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/thermo/LatticePhase.h"
#include "cantera/thermo/HMWSoln.h"
#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/thermo/IdealSolnGasVPSS.h"
#include "cantera/thermo/WaterSSTP.h"
#include "cantera/thermo/BinarySolutionTabulatedThermo.h"
#include "cantera/base/stringUtils.h"

#include <boost/algorithm/string.hpp>

namespace Cantera
{

ThermoFactory* ThermoFactory::s_factory = 0;
std::mutex ThermoFactory::thermo_mutex;

ThermoFactory::ThermoFactory()
{
    reg("none", []() { return new ThermoPhase(); });
    addDeprecatedAlias("none", "ThermoPhase");
    addDeprecatedAlias("none", "None");
    reg("ideal-gas", []() { return new IdealGasPhase(); });
    addDeprecatedAlias("ideal-gas", "IdealGas");
    reg("plasma", []() { return new PlasmaPhase(); });
    reg("ideal-surface", []() { return new SurfPhase(); });
    addDeprecatedAlias("ideal-surface", "Surface");
    addDeprecatedAlias("ideal-surface", "Surf");
    reg("coverage-dependent-surface", []() { return new CoverageDependentSurfPhase(); });
    reg("edge", []() { return new EdgePhase(); });
    addDeprecatedAlias("edge", "Edge");
    reg("electron-cloud", []() { return new MetalPhase(); });
    addDeprecatedAlias("electron-cloud", "Metal");
    reg("fixed-stoichiometry", []() { return new StoichSubstance(); });
    addDeprecatedAlias("fixed-stoichiometry", "StoichSubstance");
    reg("pure-fluid", []() { return new PureFluidPhase(); });
    addDeprecatedAlias("pure-fluid", "PureFluid");
    reg("compound-lattice", []() { return new LatticeSolidPhase(); });
    addDeprecatedAlias("compound-lattice", "LatticeSolid");
    reg("lattice", []() { return new LatticePhase(); });
    addDeprecatedAlias("lattice", "Lattice");
    reg("HMW-electrolyte", []() { return new HMWSoln(); });
    addDeprecatedAlias("HMW-electrolyte", "HMW");
    addDeprecatedAlias("HMW-electrolyte", "HMWSoln");
    reg("ideal-condensed", []() { return new IdealSolidSolnPhase(); });
    addDeprecatedAlias("ideal-condensed", "IdealSolidSolution");
    addDeprecatedAlias("ideal-condensed", "IdealSolidSoln");
    reg("Debye-Huckel", []() { return new DebyeHuckel(); });
    addDeprecatedAlias("Debye-Huckel", "DebyeHuckel");
    reg("ideal-molal-solution", []() { return new IdealMolalSoln(); });
    addDeprecatedAlias("ideal-molal-solution", "IdealMolalSolution");
    addDeprecatedAlias("ideal-molal-solution", "IdealMolalSoln");
    reg("ideal-solution-VPSS", []() { return new IdealSolnGasVPSS(); });
    reg("ideal-gas-VPSS", []() { return new IdealSolnGasVPSS(); });
    addDeprecatedAlias("ideal-solution-VPSS", "IdealSolnVPSS");
    addDeprecatedAlias("ideal-solution-VPSS", "IdealSolnGas");
    addDeprecatedAlias("ideal-gas-VPSS", "IdealGasVPSS");
    reg("Margules", []() { return new MargulesVPSSTP(); });
    reg("ions-from-neutral-molecule", []() { return new IonsFromNeutralVPSSTP(); });
    addAlias("ions-from-neutral-molecule", "IonsFromNeutralMolecule");
    addAlias("ions-from-neutral-molecule", "IonsFromNeutral");
    reg("Redlich-Kister", []() { return new RedlichKisterVPSSTP(); });
    addDeprecatedAlias("Redlich-Kister", "RedlichKister");
    reg("Redlich-Kwong", []() { return new RedlichKwongMFTP(); });
    addDeprecatedAlias("Redlich-Kwong", "RedlichKwongMFTP");
    addDeprecatedAlias("Redlich-Kwong", "RedlichKwong");
    reg("Maskell-solid-solution", []() { return new MaskellSolidSolnPhase(); });
    addAlias("Maskell-solid-solution", "MaskellSolidSolnPhase");
    addAlias("Maskell-solid-solution", "MaskellSolidsoln");
    reg("liquid-water-IAPWS95", []() { return new WaterSSTP(); });
    addDeprecatedAlias("liquid-water-IAPWS95", "PureLiquidWater");
    addDeprecatedAlias("liquid-water-IAPWS95", "Water");
    reg("binary-solution-tabulated", []() { return new BinarySolutionTabulatedThermo(); });
    addDeprecatedAlias("binary-solution-tabulated", "BinarySolutionTabulatedThermo");
    reg("Peng-Robinson", []() { return new PengRobinson(); });
}

ThermoFactory* ThermoFactory::factory()
{
    std::unique_lock<std::mutex> lock(thermo_mutex);
    if (!s_factory) {
        s_factory = new ThermoFactory;
    }
    return s_factory;
}

void ThermoFactory::deleteFactory()
{
    std::unique_lock<std::mutex> lock(thermo_mutex);
    delete s_factory;
    s_factory = 0;
}

ThermoPhase* ThermoFactory::newThermoPhase(const std::string& model)
{
    warn_deprecated("newThermoPhase",
        "To be removed after Cantera 3.0; superseded by newThermoModel.");
    return create(model);
}

ThermoPhase* newThermoPhase(const string& model)
{
    warn_deprecated("newThermoPhase",
        "To be removed after Cantera 3.0; superseded by newThermo.");
    return ThermoFactory::factory()->create(model);
}

shared_ptr<ThermoPhase> newThermoModel(const string& model)
{
    shared_ptr<ThermoPhase> tptr(ThermoFactory::factory()->create(model));
    return tptr;
}

shared_ptr<ThermoPhase> newThermo(const AnyMap& phaseNode, const AnyMap& rootNode)
{
    if (!phaseNode.hasKey("kinetics") && phaseNode.hasKey("reactions")) {
        throw InputFileError("newThermo", phaseNode["reactions"],
            "Phase entry includes a 'reactions' field but does not "
            "specify a kinetics model.");
    }
    string model = phaseNode["thermo"].asString();
    shared_ptr<ThermoPhase> t = newThermoModel(model);
    setupPhase(*t, phaseNode, rootNode);
    return t;
}

shared_ptr<ThermoPhase> newThermo(const string& infile, const string& id)
{
    size_t dot = infile.find_last_of(".");
    if (dot == npos) {
        // @todo Remove after Cantera 3.0
        warn_deprecated("newThermo",
            "Changed in Cantera 3.0. Replaced by newThermoModel.\n");
        newThermoModel(infile);
    }
    string extension;
    extension = toLowerCopy(infile.substr(dot+1));
    string id_ = id;
    if (id == "-") {
        id_ = "";
    }
    if (extension == "cti" || extension == "xml") {
        throw CanteraError("newThermo",
                           "The CTI and XML formats are no longer supported.");
    }

    AnyMap root = AnyMap::fromYamlFile(infile);
    AnyMap& phase = root["phases"].getMapWhere("name", id_);
    return newThermo(phase, root);
}

unique_ptr<ThermoPhase> newPhase(const AnyMap& phaseNode, const AnyMap& rootNode)
{
    warn_deprecated("newPhase",
        "To be removed after Cantera 3.0; superseded by\n"
        "newThermo(const AnyMap&, const AnyMap&).");
    if (!phaseNode.hasKey("kinetics") && phaseNode.hasKey("reactions")) {
        throw InputFileError("newPhase", phaseNode["reactions"],
            "Phase entry includes a 'reactions' field but does not "
            "specify a kinetics model.");
    }
    unique_ptr<ThermoPhase> t(newThermoPhase(phaseNode["thermo"].asString()));
    setupPhase(*t, phaseNode, rootNode);
    return t;
}

ThermoPhase* newPhase(const std::string& infile, std::string id)
{
    warn_deprecated("newPhase",
        "To be removed after Cantera 3.0; superseded by\n"
        "newThermo(const std::string&, const std::string&).");
    return newThermo(infile, id).get();
}

void addDefaultElements(ThermoPhase& thermo, const vector<string>& element_names) {
    for (const auto& symbol : element_names) {
        thermo.addElement(symbol);
    }
}

void addElements(ThermoPhase& thermo, const vector<string>& element_names,
                 const AnyValue& elements, bool allow_default)
{
    const auto& local_elements = elements.asMap("symbol");
    for (const auto& symbol : element_names) {
        if (local_elements.count(symbol)) {
            auto& element = *local_elements.at(symbol);
            double weight = element["atomic-weight"].asDouble();
            long int number = element.getInt("atomic-number", 0);
            double e298 = element.getDouble("entropy298", ENTROPY298_UNKNOWN);
            thermo.addElement(symbol, weight, number, e298);
        } else if (allow_default) {
            thermo.addElement(symbol);
        } else {
            throw InputFileError("addElements", elements,
                                 "Element '{}' not found", symbol);
        }
    }
}

void addSpecies(ThermoPhase& thermo, const AnyValue& names, const AnyValue& species)
{
    if (names.is<vector<string>>()) {
        // 'names' is a list of species names which should be found in 'species'
        const auto& species_nodes = species.asMap("name");
        for (const auto& name : names.asVector<string>()) {
            if (species_nodes.count(name)) {
                thermo.addSpecies(newSpecies(*species_nodes.at(name)));
            } else {
                throw InputFileError("addSpecies", names, species,
                    "Could not find a species named '{}'.", name);
            }
        }
    } else if (names == "all") {
        // The keyword 'all' means to add all species from this source
        for (const auto& item : species.asVector<AnyMap>()) {
            thermo.addSpecies(newSpecies(item));
        }
    } else {
        throw InputFileError("addSpecies", names,
            "Could not parse species declaration of type '{}'", names.type_str());
    }
}

void setupPhase(ThermoPhase& thermo, const AnyMap& phaseNode, const AnyMap& rootNode)
{
    thermo.setName(phaseNode["name"].asString());

    if (phaseNode.hasKey("deprecated")) {
        string msg = phaseNode["deprecated"].asString();
        string filename = phaseNode.getString("__file__",
            rootNode.getString("__file__", "unknown file"));
        string method = fmt::format("{}/{}", filename, phaseNode["name"].asString());
        warn_deprecated(method, phaseNode, msg);
    }

    // Add elements
    if (phaseNode.hasKey("elements")) {
        if (phaseNode.getBool("skip-undeclared-elements", false)) {
            thermo.ignoreUndefinedElements();
        } else {
            thermo.throwUndefinedElements();
        }

        if (phaseNode["elements"].is<vector<string>>()) {
            // 'elements' is a list of element symbols
            if (rootNode.hasKey("elements")) {
                addElements(thermo, phaseNode["elements"].asVector<string>(),
                            rootNode["elements"], true);
            } else {
                addDefaultElements(thermo, phaseNode["elements"].asVector<string>());
            }
        } else if (phaseNode["elements"].is<vector<AnyMap>>()) {
            // Each item in 'elements' is a map with one item, where the key is
            // a section in this file or another YAML file, and the value is a
            // list of element symbols to read from that section
            for (const auto& elemNode : phaseNode["elements"].asVector<AnyMap>()) {
                const string& source = elemNode.begin()->first;
                const auto& names = elemNode.begin()->second.asVector<string>();
                const auto& slash = boost::ifind_last(source, "/");
                if (slash) {
                    std::string fileName(source.begin(), slash.begin());
                    std::string node(slash.end(), source.end());
                    const AnyMap elements = AnyMap::fromYamlFile(fileName,
                        rootNode.getString("__file__", ""));
                    addElements(thermo, names, elements.at(node), false);
                } else if (rootNode.hasKey(source)) {
                    addElements(thermo, names, rootNode.at(source), false);
                } else if (source == "default") {
                    addDefaultElements(thermo, names);
                } else {
                    throw InputFileError("setupPhase", elemNode,
                        "Could not find elements section named '{}'", source);
                }
            }
        } else {
            throw InputFileError("setupPhase", phaseNode["elements"],
                "Could not parse elements declaration of type '{}'",
                phaseNode["elements"].type_str());
        }
    } else {
        // If no elements list is provided, just add elements as-needed from the
        // default list.
        thermo.addUndefinedElements();
    }

    // Add species
    if (phaseNode.hasKey("species")) {
        if (phaseNode["species"].is<vector<string>>()) {
            // 'species' is a list of species names to be added from the current
            // file's 'species' section
            addSpecies(thermo, phaseNode["species"], rootNode["species"]);
        } else if (phaseNode["species"].is<string>()) {
            // 'species' is a keyword applicable to the current file's 'species'
            // section
            addSpecies(thermo, phaseNode["species"], rootNode["species"]);
        } else if (phaseNode["species"].is<vector<AnyMap>>()) {
            // Each item in 'species' is a map with one item, where the key is
            // a section in this file or another YAML file, and the value is a
            // list of species names to read from that section
            for (const auto& speciesNode : phaseNode["species"].asVector<AnyMap>()) {
                const string& source = speciesNode.begin()->first;
                const auto& names = speciesNode.begin()->second;
                const auto& slash = boost::ifind_last(source, "/");
                if (slash) {
                    // source is a different input file
                    std::string fileName(source.begin(), slash.begin());
                    std::string node(slash.end(), source.end());
                    AnyMap species = AnyMap::fromYamlFile(fileName,
                        rootNode.getString("__file__", ""));
                    addSpecies(thermo, names, species[node]);
                } else if (rootNode.hasKey(source)) {
                    // source is in the current file
                    addSpecies(thermo, names, rootNode[source]);
                } else {
                    throw InputFileError("setupPhase", speciesNode,
                        "Could not find species section named '{}'", source);
                }
            }
        } else {
            throw InputFileError("setupPhase", phaseNode["species"],
                "Could not parse species declaration of type '{}'",
                phaseNode["species"].type_str());
        }
    } else if (rootNode.hasKey("species")) {
        // By default, add all species from the 'species' section
        addSpecies(thermo, AnyValue("all"), rootNode["species"]);
    }

    auto* vpssThermo = dynamic_cast<VPStandardStateTP*>(&thermo);
    if (vpssThermo) {
        for (size_t k = 0; k < thermo.nSpecies(); k++) {
            unique_ptr<PDSS> pdss;
            if (thermo.species(k)->input.hasKey("equation-of-state")) {
                // Use the first node which specifies a valid PDSS model
                auto& eos = thermo.species(k)->input["equation-of-state"];
                bool found = false;
                for (auto& node : eos.asVector<AnyMap>()) {
                    string model = node["model"].asString();
                    if (PDSSFactory::factory()->exists(model)) {
                        pdss.reset(newPDSS(model));
                        pdss->setParameters(node);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw InputFileError("setupPhase", eos,
                        "Could not find an equation-of-state specification "
                        "which defines a known PDSS model.");
                }
            } else {
                pdss.reset(newPDSS("ideal-gas"));
            }
            vpssThermo->installPDSS(k, std::move(pdss));
        }
    }

    thermo.setParameters(phaseNode, rootNode);
    thermo.initThermo();

    if (phaseNode.hasKey("state")) {
        auto node = phaseNode["state"].as<AnyMap>();
        thermo.setState(node);
    } else {
        thermo.setState_TP(298.15, OneAtm);
    }
}

}
