/**
 *  @file ElectronCrossSection.cpp Definition file for class ElectronCrossSection.
 */
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ElectronCrossSection.h"
#include "cantera/base/global.h"

namespace Cantera {

ElectronCrossSection::ElectronCrossSection()
    : threshold(0.0)
{
}

ElectronCrossSection::~ElectronCrossSection()
{
}

// void ElectronCrossSection::validate()
// {
//     if (kind == "ionization" || kind == "attachment" || kind == "excitation") {
//         if (threshold < 0.0) {
//             throw CanteraError("ElectronCrossSection::validate",
//                                "The threshold of the process",
//                                "(kind = '{}', target = '{}', product = '{}')",
//                                "cannot be negative", kind, target, product);
//         }
//     } else if (kind != "effective" && kind != "elastic") {
//         throw CanteraError("ElectronCrossSection::validate",
//             "'{}' is an unknown type of cross section data.", kind);
//     }
// }

unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node)
{
    unique_ptr<ElectronCrossSection> ecs(new ElectronCrossSection());

    ecs->kind = node["kind"].asString();
    ecs->target = node["target"].asString();

    auto& data = node["data"].asVector<vector<double>>();
    for (size_t i = 0; i < data.size(); i++) {
        ecs->energyLevel.push_back(data[i][0]);
        ecs->crossSection.push_back(data[i][1]);
    }

    if (node.hasKey("threshold")){
        ecs->threshold = node["threshold"].asDouble(); //std::stof(node.attrib("threshold"));
    } else {
        ecs->threshold = 0.0;

        if (ecs->kind == "excitation" || ecs->kind == "ionization" || ecs->kind == "attachment") {
            for (size_t i = 0; i < ecs->energyLevel.size(); i++) {
                if (ecs->energyLevel[i] > 0.0) {  // Look for first non-zero cross-section
                    ecs->threshold = ecs->energyLevel[i];
                    break;
                }
            }
        }

    }

    if (node.hasKey("products")) {
        ecs->products = node["products"].asVector<std::string>();  // Store all products
        ecs->product = ecs->products.empty() ? ecs->target : ecs->products[0];  // Keep first product for compatibility
    } else if (node.hasKey("product")) {
        ecs->product = node["product"].asString();
        ecs->products.push_back(ecs->product);  // Ensure products list is always populated
    } else {
        ecs->product = ecs->target;
        ecs->products.push_back(ecs->product);
    }

    /*if (node.hasKey("product")) {
        ecs->product = node["product"].asString();
    } else {
        ecs->product = ecs->target;
    }*/

    // Some writelog to check the datas loaded concerning the cross section
    //writelog("Kind :    {:s}\n",ecs->kind);
    //writelog("Target :    {:s}\n",ecs->target);
    //writelog("Product :    {:s}\n",ecs->product);
    //writelog("Threshold :    {:14.5g} eV\n",ecs->threshold);
    //writelog("Energy :  \n");
    //for (size_t i = 0; i < ecs->energyLevel.size(); i++){
    //   writelog("{:9.4g} {:9.4g} \n",ecs->energyLevel[i], ecs->crossSection[i]);
    //}

    return ecs;
}

}