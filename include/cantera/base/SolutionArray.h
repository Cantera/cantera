//! @file SolutionArray.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SOLUTIONARRAY_H
#define CT_SOLUTIONARRAY_H

#include "cantera/base/global.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

class Solution;
class ThermoPhase;

//! A container class holding arrays of state information.
/*!
 * SolutionArray objects provide a convenient interface for representing many
 * thermodynamic states using the same Solution object. C++ SolutionArray objects are
 * one-dimensional by design; while shape information for multi-dimensional arrays is
 * stored, reshaping operations need to be implemented in high-level API's.
 *
 * The SolutionArray class implements the main interface for saving and restoring of
 * %Cantera simulation data. SolutionArray objects can be serialized to and from YAML and
 * HDF container files using the save() and restore() methods. In addition, there is
 * limited support for CSV files.
 * @since New in %Cantera 3.0.
 * @ingroup solnGroup
 */
class SolutionArray
{
private:
    SolutionArray(const shared_ptr<Solution>& sol,
                  int size,
                  const AnyMap& meta);

    SolutionArray(const SolutionArray& arr, const vector<int>& indices);

public:
    virtual ~SolutionArray();

    /**
     *  Instantiate a new SolutionArray reference.
     *
     *  @param sol  Solution object defining phase definitions
     *  @param size  Number of SolutionArray entries
     *  @param meta  AnyMap holding SolutionArray meta data
     */
    static shared_ptr<SolutionArray> create(const shared_ptr<Solution>& sol,
                                            int size=0,
                                            const AnyMap& meta={})
    {
        return shared_ptr<SolutionArray>(new SolutionArray(sol, size, meta));
    }

    /**
     *  Share locations from an existing SolutionArray and return new reference.
     *
     *  Both SolutionArray object share common data. The method is used for slicing
     *  of SolutionArrays from high-level API's. Note that meta data are not inherited.
     *  @param selected  List of locations for shared entries
     */
    shared_ptr<SolutionArray> share(const vector<int>& selected)
    {
        return shared_ptr<SolutionArray>(new SolutionArray(*this, selected));
    }

    //! Reset all entries of the SolutionArray to the current Solution state
    void reset();

    //! Size of SolutionArray (number of entries).
    int size() const {
        return static_cast<int>(m_size);
    }

    //! Resize SolutionArray objects with a single dimension (default).
    void resize(int size);

    //! SolutionArray shape information used by high-level API's.
    vector<long int> apiShape() const {
        return m_apiShape;
    }

    //! Set SolutionArray shape information used by high-level API's.
    //! The size of the SolutionArray is adjusted automatically.
    void setApiShape(const vector<long int>& shape);

    //! Number of SolutionArray dimensions used by high-level API's.
    int apiNdim() const {
        return static_cast<int>(m_apiShape.size());
    }

    /**
     *  Print a concise summary of a SolutionArray.
     *  @param keys  List of components to be displayed; if empty, all components are
     *      considered.
     *  @param rows  Maximum number of rendered rows.
     *  @param width  Maximum width of rendered output.
     */
    string info(const vector<string>& keys, int rows=10, int width=80);

    //! SolutionArray meta data.
    AnyMap& meta() {
        return m_meta;
    }

    //! Set SolutionArray meta data.
    void setMeta(const AnyMap& meta) {
        m_meta = meta;
    }

    //! Retrieve associated Solution object
    shared_ptr<Solution> solution() {
        return m_sol;
    }

    //! Retrieve associated ThermoPhase object
    shared_ptr<ThermoPhase> thermo();

    //! Retrieve list of component names
    vector<string> componentNames() const;

    /**
     *  Check whether SolutionArray contains a component.
     *  A component is a property defining state or auxiliary variable.
     */
    bool hasComponent(const string& name) const;

    /**
     *  Retrieve a component of the SolutionArray by name.
     *  Returns an AnyValue containing an array with length size() with a type
     *  specific to the component; in most cases, the type is double, but may differ
     *  for auxiliary data.
     */
    AnyValue getComponent(const string& name) const;

    /**
     *  Set a component of the SolutionArray by name.
     *  The passed AnyValue should containing an array with length size() with a type
     *  specific to the component; in most cases, the type is double, but may differ
     *  for auxiliary data.
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param data  Component data
     */
    void setComponent(const string& name, const AnyValue& data);

    /**
     *  Update the buffered location used to access SolutionArray entries.
     */
    void setLoc(int loc, bool restore=true);

    /**
     *  Update state at given location to state of associated Solution object.
     */
    void updateState(int loc);

    //! Retrieve the state vector for a given location.
    vector<double> getState(int loc);

    //! Set the state vector for a given location.
    void setState(int loc, const vector<double>& state);

    //! Normalize mass/mole fractions
    void normalize();

    /**
     *  Add auxiliary component to SolutionArray. Initialization requires a subsequent
     *  call of setComponent().
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param back  If `true` (default), add name after components representing the
     *      state, otherwise add to front of list. Front and back components are
     *      populated left to right.
     */
    void addExtra(const string& name, bool back=true);

    //! Check whether SolutionArray contains an extra component
    bool hasExtra(const string& name) const {
        return m_extra->count(name);
    }

    //! Retrieve list of extra component names
    vector<string> listExtra(bool all=true) const;

    //! Retrieve auxiliary data for a given location.
    AnyMap getAuxiliary(int loc);

    //! Set auxiliary data for a given location.
    void setAuxiliary(int loc, const AnyMap& data);

    //! Append location entry at end of SolutionArray.
    void append(const vector<double>& state, const AnyMap& extra);

    /**
     *  Write header data to a HDF container file.
     *
     *  @param fname  Name of HDF container file
     *  @param name  Identifier of group holding header information
     *  @param desc  Custom comment describing dataset
     *  @param overwrite  Force overwrite if file/group exists;
     *      optional (default=`false`)
     */
    static void writeHeader(const string& fname, const string& name, const string& desc,
                            bool overwrite=false);

    /**
     *  Write header data to AnyMap. Used by YAML serialization.
     *
     *  @param root  Root node of AnyMap structure
     *  @param name  Identifier of node holding header information
     *  @param desc  Custom comment describing dataset
     *  @param overwrite  Force overwrite if node exists; optional (default=`false`)
     */
    static void writeHeader(AnyMap& root, const string& name, const string& desc,
                            bool overwrite=false);

    /**
     *  Write SolutionArray data to a CSV file.
     *
     *  @param fname  Name of CSV file
     *  @param overwrite  Force overwrite if file exists; optional (default=`false`)
     *  @param basis  Output mass (`"Y"`/`"mass"`) or mole (`"X"`/`"mole"`) fractions;
     *      if omitted (default=`""`), the native basis of the underlying ThermoPhase
     *      manager is used - see Phase::nativeState
     */
    void writeEntry(const string& fname, bool overwrite=false, const string& basis="");

    /**
     *  Write SolutionArray data to a HDF container file.
     *
     *  @param fname  Name of HDF container file
     *  @param name  Identifier of group holding header information
     *  @param sub  Name identifier of subgroup holding SolutionArray data
     *  @param overwrite  Force overwrite if subgroup exists; optional (default=`false`)
     *  @param compression  Compression level; optional (default=0; HDF only)
     */
    void writeEntry(const string& fname, const string& name, const string& sub,
                    bool overwrite=false, int compression=0);

    /**
     *  Write SolutionArray data to AnyMap. Used by YAML serialization.
     *
     *  @param root  Root node of AnyMap structure
     *  @param name  Identifier of node holding header information and subgroup
     *  @param sub  Name identifier of subgroup holding SolutionArray data
     *  @param overwrite  Force overwrite if subgroup exists; optional (default=`false`)
     */
    void writeEntry(AnyMap& root, const string& name, const string& sub,
                    bool overwrite=false);

    /**
     *  Save current SolutionArray contents to a data file.
     *
     *  Data can be saved either in CSV format (extension `*.csv`), YAML container
     *  format (extension `*.yaml`/`*.yml`) or HDF container format (extension
     *  `*.h5`/`*.hdf5`/`*.hdf`). The output format is automatically inferred from the
     *  file extension.
     *
     *  CSV files preserve state data and auxiliary data for a single SolutionArray in a
     *  comma-separated text format, container files may hold multiple SolutionArray
     *  entries in an internal hierarchical structure. While YAML is a human-readable
     *  text format, HDF is a binary format that supports compression and is recommended
     *  for large datasets.
     *
     *  For container files (YAML and HDF), header information contains automatically
     *  generated time stamps, version information and an optional description.
     *  Container files also preserve SolutionArray metadata (example: SolutionArray
     *  objects generated by Sim1D hold simulation settings).
     *
     *  @param fname  Name of output file (CSV, YAML or HDF)
     *  @param name  Identifier of location within the container file; this node/group
     *      contains header information and a subgroup holding actual SolutionArray data
     *      (YAML/HDF only)
     *  @param sub  Name identifier for the subgroup holding the SolutionArray data and
     *      metadata objects. If omitted (`""`), the subgroup name defaults to `"data"`
     *      (YAML/HDF only)
     *  @param desc  Custom comment describing dataset to be stored (YAML/HDF only)
     *  @param overwrite  Force overwrite if file and/or data entry exists; optional
     *      (default=`false`)
     *  @param compression  Compression level (0-9); (default=0; HDF only)
     *  @param basis  Output mass (`"Y"`/`"mass"`) or mole (`"X"`/`"mole"`) fractions;
     *      if not specified (default=`""`), the native basis of the underlying
     *      ThermoPhase manager is used - see Phase::nativeState (CSV only)
     */
    void save(const string& fname, const string& name="", const string& sub="",
              const string& desc="", bool overwrite=false, int compression=0,
              const string& basis="");

    /**
     *  Read header information from a HDF container file.
     *
     *  @param fname  Name of HDF container file
     *  @param name  Identifier of group holding header information
     */
    static AnyMap readHeader(const string& fname, const string& name);

    /**
     *  Read header information from AnyMap. Used by YAML serialization.
     *
     *  @param root  Root node of AnyMap structure
     *  @param name  Identifier of node holding header information
     */
    static AnyMap readHeader(const AnyMap& root, const string& name);

    /**
     *  Restore SolutionArray data from a HDF container file.
     *
     *  @param fname  Name of HDF container file
     *  @param name  Identifier of group holding header information
     *  @param sub  Name identifier of subgroup holding SolutionArray data
     */
    void readEntry(const string& fname, const string& name, const string& sub);

    /**
     *  Restore SolutionArray data from AnyMap. Used by YAML serialization.
     *
     *  @param root  Root node of AnyMap structure
     *  @param name  Identifier of node holding header information
     *  @param sub  Name identifier of subgroup holding SolutionArray data
     */
    void readEntry(const AnyMap& root, const string& name, const string& sub);

    /**
     *  Restore SolutionArray data and header information from a container file.
     *
     *  This method retrieves data from a YAML or HDF files that were previously saved
     *  using the save() method.
     *
     *  @param fname  Name of container file (YAML or HDF)
     *  @param name  Identifier of location within the container file; this node/group
     *      contains header information and a subgroup holding actual SolutionArray data
     *  @param sub  Name identifier for the subgroup holding the SolutionArray data and
     *      metadata objects. If omitted (`""`), the subgroup name defaults to "data"
     *  @return  AnyMap containing header information
     */
    AnyMap restore(const string& fname, const string& name, const string& sub="");

protected:
    //! Service function used to resize SolutionArray
    void _resize(size_t size);

    /**
     *  Initialize extra SolutionArray component.
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param value  Default value; used to determine type of component
     */
    void _initExtra(const string& name, const AnyValue& value);

    /**
     *  Resize extra SolutionArray component.
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param value  Default value
     */
    void _resizeExtra(const string& name, const AnyValue& value=AnyValue());

    /**
     *  Set extra SolutionArray component
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param data  Value to be set
     */
    void _setExtra(const string& name, const AnyValue& data=AnyValue());

    /**
     *  Identify storage mode of state data. The storage mode is a combination of
     *  properties defining state); valid modes include Phase::nativeState (`"native"`)
     *  or other property combinations defined by Phase::fullStates (three-letter
     *  acronyms, for example `"TDY"`, `"TPX"`).
     */
    string _detectMode(const set<string>& names, bool native=true);

    //! Retrieve set containing list of properties defining state
    set<string> _stateProperties(const string& mode, bool alias=false);

    shared_ptr<Solution> m_sol; //!< Solution object associated with state data
    size_t m_size; //!< Number of entries in SolutionArray
    size_t m_dataSize; //!< Total size of unsliced data
    size_t m_stride; //!< Stride between SolutionArray entries
    AnyMap m_meta; //!< Metadata
    size_t m_loc = npos; //!< Buffered location within data vector
    vector<long int> m_apiShape; //!< Shape information used by high-level API's

    shared_ptr<vector<double>> m_data; //!< Work vector holding states

    //! Auxiliary (extra) components; size of first dimension has to match m_dataSize
    shared_ptr<map<string, AnyValue>> m_extra;

    //! Mapping of auxiliary component names, where the index is used as the
    //! mapping key. Names with index >= zero are listed before state components, while
    //! names with index < zero are added at end. The name with the most negative index
    //! corresponds to the last entry (different from Python index convention).
    shared_ptr<map<int, string>> m_order;

    bool m_shared = false; //!< `true` if data are shared from another object
    vector<int> m_active; //!< Vector of locations referencing active entries
};

}

#endif
