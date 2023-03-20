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

/*!
 *  A container class providing a convenient interface for representing many
 *  thermodynamic states using the same Solution object. C++ SolutionArray objects are
 *  one-dimensional by design; while shape information for multi-dimensional arrays is
 *  stored, reshaping operations need to be implemented in high-level API's.
 *
 *  @since  New in Cantera 3.0.
 *  @warning This class is an experimental part of the %Cantera API and may be
 *      changed or removed without notice.
 */
class SolutionArray
{
private:
    SolutionArray(const shared_ptr<Solution>& sol,
                  size_t size,
                  const AnyMap& meta);

    SolutionArray(const SolutionArray& arr, const vector<int>& indices);

public:
    virtual ~SolutionArray() {}

    /*!
     *  Instantiate a new SolutionArray reference
     *
     *  @param sol  Solution object defining phase definitions
     *  @param size  Number of SolutionArray entries
     *  @param meta  AnyMap holding SolutionArray meta data
     */
    static shared_ptr<SolutionArray> create(const shared_ptr<Solution>& sol,
                                            size_t size=0,
                                            const AnyMap& meta={})
    {
        return shared_ptr<SolutionArray>(new SolutionArray(sol, size, meta));
    }

    /*!
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

    //! Size of SolutionArray (number of entries)
    int size() const {
        return m_size;
    }

    //! Resize SolutionArray objects with a single dimension (default).
    void resize(size_t size);

    //! SolutionArray shape information used by high-level API's.
    vector<long int> apiShape() const {
        return m_apiShape;
    }

    //! Set SolutionArray shape information used by high-level API's.
    //! The size of the SolutionArray is adjusted automatically.
    void setApiShape(const vector<long int>& shape);

    //! Number of SolutionArray dimensions used by high-level API's.
    size_t apiNdim() const {
        return m_apiShape.size();
    }

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

    //! Check whether SolutionArray contains a component (property defining state or
    //! auxiliary variable)
    bool hasComponent(const string& name) const;

    /*!
     *  Retrieve a component of the SolutionArray by name.
     *  Returns an AnyValue containing an array with length size() with a type
     *  specific to the component; in most cases, the type is double, but may differ
     *  for auxiliary data.
     */
    AnyValue getComponent(const string& name) const;

    /*!
     *  Set a component of the SolutionArray by name.
     *  The passed AnyValue should containing an array with length size() with a type
     *  specific to the component; in most cases, the type is double, but may differ
     *  for auxiliary data.
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param data  Component data
     */
    void setComponent(const string& name, const AnyValue& data);

    /*!
     *  Update the buffered location used to access SolutionArray entries.
     */
    void setLoc(size_t loc, bool restore=true);

    /*!
     *  Update state at given location to state of associated Solution object.
     */
    void updateState(size_t loc);

    //! Retrieve the state vector for a given location.
    vector<double> getState(size_t loc);

    //! Set the state vector for a given location.
    void setState(size_t loc, const vector<double>& state);

    /*!
     *  Add auxiliary component to SolutionArray. Initialization requires a subsequent
     *  call of setComponent.
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param back  If true (default), add name after components representing the
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
    AnyMap getAuxiliary(size_t loc);

    //! Set auxiliary data for a given location.
    void setAuxiliary(size_t loc, const AnyMap& data);

    //! Append location entry at end of SolutionArray.
    void append(const vector<double>& state, const AnyMap& extra);

    /*!
     *  Write header data to container file.
     *
     *  @param fname  Name of HDF container file
     *  @param id  Identifier of root location within the container file
     *  @param desc  Description
     *  @param overwrite  Force overwrite if id exists; optional (default=false)
     */
    static void writeHeader(const string& fname, const string& id, const string& desc,
                            bool overwrite=false);

    /*!
     *  Write header data to AnyMap.
     *
     *  @param root  Root node of AnyMap structure
     *  @param id  Identifier of root location within the container file
     *  @param desc  Description
     *  @param overwrite  Force overwrite if id exists; optional (default=false)
     */
    static void writeHeader(AnyMap& root, const string& id, const string& desc,
                            bool overwrite=false);

    /*!
     *  Write SolutionArray data to container file.
     *
     *  @param fname  Name of HDF container file
     *  @param id  Identifier of root location within the container file
     *  @param sub  Name identifier for the subgroup holding actual data
     *  @param overwrite  Force overwrite if sub exists; optional (default=false)
     *  @param compression  Compression level; optional (default=0; HDF only)
     */
    void writeEntry(const string& fname, const string& id, const string& sub,
                    bool overwrite=false, int compression=0);

    /*!
     *  Write SolutionArray data to AnyMap.
     *
     *  @param root  Root node of AnyMap structure
     *  @param id  Identifier of root location within the container file
     *  @param sub  Name identifier for the subgroup holding actual data
     *  @param overwrite  Force overwrite if sub exists; optional (default=false)
     */
    void writeEntry(AnyMap& root, const string& id, const string& sub,
                    bool overwrite=false);

    /*!
     *  Save current SolutionArray and header to a container file.
     *
     *  @param fname  Name of output container file (YAML or HDF)
     *  @param id  Identifier of root location within the container file
     *  @param sub  Name identifier for the subgroup holding actual data
     *  @param desc  Custom comment describing the dataset to be stored
     *  @param overwrite  Force overwrite if sub exists; optional (default=false)
     *  @param compression  Compression level; optional (default=0; HDF only)
     */
    void save(const string& fname, const string& id, const string& sub,
              const string& desc, bool overwrite=false, int compression=0);

    /*!
     *  Read header data from container file.
     *
     *  @param fname  Name of HDF container file
     *  @param id  Identifier of root location within the container file
     */
    static AnyMap readHeader(const string& fname, const string& id);

    /*!
     *  Read header data from AnyMap.
     *
     *  @param root  Root node of AnyMap structure
     *  @param id  Identifier of root location within the container file
     */
    static AnyMap readHeader(const AnyMap& root, const string& id);

    /*!
     *  Restore SolutionArray entry from a container file.
     *
     *  @param fname  Name of HDF container file
     *  @param id  Identifier of root location within the container file
     *  @param sub  Name of the subgroup holding actual data
     */
    void readEntry(const string& fname, const string& id, const string& sub);

    /*!
     *  Restore SolutionArray entry from AnyMap.
     *
     *  @param root  Root node of AnyMap structure
     *  @param id  Identifier of root location within the container file
     *  @param sub  Name of the subgroup holding actual data
     */
    void readEntry(const AnyMap& root, const string& id, const string& sub);

    /*!
     *  Restore SolutionArray entry and header from a container file.
     *
     *  @param fname  Name of container file (YAML or HDF)
     *  @param id  Identifier of SolutionArray within the container file
     */
    AnyMap restore(const string& fname, const string& id, const string& sub);

protected:
    //! Service function used to resize SolutionArray
    void _resize(size_t size);

    /*!
     *  Initialize extra SolutionArray component
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param value  Default value; used to determine type of component
     */
    void _initExtra(const string& name, const AnyValue& value);

    /*!
     *  Resize extra SolutionArray component
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param value  Default value
     */
    void _resizeExtra(const string& name, const AnyValue& value=AnyValue());

    /*!
     *  Set extra SolutionArray component
     *
     *  @param name  Name of component (property defining auxiliary variable)
     *  @param data  Value to be set
     */
    void _setExtra(const string& name, const AnyValue& data=AnyValue());

    /*!
     *  Identify storage mode of state data (combination of properties defining state);
     *  valid modes include Phase::nativeState ("native") or other property combinations
     *  defined by Phase::fullStates (three-letter acronyms, for example "TDY", "TPX").
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

    bool m_shared = false; //!< True if data are shared from another object
    vector<int> m_active; //!< Vector of locations referencing active entries
};

}

#endif
