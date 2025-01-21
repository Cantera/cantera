#include "gtest/gtest.h"
#include "cantera/base/Interface.h"
#include "cantera/base/SolutionArray.h"

using namespace Cantera;

TEST(Interface, incompatible_phase)
{
    ASSERT_THROW(newInterface("h2o2.yaml", "ohmech"), CanteraError);
    AnyMap h2o2 = AnyMap::fromYamlFile("h2o2.yaml");
    auto& ohmech = h2o2["phases"].getMapWhere("name", "ohmech");
    ASSERT_THROW(newInterface(ohmech, h2o2), CanteraError);
}

TEST(Interface, from_string_and_Solution)
{
    auto gas = newSolution("ptcombust.yaml", "gas");
    auto surf = newInterface("ptcombust.yaml", "Pt_surf", {gas});
    ASSERT_EQ(gas.get(), surf->adjacent(0).get());
}

TEST(Interface, from_AnyMap_and_Solution)
{
    auto gas = newSolution("ptcombust.yaml", "gas");
    AnyMap root = AnyMap::fromYamlFile("ptcombust.yaml");
    auto& phase = root["phases"].getMapWhere("name", "Pt_surf");
    auto surf = newInterface(phase, root, {gas});
    ASSERT_EQ(gas.get(), surf->adjacent(0).get());
    ASSERT_EQ(surf->kinetics()->nReactions(), 24u);
}

TEST(SolutionArray, empty)
{
    shared_ptr<Solution> gas;
    // gas must be initialized
    ASSERT_THROW(SolutionArray::create(gas, 0, AnyMap()), CanteraError);

    gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas);

    ASSERT_EQ(arr->size(), 0);
    ASSERT_EQ(arr->meta(), AnyMap());

    AnyValue vec;
    vec = arr->getComponent("T");
    ASSERT_TRUE(vec.isVector<double>());
    ASSERT_EQ(vec.vectorSize(), 0u);

    vec = arr->getComponent("H2");
    ASSERT_TRUE(vec.isVector<double>());
    ASSERT_EQ(vec.vectorSize(), 0u);

    ASSERT_THROW(arr->getAuxiliary(0), CanteraError);
    ASSERT_THROW(arr->getComponent("foo"), CanteraError);
}

TEST(SolutionArray, simple)
{
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 5);

    ASSERT_EQ(gas->transportModel(), "none");
    ASSERT_EQ(arr->transportModel(), "none");

    ASSERT_EQ(arr->size(), 5);
    ASSERT_EQ(arr->meta(), AnyMap());

    AnyValue vec;
    vec = arr->getComponent("T");
    ASSERT_TRUE(vec.isVector<double>());
    ASSERT_EQ(vec.vectorSize(), 5u);

    vec = arr->getComponent("H2");
    ASSERT_TRUE(vec.isVector<double>());
    ASSERT_EQ(vec.vectorSize(), 5u);

    ASSERT_THROW(arr->getComponent("foo"), CanteraError);

    ASSERT_EQ(arr->getAuxiliary(0), AnyMap());
    AnyMap m;
    m["spam"] = "eggs";
    ASSERT_THROW(arr->setAuxiliary(0, m), CanteraError);
}

TEST(SolutionArray, normalize)
{
    auto gas = newSolution("h2o2.yaml",  "", "none");

    auto arr = SolutionArray::create(gas, 10, AnyMap());
    vector<double> state = arr->getState(0); // size is 12
    vector<double> newState = {state[0], state[1], 1, 1, 1, 1, 1, 1, 1, 1, 1, -1e12};
    ASSERT_EQ(newState.size(), state.size());
    for (int loc = 0; loc < arr->size(); loc++) {
        arr->setState(loc, newState);
        state = arr->getState(loc);
        for (size_t i = 0; i < state.size(); i++) {
            ASSERT_EQ(state[i], newState[i]);
        }
    }
    arr->normalize();
    for (int loc = 0; loc < arr->size(); loc++) {
        state = arr->getState(loc);
        for (size_t i = 2; i < state.size() - 1; i++) {
            // test normalized species mass fractions - all are equal except last
            ASSERT_NEAR(state[i], state[2], 1e-12);
            ASSERT_NE(state[i], newState[i]);
            ASSERT_LE(state[i], 1.);
            ASSERT_GE(state[i], 0.);
        }
        ASSERT_EQ(state.back(), 0.); // unnormalized negative value is set to zero
    }

    newState = {state[0], state[1], 100, 0, 0, 0, 0, 0, 0, 0, 0, -1e12};
    ASSERT_EQ(newState.size(), state.size());
    for (int loc = 0; loc < arr->size(); loc++) {
        arr->setState(loc, newState);
        state = arr->getState(loc);
        for (size_t i = 0; i < state.size(); i++) {
            ASSERT_EQ(state[i], newState[i]);
        }
    }
    arr->normalize();
    for (int loc = 0; loc < arr->size(); loc++) {
        state = arr->getState(loc);
        ASSERT_EQ(state[2], 1.);
        for (size_t i = 3; i < state.size(); i++) {
            ASSERT_EQ(state[i], 0.);
        }
    }
}

TEST(SolutionArray, meta)
{
    auto gas = newSolution("h2o2.yaml",  "", "none");

    AnyMap meta;
    meta["spam"] = "eggs";
    meta["foo"]["bar"] = 42;
    meta["flag"] = true;
    auto arr = SolutionArray::create(gas, 0, meta);
    ASSERT_EQ(arr->size(), 0);
    ASSERT_EQ(arr->meta(), meta);
}

TEST(SolutionArray, extraEmpty) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas);
    ASSERT_EQ(arr->size(), 0);

    arr->addExtra("spam");
    AnyValue any;
    any = arr->getComponent("spam");
    ASSERT_TRUE(any.is<void>());

    any = 1;
    ASSERT_THROW(arr->setComponent("spam", any), CanteraError);
    any = vector<double>({1});
    ASSERT_THROW(arr->setComponent("spam", any), CanteraError);
    any = vector<double>();
    arr->setComponent("spam", any);
    any = arr->getComponent("spam");
    ASSERT_TRUE(any.isVector<double>());
    ASSERT_EQ(any.asVector<double>().size(), 0u);
    arr->setComponent("spam", AnyValue());
    any = arr->getComponent("spam");
    ASSERT_TRUE(any.is<void>());

    // component must exist
    ASSERT_THROW(arr->setComponent("eggs", any), CanteraError);
}

TEST(SolutionArray, extraSingle) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 1);
    ASSERT_EQ(arr->size(), 1);

    arr->addExtra("spam");
    AnyValue any;

    // set entire component
    any = vector<double>({3.1415});
    ASSERT_TRUE(any.isVector<double>());
    ASSERT_EQ(any.asVector<double>().size(), 1u);

    // reset
    arr->setComponent("spam", AnyValue());
    any = arr->getComponent("spam");
    ASSERT_TRUE(any.is<void>());

    // initialize using default values: scalar
    any = 1;
    arr->setComponent("spam", any);
    any = arr->getComponent("spam");
    ASSERT_TRUE(any.isVector<long int>());
    ASSERT_EQ(any.asVector<long int>().size(), 1u);
    any = 2;
    arr->setComponent("spam", any);
    any = vector<double>({2.});
    arr->setComponent("spam", any);

    // initialize using default values: vector with matching size
    arr->setComponent("spam", AnyValue());
    any = vector<string>({"eggs"});
    arr->setComponent("spam", any);
    any = arr->getComponent("spam");
    ASSERT_TRUE(any.isVector<string>());
    ASSERT_EQ(any.asVector<string>().size(), 1u);

    // replace/initialize using default values: vectors
    any = vector<long int>({3, 4, 5});
    ASSERT_THROW(arr->setComponent("spam", any), CanteraError);
    arr->setComponent("spam", AnyValue()); // needs reset first
    arr->setComponent("spam", any);
    any = arr->getComponent("spam");
    ASSERT_TRUE(any.isMatrix<long int>());
    ASSERT_EQ(any.matrixShape().first, 1u);
    ASSERT_EQ(any.matrixShape().second, 3u);

    // attempt to reset with empty vector
    arr->setComponent("spam", AnyValue());
    any = vector<double>();
    ASSERT_THROW(arr->setComponent("spam", any), CanteraError);
}

template<class T>
void testSingleCol(SolutionArray& arr, const T& value, bool sliced=false)
{
    int size = arr.size();
    T defaultValue = vector<T>(1)[0];

    AnyValue any;
    any = value;
    ASSERT_TRUE(any.is<T>());

    if (sliced) {
        AnyValue slice;
        slice = vector<T>(size, value);
        arr.setComponent("spam", slice);
    } else {
        // initialize
        arr.setComponent("spam", any);
    }
    AnyValue extra;
    extra = arr.getComponent("spam");
    ASSERT_TRUE(extra.isVector<T>());
    ASSERT_EQ(len(extra.asVector<T>()), arr.size());
    for (int i = 0; i < size; ++i) {
        ASSERT_EQ(extra.asVector<T>()[i], value);
    }

    // set directly; extra contains values with correct dimensions
    arr.setComponent("spam", extra);

    // remove one entry and test failure
    extra.as<vector<T>>().pop_back();
    ASSERT_THROW(arr.setComponent("spam", extra), CanteraError);

    if (sliced) {
        ASSERT_THROW(arr.resize(2 * size), CanteraError);
    } else {
        arr.resize(2 * size);
        extra = arr.getComponent("spam");
        ASSERT_TRUE(extra.isVector<T>());
        ASSERT_EQ(arr.size(), 2 * size);
        for (int i = 0; i < size; ++i) {
            ASSERT_EQ(extra.asVector<T>()[i], value);
            ASSERT_EQ(extra.asVector<T>()[i + size], defaultValue);
        }
    }

    arr.reset();
    extra = arr.getComponent("spam");
    AnyMap m;
    m["spam"] = value;
    for (int i = 0; i < arr.size(); ++i) {
        ASSERT_EQ(extra.asVector<T>()[i], defaultValue);
        arr.setAuxiliary(i, m);
        ASSERT_EQ(arr.getAuxiliary(i)["spam"].as<T>(), value);
    }

    // replace all data with a scalar
    any = value;
    arr.setComponent("spam", any);
    extra = arr.getComponent("spam");
    for (int i = 0; i < arr.size(); ++i) {
        ASSERT_EQ(extra.asVector<T>()[i], value);
    }
}

TEST(SolutionArray, extraHoldsDoubles) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 3);
    arr->addExtra("spam");

    testSingleCol<double>(*arr, 3.1415);
}

TEST(SolutionArray, extraHoldsIntegers) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 5);
    arr->addExtra("spam");

    testSingleCol<long int>(*arr, 42);
    ASSERT_EQ(arr->getAuxiliary(1)["spam"].asInt(), 42);
}

TEST(SolutionArray, extraHoldsStrings) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 7);
    arr->addExtra("spam");

    testSingleCol<string>(*arr, "foo");
}

TEST(SolutionArray, extraSlicedDoubles) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 5);
    AnyValue any;
    any = 3.1415;
    arr->addExtra("spam");
    arr->setComponent("spam", any);

    auto indices = vector<int>({1, 2, 3});
    auto sliced = arr->share(indices);
    testSingleCol<double>(*sliced, 2.71828, true);
}

TEST(SolutionArray, extraSlicedIntegers) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 7);
    AnyValue any;
    any = 42;
    arr->addExtra("spam");
    arr->setComponent("spam", any);

    AnyMap m;
    m["spam"] = 3;
    ASSERT_EQ(arr->getAuxiliary(1)["spam"].asInt(), 42);
    arr->setAuxiliary(1, m);
    ASSERT_EQ(arr->getAuxiliary(1)["spam"].asInt(), 3);

    auto indices = vector<int>({1, 2, 3});
    auto sliced = arr->share(indices);
    ASSERT_EQ(sliced->getAuxiliary(0)["spam"].asInt(), 3);
    ASSERT_EQ(sliced->getAuxiliary(1)["spam"].asInt(), 42);
    testSingleCol<long int>(*sliced, 101, true);
    ASSERT_EQ(sliced->getAuxiliary(0)["spam"].asInt(), 101);
    ASSERT_EQ(arr->getAuxiliary(1)["spam"].asInt(), 101);
}

TEST(SolutionArray, extraSlicedStrings) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 6);
    AnyValue any;
    any = "foo";
    arr->addExtra("spam");
    arr->setComponent("spam", any);

    auto indices = vector<int>({1, 2, 3});
    auto sliced = arr->share(indices);
    testSingleCol<string>(*sliced, "bar", true);
}

template<class T>
void testMultiCol(SolutionArray& arr, const vector<T>& value, bool sliced=false)
{
    int size = arr.size();
    T defaultValue = vector<T>(1)[0];

    AnyValue any;
    any = value;
    ASSERT_TRUE(any.isVector<T>());

    // initialize
    if (sliced) {
        AnyValue slice;
        slice = vector<vector<T>>(size, value);
        arr.setComponent("spam", slice);
    } else {
        // initialize
        arr.setComponent("spam", any);
    }
    AnyValue extra;
    extra = arr.getComponent("spam");
    ASSERT_TRUE(extra.isMatrix<T>());
    ASSERT_EQ(extra.matrixShape().first, static_cast<size_t>(arr.size()));
    ASSERT_EQ(extra.matrixShape().second, any.vectorSize());
    for (int i = 0; i < size; ++i) {
        for (size_t j = 0; j < value.size(); ++j) {
            ASSERT_EQ(extra.asVector<vector<T>>()[i][j], value[j]);
        }
    }

    // set directly; extra contains values with correct dimensions
    arr.setComponent("spam", extra);

    // remove one entry and test failure
    extra.asVector<vector<T>>().pop_back();
    ASSERT_THROW(arr.setComponent("spam", extra), CanteraError);

    // remove one nested entry (make matrix invalid) and test failure
    extra = arr.getComponent("spam");
    extra.asVector<vector<T>>()[0].pop_back();
    ASSERT_THROW(arr.setComponent("spam", extra), CanteraError);

    if (sliced) {
        ASSERT_THROW(arr.resize(2 * size), CanteraError);
    } else {
        arr.resize(2 * size);
        extra = arr.getComponent("spam");
        ASSERT_EQ(arr.size(), 2 * size);
        ASSERT_EQ(len(extra.asVector<vector<T>>()), 2 * size);
        for (int i = 0; i < size; ++i) {
            for (size_t j = 0; j < value.size(); ++j) {
                ASSERT_EQ(extra.asVector<vector<T>>()[i][j], value[j]);
                ASSERT_EQ(extra.asVector<vector<T>>()[i + size][j], defaultValue);
            }
        }
    }

    AnyMap m;
    m["spam"] = value;
    arr.reset();
    extra = arr.getComponent("spam");
    ASSERT_TRUE(extra.isMatrix<T>());
    ASSERT_EQ(extra.matrixShape().second, value.size());
    for (int i = 0; i < arr.size(); ++i) {
        for (size_t j = 0; j < value.size(); ++j) {
            ASSERT_EQ(extra.asVector<vector<T>>()[i][j], defaultValue);
        }
        arr.setAuxiliary(i, m);
        for (size_t j = 0; j < value.size(); ++j) {
            ASSERT_EQ(arr.getAuxiliary(i)["spam"].asVector<T>()[j], value[j]);
        }
    }
}

TEST(SolutionArray, extraHoldsDoubleArrays) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 5);
    arr->addExtra("spam");

    testMultiCol<double>(*arr, vector<double>({1.1, 2.2}));
}

TEST(SolutionArray, extraHoldsIntegerArrays) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 6);
    arr->addExtra("spam");

    testMultiCol<long int>(*arr, vector<long int>({3, 4, 5}));
}

TEST(SolutionArray, extraHoldsStringArrays) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 9);
    arr->addExtra("spam");

    testMultiCol<string>(*arr, vector<string>({"foo", "bar", "hello world!"}));
}

TEST(SolutionArray, extraSlicedDoubleArrays) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 5);
    AnyValue any;
    any = vector<double>({0., 10.});
    arr->addExtra("spam");
    arr->setComponent("spam", any);

    auto indices = vector<int>({1, 2, 3});
    auto sliced = arr->share(indices);

    testMultiCol<double>(*sliced, vector<double>({1.1, 2.2}), true);
}

TEST(SolutionArray, extraSlicedIntegerArrays) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 6);
    AnyValue any;
    any = vector<long int>({0, 1, 2});
    arr->addExtra("spam");
    arr->setComponent("spam", any);

    auto indices = vector<int>({1, 2, 3});
    auto sliced = arr->share(indices);

    testMultiCol<long int>(*sliced, vector<long int>({3, 4, 5}), true);
}

TEST(SolutionArray, extraSlicedStringArrays) {
    auto gas = newSolution("h2o2.yaml",  "", "none");
    auto arr = SolutionArray::create(gas, 9);
    AnyValue any;
    arr->addExtra("spam");
    any = vector<string>({"a", "b", "hello world!"});
    arr->setComponent("spam", any);

    AnyMap m;
    ASSERT_EQ(arr->getAuxiliary(1)["spam"].asVector<string>()[0], "a");
    m["spam"] = vector<string>({"abc", "def"});
    ASSERT_THROW(arr->setAuxiliary(1, m), CanteraError);
    m["spam"] = vector<string>({"abc", "def", "ghi"});
    arr->setAuxiliary(1, m);
    ASSERT_EQ(arr->getAuxiliary(1)["spam"].asVector<string>()[0], "abc");

    auto indices = vector<int>({1, 2, 3});
    auto sliced = arr->share(indices);
    ASSERT_EQ(sliced->getAuxiliary(0)["spam"].asVector<string>()[0], "abc");
    ASSERT_EQ(sliced->getAuxiliary(1)["spam"].asVector<string>()[0], "a");
    testMultiCol<string>(*sliced, vector<string>({"foo", "bar", "baz"}), true);
}
