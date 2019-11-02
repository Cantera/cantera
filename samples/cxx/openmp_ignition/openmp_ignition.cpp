// Ignition delay calculation with OpenMP. This example shows how to use OpenMP
// to run multiple reactor network calculations in parallel by using separate
// Cantera objects for each thread.

#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"

#include <omp.h>

using namespace Cantera;

void run()
{
    // The number of threads can be set by setting the OMP_NUM_THREADS
    // environment variable before running the code.
    int nThreads = omp_get_max_threads();
    writelog("Running on {} threads\n\n", nThreads);

    // Containers for Cantera objects to be used in different. Each thread needs
    // to have its own set of linked Cantera objects. Multiple threads accessing
    // the same objects at the same time will cause errors.
    std::vector<std::shared_ptr<Solution>> sols;
    std::vector<std::unique_ptr<IdealGasConstPressureReactor>> reactors;
    std::vector<std::unique_ptr<ReactorNet>> nets;

    // Create and link the Cantera objects for each thread. This step should be
    // done in serial
    for (int i = 0; i < nThreads; i++) {
        auto sol = newSolution("gri30.yaml", "gri30", "None");
        sols.emplace_back(sol);
        reactors.emplace_back(new IdealGasConstPressureReactor());
        nets.emplace_back(new ReactorNet());
        reactors.back()->insert(sol);
        nets.back()->addReactor(*reactors.back());
    }

    // Points at which to compute ignition delay time
    int nPoints = 50;
    vector_fp T0(nPoints);
    vector_fp ignition_time(nPoints);
    for (int i = 0; i < nPoints; i++) {
        T0[i] = 1000 + 500 * ((float) i) / ((float) nPoints);
    }

    // Calculate the ignition delay at each initial temperature using multiple
    // threads.
    //
    // Note on 'schedule(static, 1)':
    // This option causes points [0, nThreads, 2*nThreads, ...] to be handled by
    // the same thread, rather than the default behavior of one thread handling
    // points [0 ... nPoints/nThreads]. This helps balance the workload for each
    // thread in cases where the workload is biased, e.g. calculations for low
    // T0 take longer than calculations for high T0.
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < nPoints; i++) {
        // Get the Cantera objects that were initialized for this thread
        size_t j = omp_get_thread_num();
        auto gas = sols[j]->thermo();
        Reactor& reactor = *reactors[j];
        ReactorNet& net = *nets[j];

        // Set up the problem
        gas->setState_TPX(T0[i], OneAtm, "CH4:0.5, O2:1.0, N2:3.76");
        reactor.syncState();
        net.setInitialTime(0.0);

        // Integrate until we satisfy a crude estimate of the ignition delay
        // time: time for T to increase by 500 K
        while (reactor.temperature() < T0[i] + 500) {
            net.step();
        }

        // Save the ignition delay time for this temperature
        ignition_time[i] = net.time();
    }

    // Print the computed ignition delays
    writelog("  T (K)    t_ig (s)\n");
    writelog("--------  ----------\n");
    for (int i = 0; i < nPoints; i++) {
        writelog("{: 8.1f}  {: 10.3e}\n", T0[i], ignition_time[i]);
    }
}

int main()
{
    try {
        run();
        appdelete();
        return 0;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        std::cout << " terminating... " << std::endl;
        appdelete();
        return 1;
    }
}
