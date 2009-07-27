/// MorseGas.h

namespace Cantera {

  class Zvib {
  public:
    Zvib(const MorseVib& m);
    Real operator() {
      return 
  class MorseVib {
  public:

    /// Constructor
    MorseVib(Real freq0, Real De);

    /// Destructor
    virtual ~MorseVib() {}

    Real operator()(Real T);

  private:
    
  };
  
