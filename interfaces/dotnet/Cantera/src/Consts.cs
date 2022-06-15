namespace Cantera;

public static class Consts
{
    //! Avogadro's Number \f$ N_{\mathrm{A}} \f$ [number/kmol]
    public const double Avogadro = 6.02214076e26;

    //! Boltzmann constant \f$ k \f$ [J/K]
    public const double Boltzmann = 1.380649e-23;

    //! Planck constant \f$ h \f$ [J-s]
    public const double Planck = 6.62607015e-34;

    //! Elementary charge \f$ e \f$ [C]
    public const double ElectronCharge = 1.602176634e-19;

    //! Speed of Light in a vacuum \f$ c \f$ [m/s]
    public const double LightSpeed = 299792458.0;

    //! Electron Mass \f$ m_e \f$ [kg]
    public const double ElectronMass = 9.1093837015e-31;

    //! Universal Gas Constant \f$ R_u \f$ [J/kmol/K]
    public const double GasConstant = Avogadro * Boltzmann;

    //! Faraday constant \f$ F \f$ [C/kmol]
    public const double Faraday = ElectronCharge * Avogadro;

    //! Stefan-Boltzmann constant \f$ \sigma \f$ [W/m2/K4]
    public const double StefanBoltz = 5.670374419e-8;

    //! One atmosphere [Pa]
    public const double OneAtm = 1.01325e5;

    //! One bar [Pa]
    public const double OneBar = 1.0E5;
}
