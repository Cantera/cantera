#include "gtest/gtest.h"
#include "cantera/numerics/polyfit.h"

using namespace Cantera;

double polyval(vector_fp& coeffs, double x) {
    double sum = 0;
    double xn = 1;
    for (size_t i = 0; i < coeffs.size(); i++) {
        sum += coeffs[i] * xn;
        xn *= x;
    }
    return sum;
}

TEST(Polyfit, exact_fit)
{
    vector_fp x{0, 0.3, 1.0, 1.5, 2.0, 2.5};
    vector_fp p(6);
    vector_fp w(6, -1.0);
    for (int i = 0; i < 20; i++) {
        vector_fp y{-1.1*i, cos(i), pow(-1,i), 3.2/(i+1), 0.1*i*i, sin(i)};
        polyfit(6, 5, x.data(), y.data(), w.data(), p.data());
        for (size_t j = 0; j < 6; j++) {
            EXPECT_NEAR(polyval(p, x[j]), y[j], 1e-10);
        }
    }
}

TEST(Polyfit, sequential)
{
    vector_fp x{-1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0};
    vector_fp y{0.6, 1.0, 0.8, 0.4, -0.1, -0.5, -1};

    // Coefficients calculated using Numpy's polyfit function for polynomials
    // of degrees 0 - 5.
    std::vector<vector_fp> PP{
        {0.17142857142857154},
        {0.66190476190476177, -0.49047619047619029},
        {0.73605442176870761, -0.19387755102040838, -0.14829931972789107},
        {1.0095838335334129, -0.22426970788315401, -0.51300520208083311,
         0.12156862745098072},
        {1.0121336003688943, -0.23102395749454527, -0.51552488317194212,
         0.12746543334778632, -0.0014742014742014889},
        {0.99812799812799835, -0.093488943488944404, -0.61193011193011071,
         0.011452361452361514, 0.10963690963690906, -0.022222222222222105}
    };

    double rms_prev = 1e10;
    for (size_t i = 0; i < PP.size(); i++) {
        size_t N = i + 1;
        vector_fp p(N);
        double rms = polyfit(7, i, x.data(), y.data(), nullptr, p.data());
        EXPECT_LT(rms, rms_prev);
        rms_prev = rms;
        for (size_t j = 0; j < N; j++) {
            EXPECT_NEAR(PP[i][j], p[j], 1e-14);
        }
    }
}

TEST(Polyfit, weighted)
{
    vector_fp x{-1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0};
    vector_fp y{0.6, 1.0, 0.8, 0.4, -0.1, -0.5, -1};
    vector_fp w{25, 1, 1, 1, 1, 1, 100}; // these are the squares of Numpy's weights

    // Coefficients calculated using Numpy's polyfit function for polynomials
    // of degrees 0 - 5.
    std::vector<vector_fp> PP{
        {-0.64153846153846139},
        {0.24582603619381152, -0.41199065966141246},
        {0.64897277949822718, -0.10796777523450461, -0.14749113594542437},
        {1.0095165556633916, -0.22435606362053356, -0.51254844673169053,
         0.12135217568551074},
        {1.0121717322829622, -0.23147507683766383, -0.51492677362711337,
         0.12728869689006062, -0.0014837700620763492},
        {0.998127784554808, -0.093474983983779111, -0.61196784469972776,
         0.011482911646053995, 0.10962944760868476, -0.022222284629403764}
    };

    double rms_prev = 1e10;
    for (size_t i = 0; i < PP.size(); i++) {
        size_t N = i + 1;
        vector_fp p(N);
        double rms = polyfit(7, i, x.data(), y.data(), w.data(), p.data());
        EXPECT_LT(rms, rms_prev);
        rms_prev = rms;
        for (size_t j = 0; j < N; j++) {
             EXPECT_NEAR(PP[i][j], p[j], 1e-14);
        }
    }
}
