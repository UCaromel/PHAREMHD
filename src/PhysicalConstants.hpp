#ifndef PHYSICAL_CONSTANTS_HPP_
#define PHYSICAL_CONSTANTS_HPP_

struct Consts {
    const double sigmaCFL;
    const double gam;
    const double eta; // Resistivity
    const double nu;  // Hyper resistivity

    Consts(double sigmaCFL = 0.8, double gam = 5.0 / 3.0, double eta = 0.0, double nu = 0.0)
        : sigmaCFL(sigmaCFL), gam(gam), eta(eta), nu(nu) {}
};

class PhysicalConstants {
public:
    static PhysicalConstants& getInstance();

    void init(const Consts& consts) {
        sigmaCFL = consts.sigmaCFL;
        gam = consts.gam;
        eta = consts.eta;
        nu = consts.nu;
    }

    double sigmaCFL;
    double gam;
    double eta;
    double nu;

private:
    PhysicalConstants() = default;
};

#endif // PHYSICAL_CONSTANTS_HPP_
