#ifndef RECONSTRUCTED_VALUES_HPP_
#define RECONSTRUCTED_VALUES_HPP_

struct ReconstructedValues {
    double rho;
    double vx;
    double vy;
    double vz;
    double Bx;
    double By;
    double Bz;
    double P;

    ReconstructedValues operator+(const ReconstructedValues& rhs) const {
        ReconstructedValues result;
        result.rho = rho + rhs.rho;
        result.vx = vx + rhs.vx;
        result.vy = vy + rhs.vy;
        result.vz = vz + rhs.vz;
        result.Bx = Bx + rhs.Bx;
        result.By = By + rhs.By;
        result.Bz = Bz + rhs.Bz;
        result.P = P + rhs.P;
        return result;
    }

    ReconstructedValues operator-(const ReconstructedValues& rhs) const {
        ReconstructedValues result;
        result.rho = rho - rhs.rho;
        result.vx = vx - rhs.vx;
        result.vy = vy - rhs.vy;
        result.vz = vz - rhs.vz;
        result.Bx = Bx - rhs.Bx;
        result.By = By - rhs.By;
        result.Bz = Bz - rhs.Bz;
        result.P = P - rhs.P;
        return result;
    }

    ReconstructedValues operator*(double scalar) const {
        ReconstructedValues result;
        result.rho = rho * scalar;
        result.vx = vx * scalar;
        result.vy = vy * scalar;
        result.vz = vz * scalar;
        result.Bx = Bx * scalar;
        result.By = By * scalar;
        result.Bz = Bz * scalar;
        result.P = P * scalar;
        return result;
    }

    friend ReconstructedValues operator*(double scalar, const ReconstructedValues& rhs) {
        return rhs * scalar;
    }

    ReconstructedValues operator=(const ReconstructedValues& other) {
        rho = other.rho;
        vx = other.vx;
        vy = other.vy;
        vz = other.vz;
        Bx = other.Bx;
        By = other.By;
        Bz = other.Bz;
        P = other.P;
        return *this;
    }
};

#endif //RECONSTRUCTED_VALUES_HPP_