#ifndef RECONSTRUCTED_VALUES_HPP_
#define RECONSTRUCTED_VALUES_HPP_

struct ReconstructedValues {
    double rho;
    double vx;
    double vy;
    double vz;
    double P;
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
        result.P = P + rhs.P;
        result.Bx = Bx + rhs.Bx;
        result.By = By + rhs.By;
        result.Bz = Bz + rhs.Bz;
        return result;
    }

    ReconstructedValues operator-(const ReconstructedValues& rhs) const {
        ReconstructedValues result;
        result.rho = rho - rhs.rho;
        result.vx = vx - rhs.vx;
        result.vy = vy - rhs.vy;
        result.vz = vz - rhs.vz;
        result.P = P - rhs.P;
        result.Bx = Bx - rhs.Bx;
        result.By = By - rhs.By;
        result.Bz = Bz - rhs.Bz;
        return result;
    }

    ReconstructedValues operator*(double scalar) const {
        ReconstructedValues result;
        result.rho = rho * scalar;
        result.vx = vx * scalar;
        result.vy = vy * scalar;
        result.vz = vz * scalar;
        result.P = P * scalar;
        result.Bx = Bx * scalar;
        result.By = By * scalar;
        result.Bz = Bz * scalar;
        return result;
    }
};

ReconstructedValues operator*(double scalar, const ReconstructedValues& rhs) {
    return rhs * scalar;
}

#endif //RECONSTRUCTED_VALUES_HPP_