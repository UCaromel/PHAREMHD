#include "Interface.hpp"

static const PhysicalConstants& pc = PhysicalConstants::getInstance(); 

ReconstructedValues ComputeFluxVector(const ReconstructedValues& u, Dir dir) {
    ReconstructedValues flux;
    double GeneralisedPressure = u.P + 0.5*(u.Bx * u.Bx + u.By * u.By + u.Bz * u.Bz);

    if (dir == Dir::X) {
        flux.rho = u.rho * u.vx;
        flux.vx = u.rho * u.vx * u.vx + GeneralisedPressure - u.Bx * u.Bx;
        flux.vy = u.rho * u.vx * u.vy - u.Bx * u.By;
        flux.vz = u.rho * u.vx * u.vz - u.Bx * u.Bz;
        flux.Bx = 0.0;
        flux.By = u.By * u.vx - u.vy * u.Bx;
        flux.Bz = u.Bz * u.vx - u.vz * u.Bx;
        flux.P = (EosEtot(u) + GeneralisedPressure)*u.vx - u.Bx*(u.vx*u.Bx + u.vy*u.By + u.vz*u.Bz);
    } else if (dir == Dir::Y) {
        flux.rho = u.rho * u.vy;
        flux.vx = u.rho * u.vy * u.vx - u.By * u.Bx;
        flux.vy = u.rho * u.vy * u.vy + GeneralisedPressure - u.By * u.By;
        flux.vz = u.rho * u.vy * u.vz - u.By * u.Bz;
        flux.Bx = u.Bx * u.vy - u.By * u.vx;
        flux.By = 0.0;
        flux.Bz = u.Bz * u.vy - u.vz * u.By;
        flux.P = (EosEtot(u) + GeneralisedPressure)*u.vy - u.By*(u.vx*u.Bx + u.vy*u.By + u.vz*u.Bz);
    }
    return flux;
}

void AddNonIdealFlux(ReconstructedValues& f, const ReconstructedValues& u, double Jx, double Jy, double Jz, double LaplJx, double LaplJy, double LaplJz, OptionalPhysics OptP, Dir dir){
    if (OptP == OptionalPhysics::HallRes) {
        if (dir == Dir::X) {
            f.By += - (1.0/u.rho) * (Jx * u.By - Jy * u.Bx) *1.0;
            //f.By += - pc.eta * Jz;
            //f.By -= - pc.nu * LaplJz;
            f.Bz += (1.0/u.rho) * (Jz * u.Bx - Jx * u.Bz) *1.0;
            //f.Bz += pc.eta * Jy;
            //f.Bz -= pc.nu * LaplJy;

            f.P += (1.0/u.rho) * ((u.Bx * Jx + u.By * Jy + u.Bz * Jz)*u.Bx - (u.Bx * u.Bx + u.By * u.By + u.Bz * u.Bz) * Jx) *1.0;
            //f.P += pc.eta * (Jy * u.Bz - Jz * u.By);
            //f.P -= pc.nu * (LaplJy * u.Bz - LaplJz * u.By);
        } else if (dir == Dir::Y) {
            f.Bx += (1.0/u.rho) * (Jx * u.By - Jy * u.Bx);
            //f.Bx += pc.eta * Jz;
            //f.Bx -= pc.nu * LaplJz;
            f.Bz += - (1.0/u.rho) * (Jy * u.Bz - Jz * u.By);
            //f.Bz += - pc.eta * Jx;
            //f.Bz -= - pc.nu * LaplJx;

            f.P += (1.0/u.rho) * ((u.Bx * Jx + u.By * Jy + u.Bz * Jz)*u.By - (u.Bx * u.Bx + u.By * u.By + u.Bz * u.Bz) * Jy);
            //f.P += pc.eta * (Jz * u.Bx - Jx * u.Bz);
            //f.P -= pc.nu * (LaplJz * u.Bx - LaplJx * u.Bz);
        }
    }
}

template<typename Func>
std::pair<std::pair<double, double>, std::pair<double, double>> ComputeRiemannJ(const std::vector<std::vector<double>>& J, Func AVERAGING, int i, int j, double Dx, double Dy, Reconstruction rec, Slope sl, Dir Dir) {
    double JL;
    double JR;
    double LaplJL;
    double LaplJR;

    // J has one more ghost cell
    i++;
    j++;

    if (rec == Reconstruction::Constant) {
        // for constant reconstruction, the reconstructed values are the left and right cell centered values next to the interface
        // for x derivatives
        double J_2x = AVERAGING(J, i-2, j);
        double J_1x = AVERAGING(J, i-1, j);
        double Jx = AVERAGING(J, i, j);
        double J1x = AVERAGING(J, i+1, j);

        // for y derivatives
        double J_2y = AVERAGING(J, i, j-2);
        double J_1y = AVERAGING(J, i, j-1);
        double Jy = AVERAGING(J, i, j);
        double J1y = AVERAGING(J, i, j+1);

        LaplJL = (J_2x - 2.0 * J_1x + Jx) / (Dx * Dx) + (J_2y - 2.0 * J_1y + Jy) / (Dy * Dy);
        LaplJR = (J_1x - 2.0 * Jx + J1x) / (Dx * Dx) + (J_1y - 2.0 * Jy + J1y) / (Dy * Dy);

        if (Dir == Dir::X) {
            JL = J_1x;
            JR = Jx;
        } else if (Dir == Dir::Y) {
            JL = J_1y;
            JR = Jy;
        }

    } else if (rec == Reconstruction::Linear) {
        SLFunctionDouble ChosenSLDouble = getSlopeLimiterDouble(sl);

        // for x derivatives
        double J_3x = AVERAGING(J, i-3, j);
        double J_2x = AVERAGING(J, i-2, j);
        double J_1x = AVERAGING(J, i-1, j);
        double Jx = AVERAGING(J, i, j);
        double J1x = AVERAGING(J, i+1, j);
        double J2x = AVERAGING(J, i+2, j);

        // for y derivatives
        double J_3y = AVERAGING(J, i, j-3);
        double J_2y = AVERAGING(J, i, j-2);
        double J_1y = AVERAGING(J, i, j-1);
        double Jy = AVERAGING(J, i, j);
        double J1y = AVERAGING(J, i, j+1);
        double J2y = AVERAGING(J, i, j+2);

        // Slopes
        double Di_1x = ChosenSLDouble(J_1x - J_2x, J_2x - J_3x);
        double Dix = ChosenSLDouble(Jx - J_1x, J_1x - J_2x);
        double Di1x = ChosenSLDouble(J1x - Jx, Jx - J_1x);
        double Di2x = ChosenSLDouble(J2x - J1x, J1x - Jx);

        double Di_1y = ChosenSLDouble(J_1y - J_2y, J_2y - J_3y);
        double Diy = ChosenSLDouble(Jy - J_1y, J_1y - J_2y);
        double Di1y = ChosenSLDouble(J1y - Jy, Jy - J_1y);
        double Di2y = ChosenSLDouble(J2y - J1y, J1y - Jy);

        // for L interface Laplacian
        double J_1xR = J_2x + 0.5 * Di_1x;
        double JxR = J_1x + 0.5 * Dix;
        double J1xR = Jx + 0.5 * Di1x;

        double J_1yR = J_2y + 0.5 * Di_1y;
        double JyR = J_1y + 0.5 * Diy;
        double J1yR = Jy + 0.5 * Di1y;

        // for R interface Laplacian
        double JxL = J_1x - 0.5 * Dix;
        double J1xL = Jx - 0.5 * Di1x;
        double J2xL = J1x - 0.5 * Di2x;

        double JyL = J_1y - 0.5 * Diy;
        double J1yL = Jy - 0.5 * Di1y;
        double J2yL = J1y - 0.5 * Di2y;

        LaplJL = (J_1xR - 2.0 * JxR + J1xR) / (Dx * Dx) + (J_1yR - 2.0 * JyR + J1yR) / (Dy * Dy);
        LaplJR = (JxL - 2.0 * J1xL + J2xL) / (Dx * Dx) + (JyL - 2.0 * J1yL + J2yL) / (Dy * Dy);

        if (Dir == Dir::X) {
            JL = JxR;
            JR = J1xL;
        } else if (Dir == Dir::Y) {
            JL = JyR;
            JR = J1yL;
        }

    } else {
        throw std::invalid_argument("Invalid reconstruction type");
    }

    return std::make_pair(std::make_pair(JL, LaplJL), std::make_pair(JR, LaplJR));
}

Interface::Interface() = default;

Interface::Interface(const PrimitiveVariables& P_cc /* Assuming ghost cells are added */, int i /* (0 to nx) + nghost */, int j /* (0 to ny) + nghost */, double Dx, double Dy, int nghost, Reconstruction rec, Slope sl, OptionalPhysics OptP,  Dir dir) {
    // For riemann solver
    OP = OptP;
        
    if (rec == Reconstruction::Constant) {
        if (dir == Dir::X) {
            uL = P_cc(i-1,j);
            uR = P_cc(i,j);
        } else if (dir == Dir::Y) {
            uL = P_cc(i,j-1);
            uR = P_cc(i,j);
        }
    }

    if (rec == Reconstruction::Linear){
        if (nghost < 2) {
            throw std::invalid_argument("nghost must be at least 2 for linear reconstruction");
        }

        SLFunction ChosenSL = getSlopeLimiter(sl);

        ReconstructedValues ui_1;
        ReconstructedValues ui;
        ReconstructedValues ui1;
        ReconstructedValues ui2;

        if (dir == Dir::X){
            ui_1 = P_cc(i-2,j);
            ui = P_cc(i-1,j);
            ui1 = P_cc(i,j);
            ui2 = P_cc(i+1,j);
        } else if (dir == Dir::Y){
            ui_1 = P_cc(i,j-2);
            ui = P_cc(i,j-1);
            ui1 = P_cc(i,j);
            ui2 = P_cc(i,j+1);
        }

        ReconstructedValues Dui_12 = ui - ui_1;
        ReconstructedValues Dui12 = ui1 - ui;
        ReconstructedValues Di = ChosenSL(Dui12, Dui_12);

        ReconstructedValues Dui1_12 = Dui12;
        ReconstructedValues Dui112 = ui2 - ui1;
        ReconstructedValues Di1 = ChosenSL(Dui112, Dui1_12);

        ReconstructedValues uiR = ui + 0.5 * Di;
        ReconstructedValues ui1L = ui1 - 0.5 * Di1;

        uL = uiR;
        uR = ui1L;
    }

    fL = ComputeFluxVector(uL, dir);
    fR = ComputeFluxVector(uR, dir);

    if (OptP == OptionalPhysics::HallResHyper) {
        auto [Lx, Rx] = ComputeRiemannJ(P_cc.Jx, AVERAGEY, i, j, Dx, Dy, rec, sl, dir);
        auto [JxL, LaplJxL] = Lx;
        auto [JxR, LaplJxR] = Rx;

        auto [Ly, Ry] = ComputeRiemannJ(P_cc.Jy, AVERAGEX, i, j, Dx, Dy, rec, sl, dir);
        auto [JyL, LaplJyL] = Ly;
        auto [JyR, LaplJyR] = Ry;

        auto [Lz, Rz] = ComputeRiemannJ(P_cc.Jz, AVERAGEXY, i, j, Dx, Dy, rec, sl, dir);
        auto [JzL, LaplJzL] = Lz;
        auto [JzR, LaplJzR] = Rz;

        AddNonIdealFlux(fL, uL, JxL, JyL, JzL, LaplJxL, LaplJyL, LaplJzL, OptP, dir);
        AddNonIdealFlux(fR, uR, JxR, JyR, JzR, LaplJxR, LaplJyR, LaplJzR, OptP, dir);
    }

    double c0L = std::sqrt((pc.gam*uL.P)/uL.rho); // Sound speeds
    double c0R = std::sqrt((pc.gam*uR.P)/uR.rho); 

    double caxL = std::sqrt((uL.Bx*uL.Bx)/(uL.rho)); // Alfven speeds in x
    double caxR = std::sqrt((uR.Bx*uR.Bx)/(uR.rho));
    double cayL = std::sqrt((uL.By*uL.By)/(uL.rho)); // Alfven speeds in y
    double cayR = std::sqrt((uR.By*uR.By)/(uR.rho));

    double caL = std::sqrt((uL.Bx*uL.Bx + uL.By*uL.By + uL.Bz*uL.Bz)/(uL.rho)); // Alfven speeds
    double caR = std::sqrt((uR.Bx*uR.Bx + uR.By*uR.By + uR.Bz*uR.Bz)/(uR.rho));

    cfastxL = std::sqrt((c0L*c0L + caL*caL)*0.5 + (std::sqrt((c0L*c0L + caL*caL)*(c0L*c0L + caL*caL) - 4*c0L*c0L*caxL*caxL))*0.5); // Fast magnetosonic speeds in x
    cfastxR = std::sqrt((c0R*c0R + caR*caR)*0.5 + (std::sqrt((c0R*c0R + caR*caR)*(c0R*c0R + caR*caR) - 4*c0R*c0R*caxR*caxR))*0.5);
    cfastyL = std::sqrt((c0L*c0L + caL*caL)*0.5 + (std::sqrt((c0L*c0L + caL*caL)*(c0L*c0L + caL*caL) - 4*c0L*c0L*cayL*cayL))*0.5); // Fast magnetosonic speeds in y
    cfastyR = std::sqrt((c0R*c0R + caR*caR)*0.5 + (std::sqrt((c0R*c0R + caR*caR)*(c0R*c0R + caR*caR) - 4*c0R*c0R*cayR*cayR))*0.5);

    if (OptP == OptionalPhysics::HallResHyper) {
        double vwx = M_PI * (std::sqrt(1 + 0.25/(Dx*Dx)) + 0.5/Dx);
        double vwy = M_PI * (std::sqrt(1 + 0.25/(Dy*Dy)) + 0.5/Dy);
        cwxL = std::sqrt(uL.Bx*uL.Bx + uL.By*uL.By + uL.Bz*uL.Bz) / (uL.rho) * vwx;
        cwyL = std::sqrt(uL.Bx*uL.Bx + uL.By*uL.By + uL.Bz*uL.Bz) / (uL.rho) * vwy;
        cwxR = std::sqrt(uR.Bx*uR.Bx + uR.By*uR.By + uR.Bz*uR.Bz) / (uR.rho) * vwx;
        cwyR = std::sqrt(uR.Bx*uR.Bx + uR.By*uR.By + uR.Bz*uR.Bz) / (uR.rho) * vwy;
        //double hallxL = std::sqrt(uL.Bx*uL.Bx + uL.By*uL.By + uL.Bz*uL.Bz) / (2.0*(uL.rho)*Dx);
        //double hallyL = std::sqrt(uL.Bx*uL.Bx + uL.By*uL.By + uL.Bz*uL.Bz) / (2.0*(uL.rho)*Dy);
        //double hallxR = std::sqrt(uR.Bx*uR.Bx + uR.By*uR.By + uR.Bz*uR.Bz) / (2.0*(uR.rho)*Dx);
        //double hallyR = std::sqrt(uR.Bx*uR.Bx + uR.By*uR.By + uR.Bz*uR.Bz) / (2.0*(uR.rho)*Dy);
        //cwxL = std::fabs(hallxL) + std::sqrt(hallxL*hallxL + caxL*caxL);
        //cwyL = std::fabs(hallyL) + std::sqrt(hallyL*hallyL + cayL*cayL);
        //cwxR = std::fabs(hallxR) + std::sqrt(hallxR*hallxR + caxR*caxR);
        //cwyR = std::fabs(hallyR) + std::sqrt(hallyR*hallyR + cayR*cayR);

        if(dir == Dir::X){
            SLb = std::min(uL.vx - cfastxL - cwxL, uR.vx - cfastxR - cwxR);
            SRb = std::max(uL.vx + cfastxL + cwxL, uR.vx + cfastxR + cwxR);
            // For rusanov :
            Splusb = std::max(std::abs(uL.vx) + cfastxL + cwxL, std::abs(uR.vx) + cfastxR + cwxR);
        }else if(dir == Dir::Y){
            SLb = std::min(uL.vy - cfastyL - cwyL, uR.vy - cfastyR - cwyR);
            SRb = std::max(uL.vy + cfastyL + cwyL, uR.vy + cfastyR + cwyR);
            //  For rusanov :
            Splusb = std::max(std::abs(uL.vy) + cfastyL + cwyL, std::abs(uR.vy) + cfastyR + cwyR);
        }
    }

    // Wave speeds
    if(dir == Dir::X){
        SL = std::min(uL.vx - cfastxL, uR.vx - cfastxR);
        SR = std::max(uL.vx + cfastxL, uR.vx + cfastxR);
        // For rusanov :
        Splus = std::max(std::abs(uL.vx) + cfastxL, std::abs(uR.vx) + cfastxR);
    }else if(dir == Dir::Y){
        SL = std::min(uL.vy - cfastyL, uR.vy - cfastyR);
        SR = std::max(uL.vy + cfastyL, uR.vy + cfastyR);
        //  For rusanov :
        Splus = std::max(std::abs(uL.vy) + cfastyL, std::abs(uR.vy) + cfastyR);
    }

    // Pass uL and uR in conservative form
    uL.vx = uL.rho * uL.vx;
    uL.vy = uL.rho * uL.vy;
    uL.vz = uL.rho * uL.vz;
    uL.P = EosEtot(uL);
    uR.vx = uR.rho * uR.vx;
    uR.vy = uR.rho * uR.vy;
    uR.vz = uR.rho * uR.vz;
    uR.P = EosEtot(uR);
}

Interface::~Interface() = default;