#include "Interface.hpp"
#include "ComputeRiemannJ.hpp"
#include "Enums.hpp"
#include "ReconstructedValues.hpp"

static const PhysicalConstants &pc = PhysicalConstants::getInstance();

ReconstructedValues ComputeFluxVector(const ReconstructedValues &u, Dir dir) {
  ReconstructedValues flux;
  double GeneralisedPressure =
      u.P + 0.5 * (u.Bx * u.Bx + u.By * u.By + u.Bz * u.Bz);

  if (dir == Dir::X) {
    flux.rho = u.rho * u.vx;
    flux.vx = u.rho * u.vx * u.vx + GeneralisedPressure - u.Bx * u.Bx;
    flux.vy = u.rho * u.vx * u.vy - u.Bx * u.By;
    flux.vz = u.rho * u.vx * u.vz - u.Bx * u.Bz;
    flux.Bx = 0.0;
    flux.By = u.By * u.vx - u.vy * u.Bx;
    flux.Bz = u.Bz * u.vx - u.vz * u.Bx;
    flux.P = (EosEtot(u) + GeneralisedPressure) * u.vx -
             u.Bx * (u.vx * u.Bx + u.vy * u.By + u.vz * u.Bz);
  } else if (dir == Dir::Y) {
    flux.rho = u.rho * u.vy;
    flux.vx = u.rho * u.vy * u.vx - u.By * u.Bx;
    flux.vy = u.rho * u.vy * u.vy + GeneralisedPressure - u.By * u.By;
    flux.vz = u.rho * u.vy * u.vz - u.By * u.Bz;
    flux.Bx = u.Bx * u.vy - u.By * u.vx;
    flux.By = 0.0;
    flux.Bz = u.Bz * u.vy - u.vz * u.By;
    flux.P = (EosEtot(u) + GeneralisedPressure) * u.vy -
             u.By * (u.vx * u.Bx + u.vy * u.By + u.vz * u.Bz);
  }
  return flux;
}

void AddNonIdealFlux(ReconstructedValues &f, const ReconstructedValues &u,
                     double Jx, double Jy, double Jz, double LaplJx,
                     double LaplJy, double LaplJz, OptionalPhysics OptP,
                     Dir dir) {
  if (OptP == OptionalPhysics::HallResHyper) {
    if (dir == Dir::X) {
      f.By += -(1.0 / u.rho) * (Jx * u.By - Jy * u.Bx) * 1.0;
      // f.By += - pc.eta * Jz;
      // f.By -= - pc.nu * LaplJz;
      f.Bz += (1.0 / u.rho) * (Jz * u.Bx - Jx * u.Bz) * 1.0;
      // f.Bz += pc.eta * Jy;
      // f.Bz -= pc.nu * LaplJy;

      f.P += (1.0 / u.rho) *
             ((u.Bx * Jx + u.By * Jy + u.Bz * Jz) * u.Bx -
              (u.Bx * u.Bx + u.By * u.By + u.Bz * u.Bz) * Jx) *
             1.0;
      // f.P += pc.eta * (Jy * u.Bz - Jz * u.By);
      // f.P -= pc.nu * (LaplJy * u.Bz - LaplJz * u.By);
    } else if (dir == Dir::Y) {
      f.Bx += (1.0 / u.rho) * (Jx * u.By - Jy * u.Bx);
      // f.Bx += pc.eta * Jz;
      // f.Bx -= pc.nu * LaplJz;
      f.Bz += -(1.0 / u.rho) * (Jy * u.Bz - Jz * u.By);
      // f.Bz += - pc.eta * Jx;
      // f.Bz -= - pc.nu * LaplJx;

      f.P += (1.0 / u.rho) * ((u.Bx * Jx + u.By * Jy + u.Bz * Jz) * u.By -
                              (u.Bx * u.Bx + u.By * u.By + u.Bz * u.Bz) * Jy);
      // f.P += pc.eta * (Jz * u.Bx - Jx * u.Bz);
      // f.P -= pc.nu * (LaplJz * u.Bx - LaplJx * u.Bz);
    }
  }
}

Interface::Interface() = default;

Interface::Interface(
    const PrimitiveVariables &P_cc /* Assuming ghost cells are added */,
    int i /* (0 to nx) + nghost */, int j /* (0 to ny) + nghost */, double Dx,
    double Dy, int nghost, Reconstruction rec, Slope sl, OptionalPhysics OptP,
    Dir dir) {
  // For riemann solver
  OP = OptP;

  if (rec == Reconstruction::Constant) {
    if (dir == Dir::X) {
      uL = P_cc(i - 1, j);
      uR = P_cc(i, j);
    } else if (dir == Dir::Y) {
      uL = P_cc(i, j - 1);
      uR = P_cc(i, j);
    }
  }

  if (rec == Reconstruction::Linear) {
    if (nghost < 2) {
      throw std::invalid_argument(
          "nghost must be at least 2 for linear reconstruction");
    }

    SLFunction ChosenSL = getSlopeLimiter(sl);

    ReconstructedValues ui_2;
    ReconstructedValues ui_1;
    ReconstructedValues ui;
    ReconstructedValues ui1;

    if (dir == Dir::X) {
      ui_2 = P_cc(i - 2, j);
      ui_1 = P_cc(i - 1, j);
      ui = P_cc(i, j);
      ui1 = P_cc(i + 1, j);
    } else if (dir == Dir::Y) {
      ui_2 = P_cc(i, j - 2);
      ui_1 = P_cc(i, j - 1);
      ui = P_cc(i, j);
      ui1 = P_cc(i, j + 1);
    }

    ReconstructedValues Dui_12 = ui_1 - ui_2;
    ReconstructedValues Dui12 = ui - ui_1;
    ReconstructedValues Di_1 = ChosenSL(Dui12, Dui_12);

    ReconstructedValues Dui1_12 = Dui12;
    ReconstructedValues Dui112 = ui1 - ui;
    ReconstructedValues Di = ChosenSL(Dui112, Dui1_12);

    ReconstructedValues ui_1R = ui_1 + 0.5 * Di_1;
    ReconstructedValues uiL = ui - 0.5 * Di;

    uL = ui_1R;
    uR = uiL;
  }

  if (rec == Reconstruction::WENO3) {
    if (nghost < 2) {
      throw std::invalid_argument(
          "nghost must be at least 2 for linear reconstruction");
    }

    ReconstructedValues ui_2;
    ReconstructedValues ui_1;
    ReconstructedValues ui;
    ReconstructedValues ui1;

    if (dir == Dir::X) {
      ui_2 = P_cc(i - 2, j);
      ui_1 = P_cc(i - 1, j);
      ui = P_cc(i, j);
      ui1 = P_cc(i + 1, j);
    } else if (dir == Dir::Y) {
      ui_2 = P_cc(i, j - 2);
      ui_1 = P_cc(i, j - 1);
      ui = P_cc(i, j);
      ui1 = P_cc(i, j + 1);
    }

    double eps = 1.e-6;

    double dL0 = 1. / 3.;
    double dL1 = 2. / 3.;
    double dR0 = 2. / 3.;
    double dR1 = 1. / 3.;

    ReconstructedValues beta0 = (ui_1 - ui_2) * (ui_1 - ui_2);
    ReconstructedValues beta1 = (ui - ui_1) * (ui - ui_1);
    ReconstructedValues beta2 = (ui1 - ui) * (ui1 - ui);

    ReconstructedValues alphaL0 = dL0 / ((beta0 + eps) * (beta0 + eps));
    ReconstructedValues alphaL1 = dL1 / ((beta1 + eps) * (beta1 + eps));
    ReconstructedValues alphaR0 = dR0 / ((beta1 + eps) * (beta1 + eps));
    ReconstructedValues alphaR1 = dR1 / ((beta2 + eps) * (beta2 + eps));

    ReconstructedValues wL0 = alphaL0 / (alphaL0 + alphaL1);
    ReconstructedValues wL1 = alphaL1 / (alphaL0 + alphaL1);
    ReconstructedValues wR0 = alphaR0 / (alphaR0 + alphaR1);
    ReconstructedValues wR1 = alphaR1 / (alphaR0 + alphaR1);

    uL = wL0 * (-0.5 * ui_2 + 1.5 * ui_1) + wL1 * (0.5 * ui_1 + 0.5 * ui);
    uR = wR0 * (0.5 * ui + 0.5 * ui_1) + wR1 * (-0.5 * ui1 + 1.5 * ui);
  }

  if (rec == Reconstruction::WENO5) {
    if (nghost < 3) {
      throw std::invalid_argument("nghost must be at least 3 for WENO5");
    }

    ReconstructedValues ui_3;
    ReconstructedValues ui_2;
    ReconstructedValues ui_1;
    ReconstructedValues ui;
    ReconstructedValues ui1;
    ReconstructedValues ui2;

    if (dir == Dir::X) {
      ui_3 = P_cc(i - 3, j);
      ui_2 = P_cc(i - 2, j);
      ui_1 = P_cc(i - 1, j);
      ui = P_cc(i, j);
      ui1 = P_cc(i + 1, j);
      ui2 = P_cc(i + 2, j);
    } else if (dir == Dir::Y) {
      ui_3 = P_cc(i, j - 3);
      ui_2 = P_cc(i, j - 2);
      ui_1 = P_cc(i, j - 1);
      ui = P_cc(i, j);
      ui1 = P_cc(i, j + 1);
      ui2 = P_cc(i, j + 2);
    }

    double eps = 1.e-6;

    double dL0 = 1. / 10.;
    double dL1 = 3. / 5.;
    double dL2 = 3. / 10.;
    double dR0 = 3. / 10.;
    double dR1 = 3. / 5.;
    double dR2 = 1. / 10.;

    ReconstructedValues betaL0 =
        (13. / 12.) * (ui_3 - 2 * ui_2 + ui_1) * (ui_3 - 2 * ui_2 + ui_1) +
        (1. / 4.) * (ui_3 - 4 * ui_2 + 3 * ui_1) * (ui_3 - 4 * ui_2 + 3 * ui_1);
    ReconstructedValues betaL1 =
        (13. / 12.) * (ui_2 - 2 * ui_1 + ui) * (ui_2 - 2 * ui_1 + ui) +
        (1. / 4.) * (ui_2 - ui) * (ui_2 - ui);
    ReconstructedValues betaL2 =
        (13. / 12.) * (ui_1 - 2 * ui + ui1) * (ui_1 - 2 * ui + ui1) +
        (1. / 4.) * (3 * ui_1 - 4 * ui + ui1) * (3 * ui_1 - 4 * ui + ui1);

    ReconstructedValues betaR0 =
        (13. / 12.) * (ui_2 - 2 * ui_1 + ui) * (ui_2 - 2 * ui_1 + ui) +
        (1. / 4.) * (ui_2 - 4 * ui_1 + 3 * ui) * (ui_2 - 4 * ui_1 + 3 * ui);
    ReconstructedValues betaR1 =
        (13. / 12.) * (ui_1 - 2 * ui + ui1) * (ui_1 - 2 * ui + ui1) +
        (1. / 4.) * (ui_1 - ui1) * (ui_1 - ui1);
    ReconstructedValues betaR2 =
        (13. / 12.) * (ui - 2 * ui1 + ui2) * (ui - 2 * ui1 + ui2) +
        (1. / 4.) * (3 * ui - 4 * ui1 + ui2) * (3 * ui - 4 * ui1 + ui2);

    ReconstructedValues alphaL0 = dL0 / ((betaL0 + eps) * (betaL0 + eps));
    ReconstructedValues alphaL1 = dL1 / ((betaL1 + eps) * (betaL1 + eps));
    ReconstructedValues alphaL2 = dL2 / ((betaL2 + eps) * (betaL2 + eps));

    ReconstructedValues alphaR0 = dR0 / ((betaR0 + eps) * (betaR0 + eps));
    ReconstructedValues alphaR1 = dR1 / ((betaR1 + eps) * (betaR1 + eps));
    ReconstructedValues alphaR2 = dR2 / ((betaR2 + eps) * (betaR2 + eps));

    ReconstructedValues wL0 = alphaL0 / (alphaL0 + alphaL1 + alphaL2);
    ReconstructedValues wL1 = alphaL1 / (alphaL0 + alphaL1 + alphaL2);
    ReconstructedValues wL2 = alphaL2 / (alphaL0 + alphaL1 + alphaL2);
    ReconstructedValues wR0 = alphaR0 / (alphaR0 + alphaR1 + alphaR2);
    ReconstructedValues wR1 = alphaR1 / (alphaR0 + alphaR1 + alphaR2);
    ReconstructedValues wR2 = alphaR2 / (alphaR0 + alphaR1 + alphaR2);

    uL = wL0 * ((1. / 3.) * ui_3 - (7. / 6.) * ui_2 + (11. / 6.) * ui_1) +
         wL1 * (-(1. / 6.) * ui_2 + (5. / 6.) * ui_1 + (1. / 3.) * ui) +
         wL2 * ((1. / 3.) * ui_1 + (5. / 6.) * ui - (1. / 6.) * ui1);
    uR = wR0 * ((1. / 3.) * ui + (5. / 6.) * ui_1 - (1. / 6.) * ui_2) +
         wR1 * (-(1. / 6.) * ui1 + (5. / 6.) * ui + (1. / 3.) * ui_1) +
         wR2 * ((1. / 3.) * ui2 - (7. / 6.) * ui1 + (11. / 6.) * ui);
  }

  if (rec == Reconstruction::WENO7) {
    if (nghost < 4) {
      throw std::invalid_argument("nghost must be at least 4 for WENO7");
    }

    ReconstructedValues ui_4;
    ReconstructedValues ui_3;
    ReconstructedValues ui_2;
    ReconstructedValues ui_1;
    ReconstructedValues ui;
    ReconstructedValues ui1;
    ReconstructedValues ui2;
    ReconstructedValues ui3;

    if (dir == Dir::X) {
      ui_4 = P_cc(i - 4, j);
      ui_3 = P_cc(i - 3, j);
      ui_2 = P_cc(i - 2, j);
      ui_1 = P_cc(i - 1, j);
      ui = P_cc(i, j);
      ui1 = P_cc(i + 1, j);
      ui2 = P_cc(i + 2, j);
      ui3 = P_cc(i + 3, j);
    } else if (dir == Dir::Y) {
      ui_4 = P_cc(i, j - 4);
      ui_3 = P_cc(i, j - 3);
      ui_2 = P_cc(i, j - 2);
      ui_1 = P_cc(i, j - 1);
      ui = P_cc(i, j);
      ui1 = P_cc(i, j + 1);
      ui2 = P_cc(i, j + 2);
      ui3 = P_cc(i, j + 3);
    }

    double eps = 1.e-6;

    double dL0 = 1. / 35.;
    double dL1 = 12. / 35.;
    double dL2 = 18. / 35.;
    double dL3 = 4. / 35.;
    double dR0 = 4. / 35.;
    double dR1 = 18. / 35.;
    double dR2 = 12. / 35.;
    double dR3 = 1. / 35.;

    ReconstructedValues betaL0 =
        ui_4 * (547 * ui_4 - 3882 * ui_3 + 4642 * ui_2 - 1854 * ui_1) +
        ui_3 * (7043 * ui_3 - 17246 * ui_2 + 7042 * ui_1) +
        ui_2 * (11003 * ui_2 - 9402 * ui_1) + ui_1 * (2107 * ui_1);
    ReconstructedValues betaL1 =
        ui_3 * (267 * ui_3 - 1642 * ui_2 + 1602 * ui_1 - 494 * ui) +
        ui_2 * (2843 * ui_2 - 5966 * ui_1 + 1922 * ui) +
        ui_1 * (3443 * ui_1 - 2522 * ui) + ui * (547 * ui);
    ReconstructedValues betaL2 =
        ui_2 * (547 * ui_2 - 2522 * ui_1 + 1922 * ui - 494 * ui1) +
        ui_1 * (3443 * ui_1 - 5966 * ui + 1602 * ui1) +
        ui * (2843 * ui - 1642 * ui1) + ui1 * (267 * ui1);
    ReconstructedValues betaL3 =
        ui_1 * (2107 * ui_1 - 9402 * ui + 7042 * ui1 - 1854 * ui2) +
        ui * (11003 * ui - 17246 * ui1 + 4642 * ui2) +
        ui1 * (7043 * ui1 - 3882 * ui2) + ui2 * (547 * ui2);

    ReconstructedValues betaR0 =
        ui_3 * (547 * ui_3 - 3882 * ui_2 + 4642 * ui_1 - 1854 * ui) +
        ui_2 * (7043 * ui_2 - 17246 * ui_1 + 7042 * ui) +
        ui_1 * (11003 * ui_1 - 9402 * ui) + ui * (2107 * ui);
    ReconstructedValues betaR1 =
        ui_2 * (267 * ui_2 - 1642 * ui_1 + 1602 * ui - 494 * ui1) +
        ui_1 * (2843 * ui_1 - 5966 * ui + 1922 * ui1) +
        ui * (3443 * ui - 2522 * ui1) + ui1 * (547 * ui1);
    ReconstructedValues betaR2 =
        ui_1 * (547 * ui_1 - 2522 * ui + 1922 * ui1 - 494 * ui2) +
        ui * (3443 * ui - 5966 * ui1 + 1602 * ui2) +
        ui1 * (2843 * ui1 - 1642 * ui2) + ui2 * (267 * ui2);
    ReconstructedValues betaR3 =
        ui * (2107 * ui - 9402 * ui1 + 7042 * ui2 - 1854 * ui3) +
        ui1 * (11003 * ui1 - 17246 * ui2 + 4642 * ui3) +
        ui2 * (7043 * ui2 - 3882 * ui3) + ui3 * (547 * ui3);

    ReconstructedValues alphaL0 = dL0 / ((betaL0 + eps) * (betaL0 + eps));
    ReconstructedValues alphaL1 = dL1 / ((betaL1 + eps) * (betaL1 + eps));
    ReconstructedValues alphaL2 = dL2 / ((betaL2 + eps) * (betaL2 + eps));
    ReconstructedValues alphaL3 = dL3 / ((betaL3 + eps) * (betaL3 + eps));

    ReconstructedValues alphaR0 = dR0 / ((betaR0 + eps) * (betaR0 + eps));
    ReconstructedValues alphaR1 = dR1 / ((betaR1 + eps) * (betaR1 + eps));
    ReconstructedValues alphaR2 = dR2 / ((betaR2 + eps) * (betaR2 + eps));
    ReconstructedValues alphaR3 = dR3 / ((betaR3 + eps) * (betaR3 + eps));

    ReconstructedValues wL0 = alphaL0 / (alphaL0 + alphaL1 + alphaL2 + alphaL3);
    ReconstructedValues wL1 = alphaL1 / (alphaL0 + alphaL1 + alphaL2 + alphaL3);
    ReconstructedValues wL2 = alphaL2 / (alphaL0 + alphaL1 + alphaL2 + alphaL3);
    ReconstructedValues wL3 = alphaL3 / (alphaL0 + alphaL1 + alphaL2 + alphaL3);

    ReconstructedValues wR0 = alphaR0 / (alphaR0 + alphaR1 + alphaR2 + alphaR3);
    ReconstructedValues wR1 = alphaR1 / (alphaR0 + alphaR1 + alphaR2 + alphaR3);
    ReconstructedValues wR2 = alphaR2 / (alphaR0 + alphaR1 + alphaR2 + alphaR3);
    ReconstructedValues wR3 = alphaR3 / (alphaR0 + alphaR1 + alphaR2 + alphaR3);

    uL = wL0 * (-(1. / 4.) * ui_4 + (13. / 12.) * ui_3 - (23. / 12.) * ui_2 +
                (25. / 12.) * ui_1) +
         wL1 * ((1. / 12.) * ui_3 - (5. / 12.) * ui_2 + (13. / 12.) * ui_1 +
                (1. / 4.) * ui) +
         wL2 * (-(1. / 12.) * ui_2 + (7. / 12.) * ui_1 + (7. / 12.) * ui -
                (1. / 12.) * ui1) +
         wL3 * ((1. / 4.) * ui_1 + (13. / 12.) * ui - (5. / 12.) * ui1 +
                (1. / 12.) * ui2);
    uR = wR0 * ((1. / 4.) * ui + (13. / 12.) * ui_1 - (5. / 12.) * ui_2 +
                (1. / 12.) * ui_3) +
         wR1 * (-(1. / 12.) * ui1 + (7. / 12.) * ui + (7. / 12.) * ui_1 -
                (1. / 12.) * ui_2) +
         wR2 * ((1. / 12.) * ui2 - (5. / 12.) * ui1 + (13. / 12.) * ui +
                (1. / 4.) * ui_1) +
         wR3 * (-(1. / 4.) * ui3 + (13. / 12.) * ui2 - (23. / 12.) * ui1 +
                (25. / 12.) * ui);
  }

  fL = ComputeFluxVector(uL, dir);
  fR = ComputeFluxVector(uR, dir);

  if (OptP == OptionalPhysics::HallResHyper) {
    auto [Lx, Rx] =
        ComputeRiemannJ(P_cc.Jx, AVERAGEY, i, j, Dx, Dy, rec, sl, dir);
    auto [JxL, LaplJxL] = Lx;
    auto [JxR, LaplJxR] = Rx;

    auto [Ly, Ry] =
        ComputeRiemannJ(P_cc.Jy, AVERAGEX, i, j, Dx, Dy, rec, sl, dir);
    auto [JyL, LaplJyL] = Ly;
    auto [JyR, LaplJyR] = Ry;

    auto [Lz, Rz] =
        ComputeRiemannJ(P_cc.Jz, AVERAGEXY, i, j, Dx, Dy, rec, sl, dir);
    auto [JzL, LaplJzL] = Lz;
    auto [JzR, LaplJzR] = Rz;

    AddNonIdealFlux(fL, uL, JxL, JyL, JzL, LaplJxL, LaplJyL, LaplJzL, OptP,
                    dir);
    AddNonIdealFlux(fR, uR, JxR, JyR, JzR, LaplJxR, LaplJyR, LaplJzR, OptP,
                    dir);
  }

  double c0L = std::sqrt((pc.gam * uL.P) / uL.rho); // Sound speeds
  double c0R = std::sqrt((pc.gam * uR.P) / uR.rho);

  double caxL = std::sqrt((uL.Bx * uL.Bx) / (uL.rho)); // Alfven speeds in x
  double caxR = std::sqrt((uR.Bx * uR.Bx) / (uR.rho));
  double cayL = std::sqrt((uL.By * uL.By) / (uL.rho)); // Alfven speeds in y
  double cayR = std::sqrt((uR.By * uR.By) / (uR.rho));

  double caL = std::sqrt((uL.Bx * uL.Bx + uL.By * uL.By + uL.Bz * uL.Bz) /
                         (uL.rho)); // Alfven speeds
  double caR =
      std::sqrt((uR.Bx * uR.Bx + uR.By * uR.By + uR.Bz * uR.Bz) / (uR.rho));

  cfastxL =
      std::sqrt((c0L * c0L + caL * caL) * 0.5 +
                (std::sqrt((c0L * c0L + caL * caL) * (c0L * c0L + caL * caL) -
                           4 * c0L * c0L * caxL * caxL)) *
                    0.5); // Fast magnetosonic speeds in x
  cfastxR =
      std::sqrt((c0R * c0R + caR * caR) * 0.5 +
                (std::sqrt((c0R * c0R + caR * caR) * (c0R * c0R + caR * caR) -
                           4 * c0R * c0R * caxR * caxR)) *
                    0.5);
  cfastyL =
      std::sqrt((c0L * c0L + caL * caL) * 0.5 +
                (std::sqrt((c0L * c0L + caL * caL) * (c0L * c0L + caL * caL) -
                           4 * c0L * c0L * cayL * cayL)) *
                    0.5); // Fast magnetosonic speeds in y
  cfastyR =
      std::sqrt((c0R * c0R + caR * caR) * 0.5 +
                (std::sqrt((c0R * c0R + caR * caR) * (c0R * c0R + caR * caR) -
                           4 * c0R * c0R * cayR * cayR)) *
                    0.5);

  if (OptP == OptionalPhysics::HallResHyper) {
    double vwx = (std::sqrt(1 + 0.25 / (Dx * Dx)) + 0.5 / Dx);
    double vwy = (std::sqrt(1 + 0.25 / (Dy * Dy)) + 0.5 / Dy);
    cwxL = std::sqrt(uL.Bx * uL.Bx + uL.By * uL.By + uL.Bz * uL.Bz) / (uL.rho) *
           vwx;
    cwyL = std::sqrt(uL.Bx * uL.Bx + uL.By * uL.By + uL.Bz * uL.Bz) / (uL.rho) *
           vwy;
    cwxR = std::sqrt(uR.Bx * uR.Bx + uR.By * uR.By + uR.Bz * uR.Bz) / (uR.rho) *
           vwx;
    cwyR = std::sqrt(uR.Bx * uR.Bx + uR.By * uR.By + uR.Bz * uR.Bz) / (uR.rho) *
           vwy;
    // double hallxL = std::sqrt(uL.Bx*uL.Bx + uL.By*uL.By + uL.Bz*uL.Bz) /
    // (2.0*(uL.rho)*Dx); double hallyL = std::sqrt(uL.Bx*uL.Bx + uL.By*uL.By +
    // uL.Bz*uL.Bz) / (2.0*(uL.rho)*Dy); double hallxR = std::sqrt(uR.Bx*uR.Bx +
    // uR.By*uR.By + uR.Bz*uR.Bz) / (2.0*(uR.rho)*Dx); double hallyR =
    // std::sqrt(uR.Bx*uR.Bx + uR.By*uR.By + uR.Bz*uR.Bz) / (2.0*(uR.rho)*Dy);
    // cwxL = std::fabs(hallxL) + std::sqrt(hallxL*hallxL + caxL*caxL);
    // cwyL = std::fabs(hallyL) + std::sqrt(hallyL*hallyL + cayL*cayL);
    // cwxR = std::fabs(hallxR) + std::sqrt(hallxR*hallxR + caxR*caxR);
    // cwyR = std::fabs(hallyR) + std::sqrt(hallyR*hallyR + cayR*cayR);

    if (dir == Dir::X) {
      SLb = std::min(uL.vx - cfastxL - cwxL, uR.vx - cfastxR - cwxR);
      SRb = std::max(uL.vx + cfastxL + cwxL, uR.vx + cfastxR + cwxR);
      // For rusanov :
      Splusb = std::max(std::abs(uL.vx) + cfastxL + cwxL,
                        std::abs(uR.vx) + cfastxR + cwxR);
    } else if (dir == Dir::Y) {
      SLb = std::min(uL.vy - cfastyL - cwyL, uR.vy - cfastyR - cwyR);
      SRb = std::max(uL.vy + cfastyL + cwyL, uR.vy + cfastyR + cwyR);
      //  For rusanov :
      Splusb = std::max(std::abs(uL.vy) + cfastyL + cwyL,
                        std::abs(uR.vy) + cfastyR + cwyR);
    }
  }

  // Wave speeds
  if (dir == Dir::X) {
    SL = std::min(uL.vx - cfastxL, uR.vx - cfastxR);
    SR = std::max(uL.vx + cfastxL, uR.vx + cfastxR);
    // For rusanov :
    Splus = std::max(std::abs(uL.vx) + cfastxL, std::abs(uR.vx) + cfastxR);
  } else if (dir == Dir::Y) {
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
