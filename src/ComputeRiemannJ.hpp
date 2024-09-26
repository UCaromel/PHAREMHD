#ifndef COMPUTE_RIEMANN_J_HPP_
#define COMPUTE_RIEMANN_J_HPP_

#include <cmath>

#include "ModularityUtils.hpp"

template <typename Func>
std::pair<std::pair<double, double>, std::pair<double, double>>
ComputeRiemannJ(const std::vector<std::vector<double>> &J, Func AVERAGING,
                int i, int j, double Dx, double Dy, Reconstruction rec,
                Slope sl, Dir Dir) {
  double JL;
  double JR;
  double LaplJL;
  double LaplJR;

  // J has one more ghost cell
  i++;
  j++;

  if (rec == Reconstruction::Constant) {
    // for constant reconstruction, the reconstructed values are the left and
    // right cell centered values next to the interface for x derivatives
    double J_2x = AVERAGING(J, i - 2, j);
    double J_1x = AVERAGING(J, i - 1, j);
    double Jx = AVERAGING(J, i, j);
    double J1x = AVERAGING(J, i + 1, j);

    // for y derivatives
    double J_2y = AVERAGING(J, i, j - 2);
    double J_1y = AVERAGING(J, i, j - 1);
    double Jy = AVERAGING(J, i, j);
    double J1y = AVERAGING(J, i, j + 1);

    LaplJL = (J_2x - 2.0 * J_1x + Jx) / (Dx * Dx) +
             (J_2y - 2.0 * J_1y + Jy) / (Dy * Dy);
    LaplJR = (J_1x - 2.0 * Jx + J1x) / (Dx * Dx) +
             (J_1y - 2.0 * Jy + J1y) / (Dy * Dy);

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
    // double J_3x = AVERAGING(J, i-3, j);
    double J_2x = AVERAGING(J, i - 2, j);
    double J_1x = AVERAGING(J, i - 1, j);
    double Jx = AVERAGING(J, i, j);
    double J1x = AVERAGING(J, i + 1, j);
    // double J2x = AVERAGING(J, i+2, j);

    // for y derivatives
    // double J_3y = AVERAGING(J, i, j-3);
    double J_2y = AVERAGING(J, i, j - 2);
    double J_1y = AVERAGING(J, i, j - 1);
    double Jy = AVERAGING(J, i, j);
    double J1y = AVERAGING(J, i, j + 1);
    // double J2y = AVERAGING(J, i, j+2);

    // Slopes
    // double Di_1x = ChosenSLDouble(J_1x - J_2x, J_2x - J_3x);
    double Dix = ChosenSLDouble(Jx - J_1x, J_1x - J_2x);
    double Di1x = ChosenSLDouble(J1x - Jx, Jx - J_1x);
    // double Di2x = ChosenSLDouble(J2x - J1x, J1x - Jx);

    // double Di_1y = ChosenSLDouble(J_1y - J_2y, J_2y - J_3y);
    double Diy = ChosenSLDouble(Jy - J_1y, J_1y - J_2y);
    double Di1y = ChosenSLDouble(J1y - Jy, Jy - J_1y);
    // double Di2y = ChosenSLDouble(J2y - J1y, J1y - Jy);

    // for L interface Laplacian
    // double J_1xR = J_2x + 0.5 * Di_1x;
    double JxR = J_1x + 0.5 * Dix;
    // double J1xR = Jx + 0.5 * Di1x;

    // double J_1yR = J_2y + 0.5 * Di_1y;
    double JyR = J_1y + 0.5 * Diy;
    // double J1yR = Jy + 0.5 * Di1y;

    // for R interface Laplacian
    // double JxL = J_1x - 0.5 * Dix;
    double J1xL = Jx - 0.5 * Di1x;
    // double J2xL = J1x - 0.5 * Di2x;

    // double JyL = J_1y - 0.5 * Diy;
    double J1yL = Jy - 0.5 * Di1y;
    // double J2yL = J1y - 0.5 * Di2y;

    LaplJL = 0.0; //(J_1xR - 2.0 * JxR + J1xR) / (Dx * Dx) + (J_1yR - 2.0 * JyR
                  //+ J1yR) / (Dy * Dy);
    LaplJR = 0.0; //(JxL - 2.0 * J1xL + J2xL) / (Dx * Dx) + (JyL - 2.0 * J1yL +
                  // J2yL) / (Dy * Dy);

    if (Dir == Dir::X) {
      JL = JxR;
      JR = J1xL;
    } else if (Dir == Dir::Y) {
      JL = JyR;
      JR = J1yL;
    }

  } else if (rec == Reconstruction::WENO3) {
    // for x derivatives
    double J_2x = AVERAGING(J, i - 2, j);
    double J_1x = AVERAGING(J, i - 1, j);
    double Jx = AVERAGING(J, i, j);
    double J1x = AVERAGING(J, i + 1, j);

    // for y derivatives
    double J_2y = AVERAGING(J, i, j - 2);
    double J_1y = AVERAGING(J, i, j - 1);
    double Jy = AVERAGING(J, i, j);
    double J1y = AVERAGING(J, i, j + 1);

    double eps = 1.e-6;

    // Calculate beta values for WENO
    double beta0_x = (J_1x - J_2x) * (J_1x - J_2x);
    double beta1_x = (Jx - J_1x) * (Jx - J_1x);
    double beta2_x = (J1x - Jx) * (J1x - Jx);

    double beta0_y = (J_1y - J_2y) * (J_1y - J_2y);
    double beta1_y = (Jy - J_1y) * (Jy - J_1y);
    double beta2_y = (J1y - Jy) * (J1y - Jy);

    // Calculate alpha values for WENO
    double alphaL0_x = (1. / 3.) / ((beta0_x + eps) * (beta0_x + eps));
    double alphaL1_x = (2. / 3.) / ((beta1_x + eps) * (beta1_x + eps));
    double alphaR0_x = (2. / 3.) / ((beta1_x + eps) * (beta1_x + eps));
    double alphaR1_x = (1. / 3.) / ((beta2_x + eps) * (beta2_x + eps));

    double alphaL0_y = (1. / 3.) / ((beta0_y + eps) * (beta0_y + eps));
    double alphaL1_y = (2. / 3.) / ((beta1_y + eps) * (beta1_y + eps));
    double alphaR0_y = (2. / 3.) / ((beta1_y + eps) * (beta1_y + eps));
    double alphaR1_y = (1. / 3.) / ((beta2_y + eps) * (beta2_y + eps));

    // Calculate weights for WENO
    double wL0_x = alphaL0_x / (alphaL0_x + alphaL1_x);
    double wL1_x = alphaL1_x / (alphaL0_x + alphaL1_x);
    double wR0_x = alphaR0_x / (alphaR0_x + alphaR1_x);
    double wR1_x = alphaR1_x / (alphaR0_x + alphaR1_x);

    double wL0_y = alphaL0_y / (alphaL0_y + alphaL1_y);
    double wL1_y = alphaL1_y / (alphaL0_y + alphaL1_y);
    double wR0_y = alphaR0_y / (alphaR0_y + alphaR1_y);
    double wR1_y = alphaR1_y / (alphaR0_y + alphaR1_y);

    // Calculate reconstructed values for left and right states
    double JL_x =
        wL0_x * (-0.5 * J_2x + 1.5 * J_1x) + wL1_x * (0.5 * J_1x + 0.5 * Jx);
    double JR_x =
        wR0_x * (0.5 * Jx + 0.5 * J_1x) + wR1_x * (-0.5 * J1x + 1.5 * Jx);

    double JL_y =
        wL0_y * (-0.5 * J_2y + 1.5 * J_1y) + wL1_y * (0.5 * J_1y + 0.5 * Jy);
    double JR_y =
        wR0_y * (0.5 * Jy + 0.5 * J_1y) + wR1_y * (-0.5 * J1y + 1.5 * Jy);

    LaplJL = 0.0;
    LaplJR = 0.0;

    if (Dir == Dir::X) {
      JL = JL_x;
      JR = JR_x;
    } else if (Dir == Dir::Y) {
      JL = JL_y;
      JR = JR_y;
    }
  } else if (rec == Reconstruction::WENO5) {
    double J_3x = AVERAGING(J, i - 3, j);
    double J_2x = AVERAGING(J, i - 2, j);
    double J_1x = AVERAGING(J, i - 1, j);
    double Jx = AVERAGING(J, i, j);
    double J1x = AVERAGING(J, i + 1, j);
    double J2x = AVERAGING(J, i + 2, j);

    double J_3y = AVERAGING(J, i, j - 3);
    double J_2y = AVERAGING(J, i, j - 2);
    double J_1y = AVERAGING(J, i, j - 1);
    double Jy = AVERAGING(J, i, j);
    double J1y = AVERAGING(J, i, j + 1);
    double J2y = AVERAGING(J, i, j + 2);

    double eps = 1.e-6;

    double dL0 = 1. / 10.;
    double dL1 = 3. / 5.;
    double dL2 = 3. / 10.;
    double dR0 = 3. / 10.;
    double dR1 = 3. / 5.;
    double dR2 = 1. / 10.;

    double betaL0_x =
        (13. / 12.) * (J_3x - 2 * J_2x + J_1x) * (J_3x - 2 * J_2x + J_1x) +
        (1. / 4.) * (J_3x - 4 * J_2x + 3 * J_1x) * (J_3x - 4 * J_2x + 3 * J_1x);

    double betaL1_x =
        (13. / 12.) * (J_2x - 2 * J_1x + Jx) * (J_2x - 2 * J_1x + Jx) +
        (1. / 4.) * (J_2x - Jx) * (J_2x - Jx);

    double betaL2_x =
        (13. / 12.) * (J_1x - 2 * Jx + J1x) * (J_1x - 2 * Jx + J1x) +
        (1. / 4.) * (3 * J_1x - 4 * Jx + J1x) * (3 * J_1x - 4 * Jx + J1x);

    double betaL0_y =
        (13. / 12.) * (J_3y - 2 * J_2y + J_1y) * (J_3y - 2 * J_2y + J_1y) +
        (1. / 4.) * (J_3y - 4 * J_2y + 3 * J_1y) * (J_3y - 4 * J_2y + 3 * J_1y);

    double betaL1_y =
        (13. / 12.) * (J_2y - 2 * J_1y + Jy) * (J_2y - 2 * J_1y + Jy) +
        (1. / 4.) * (J_2y - Jy) * (J_2y - Jy);

    double betaL2_y =
        (13. / 12.) * (J_1y - 2 * Jy + J1y) * (J_1y - 2 * Jy + J1y) +
        (1. / 4.) * (3 * J_1y - 4 * Jy + J1y) * (3 * J_1y - 4 * Jy + J1y);

    double betaR0_x =
        (13. / 12.) * (J_2x - 2 * J_1x + Jx) * (J_2x - 2 * J_1x + Jx) +
        (1. / 4.) * (J_2x - 4 * J_1x + 3 * Jx) * (J_2x - 4 * J_1x + 3 * Jx);

    double betaR1_x =
        (13. / 12.) * (J_1x - 2 * Jx + J1x) * (J_1x - 2 * Jx + J1x) +
        (1. / 4.) * (J_1x - J1x) * (J_1x - J1x);

    double betaR2_x =
        (13. / 12.) * (Jx - 2 * J1x + J2x) * (Jx - 2 * J1x + J2x) +
        (1. / 4.) * (3 * Jx - 4 * J1x + J2x) * (3 * Jx - 4 * J1x + J2x);

    double betaR0_y =
        (13. / 12.) * (J_2y - 2 * J_1y + Jy) * (J_2y - 2 * J_1y + Jy) +
        (1. / 4.) * (J_2y - 4 * J_1y + 3 * Jy) * (J_2y - 4 * J_1y + 3 * Jy);

    double betaR1_y =
        (13. / 12.) * (J_1y - 2 * Jy + J1y) * (J_1y - 2 * Jy + J1y) +
        (1. / 4.) * (J_1y - J1y) * (J_1y - J1y);

    double betaR2_y =
        (13. / 12.) * (Jy - 2 * J1y + J2y) * (Jy - 2 * J1y + J2y) +
        (1. / 4.) * (3 * Jy - 4 * J1y + J2y) * (3 * Jy - 4 * J1y + J2y);

    double alphaL0_x = dL0 / ((betaL0_x + eps) * (betaL0_x + eps));
    double alphaL1_x = dL1 / ((betaL1_x + eps) * (betaL1_x + eps));
    double alphaL2_x = dL2 / ((betaL2_x + eps) * (betaL2_x + eps));

    double alphaL0_y = dL0 / ((betaL0_y + eps) * (betaL0_y + eps));
    double alphaL1_y = dL1 / ((betaL1_y + eps) * (betaL1_y + eps));
    double alphaL2_y = dL2 / ((betaL2_y + eps) * (betaL2_y + eps));

    double alphaR0_x = dR0 / ((betaR0_x + eps) * (betaR0_x + eps));
    double alphaR1_x = dR1 / ((betaR1_x + eps) * (betaR1_x + eps));
    double alphaR2_x = dR2 / ((betaR2_x + eps) * (betaR2_x + eps));

    double alphaR0_y = dR0 / ((betaR0_y + eps) * (betaR0_y + eps));
    double alphaR1_y = dR1 / ((betaR1_y + eps) * (betaR1_y + eps));
    double alphaR2_y = dR2 / ((betaR2_y + eps) * (betaR2_y + eps));

    double wL0_x = alphaL0_x / (alphaL0_x + alphaL1_x + alphaL2_x);
    double wL1_x = alphaL1_x / (alphaL0_x + alphaL1_x + alphaL2_x);
    double wL2_x = alphaL2_x / (alphaL0_x + alphaL1_x + alphaL2_x);

    double wR0_x = alphaR0_x / (alphaR0_x + alphaR1_x + alphaR2_x);
    double wR1_x = alphaR1_x / (alphaR0_x + alphaR1_x + alphaR2_x);
    double wR2_x = alphaR2_x / (alphaR0_x + alphaR1_x + alphaR2_x);

    double wL0_y = alphaL0_y / (alphaL0_y + alphaL1_y + alphaL2_y);
    double wL1_y = alphaL1_y / (alphaL0_y + alphaL1_y + alphaL2_y);
    double wL2_y = alphaL2_y / (alphaL0_y + alphaL1_y + alphaL2_y);

    double wR0_y = alphaR0_y / (alphaR0_y + alphaR1_y + alphaR2_y);
    double wR1_y = alphaR1_y / (alphaR0_y + alphaR1_y + alphaR2_y);
    double wR2_y = alphaR2_y / (alphaR0_y + alphaR1_y + alphaR2_y);

    double JL_x =
        wL0_x * ((1. / 3.) * J_3x - (7. / 6.) * J_2x + (11. / 6.) * J_1x) +
        wL1_x * (-(1. / 6.) * J_2x + (5. / 6.) * J_1x + (1. / 3.) * Jx) +
        wL2_x * ((1. / 3.) * J_1x + (5. / 6.) * Jx - (1. / 6.) * J1x);

    double JR_x =
        wR0_x * ((1. / 3.) * Jx + (5. / 6.) * J_1x - (1. / 6.) * J_2x) +
        wR1_x * (-(1. / 6.) * J1x + (5. / 6.) * Jx + (1. / 3.) * J_1x) +
        wR2_x * ((1. / 3.) * J_2x - (7. / 6.) * J1x + (11. / 6.) * Jx);

    double JL_y =
        wL0_y * ((1. / 3.) * J_3y - (7. / 6.) * J_2y + (11. / 6.) * J_1y) +
        wL1_y * (-(1. / 6.) * J_2y + (5. / 6.) * J_1y + (1. / 3.) * Jy) +
        wL2_y * ((1. / 3.) * J_1y + (5. / 6.) * Jy - (1. / 6.) * J1y);

    double JR_y =
        wR0_y * ((1. / 3.) * Jy + (5. / 6.) * J_1y - (1. / 6.) * J_2y) +
        wR1_y * (-(1. / 6.) * J1y + (5. / 6.) * Jy + (1. / 3.) * J_1y) +
        wR2_y * ((1. / 3.) * J_2y - (7. / 6.) * J1y + (11. / 6.) * Jy);

    LaplJL = 0.0;
    LaplJR = 0.0;

    if (Dir == Dir::X) {
      JL = JL_x;
      JR = JR_x;
    } else if (Dir == Dir::Y) {
      JL = JL_y;
      JR = JR_y;
    }
  } else if (rec == Reconstruction::WENO7) {
    double J_4x = AVERAGING(J, i - 4, j);
    double J_3x = AVERAGING(J, i - 3, j);
    double J_2x = AVERAGING(J, i - 2, j);
    double J_1x = AVERAGING(J, i - 1, j);
    double Jx = AVERAGING(J, i, j);
    double J1x = AVERAGING(J, i + 1, j);
    double J2x = AVERAGING(J, i + 2, j);
    double J3x = AVERAGING(J, i + 3, j);

    double J_4y = AVERAGING(J, i, j - 4);
    double J_3y = AVERAGING(J, i, j - 3);
    double J_2y = AVERAGING(J, i, j - 2);
    double J_1y = AVERAGING(J, i, j - 1);
    double Jy = AVERAGING(J, i, j);
    double J1y = AVERAGING(J, i, j + 1);
    double J2y = AVERAGING(J, i, j + 2);
    double J3y = AVERAGING(J, i, j + 3);

    double eps = 1.e-6;

    double dL0 = 1. / 35.;
    double dL1 = 12. / 35.;
    double dL2 = 18. / 35.;
    double dL3 = 4. / 35.;
    double dR0 = 4. / 35.;
    double dR1 = 18. / 35.;
    double dR2 = 12. / 35.;
    double dR3 = 1. / 35.;

    // Beta coefficients for left side (x-direction)
    double betaL0_x =
        J_4x * (547 * J_4x - 3882 * J_3x + 4642 * J_2x - 1854 * J_1x) +
        J_3x * (7043 * J_3x - 17246 * J_2x + 7042 * J_1x) +
        J_2x * (11003 * J_2x - 9402 * J_1x) + J_1x * (2107 * J_1x);
    double betaL1_x =
        J_3x * (267 * J_3x - 1642 * J_2x + 1602 * J_1x - 494 * Jx) +
        J_2x * (2843 * J_2x - 5966 * J_1x + 1922 * Jx) +
        J_1x * (3443 * J_1x - 2522 * Jx) + Jx * (547 * Jx);
    double betaL2_x =
        J_2x * (547 * J_2x - 2522 * J_1x + 1922 * Jx - 494 * J1x) +
        J_1x * (3443 * J_1x - 5966 * Jx + 1602 * J1x) +
        Jx * (2843 * Jx - 1642 * J1x) + J1x * (267 * J1x);
    double betaL3_x =
        J_1x * (2107 * J_1x - 9402 * Jx + 7042 * J1x - 1854 * J2x) +
        Jx * (11003 * Jx - 17246 * J1x + 4642 * J2x) +
        J1x * (7043 * J1x - 3882 * J2x) + J2x * (547 * J2x);

    // Beta coefficients for right side (x-direction)
    double betaR0_x =
        J_3x * (547 * J_3x - 3882 * J_2x + 4642 * J_1x - 1854 * Jx) +
        J_2x * (7043 * J_2x - 17246 * J_1x + 7042 * Jx) +
        J_1x * (11003 * J_1x - 9402 * Jx) + Jx * (2107 * Jx);
    double betaR1_x =
        J_2x * (267 * J_2x - 1642 * J_1x + 1602 * Jx - 494 * J1x) +
        J_1x * (2843 * J_1x - 5966 * Jx + 1922 * J1x) +
        Jx * (3443 * Jx - 2522 * J1x) + J1x * (547 * J1x);
    double betaR2_x = J_1x * (547 * J_1x - 2522 * Jx + 1922 * J1x - 494 * J2x) +
                      Jx * (3443 * Jx - 5966 * J1x + 1602 * J2x) +
                      J1x * (2843 * J1x - 1642 * J2x) + J2x * (267 * J2x);
    double betaR3_x = Jx * (2107 * Jx - 9402 * J1x + 7042 * J2x - 1854 * J3x) +
                      J1x * (11003 * J1x - 17246 * J2x + 4642 * J3x) +
                      J2x * (7043 * J2x - 3882 * J3x) + J3x * (547 * J3x);

    // Beta coefficients for left side (y-direction)
    double betaL0_y =
        J_4y * (547 * J_4y - 3882 * J_3y + 4642 * J_2y - 1854 * J_1y) +
        J_3y * (7043 * J_3y - 17246 * J_2y + 7042 * J_1y) +
        J_2y * (11003 * J_2y - 9402 * J_1y) + J_1y * (2107 * J_1y);
    double betaL1_y =
        J_3y * (267 * J_3y - 1642 * J_2y + 1602 * J_1y - 494 * Jy) +
        J_2y * (2843 * J_2y - 5966 * J_1y + 1922 * Jy) +
        J_1y * (3443 * J_1y - 2522 * Jy) + Jy * (547 * Jy);
    double betaL2_y =
        J_2y * (547 * J_2y - 2522 * J_1y + 1922 * Jy - 494 * J1y) +
        J_1y * (3443 * J_1y - 5966 * Jy + 1602 * J1y) +
        Jy * (2843 * Jy - 1642 * J1y) + J1y * (267 * J1y);
    double betaL3_y =
        J_1y * (2107 * J_1y - 9402 * Jy + 7042 * J1y - 1854 * J2y) +
        Jy * (11003 * Jy - 17246 * J1y + 4642 * J2y) +
        J1y * (7043 * J1y - 3882 * J2y) + J2y * (547 * J2y);

    // Beta coefficients for right side (y-direction)
    double betaR0_y =
        J_3y * (547 * J_3y - 3882 * J_2y + 4642 * J_1y - 1854 * Jy) +
        J_2y * (7043 * J_2y - 17246 * J_1y + 7042 * Jy) +
        J_1y * (11003 * J_1y - 9402 * Jy) + Jy * (2107 * Jy);
    double betaR1_y =
        J_2y * (267 * J_2y - 1642 * J_1y + 1602 * Jy - 494 * J1y) +
        J_1y * (2843 * J_1y - 5966 * Jy + 1922 * J1y) +
        Jy * (3443 * Jy - 2522 * J1y) + J1y * (547 * J1y);
    double betaR2_y = J_1y * (547 * J_1y - 2522 * Jy + 1922 * J1y - 494 * J2y) +
                      Jy * (3443 * Jy - 5966 * J1y + 1602 * J2y) +
                      J1y * (2843 * J1y - 1642 * J2y) + J2y * (267 * J2y);
    double betaR3_y = Jy * (2107 * Jy - 9402 * J1y + 7042 * J2y - 1854 * J3y) +
                      J1y * (11003 * J1y - 17246 * J2y + 4642 * J3y) +
                      J2y * (7043 * J2y - 3882 * J3y) + J3y * (547 * J3y);

    // Alpha coefficients for left and right sides (x-direction)
    double alphaL0_x = dL0 / ((betaL0_x + eps) * (betaL0_x + eps));
    double alphaL1_x = dL1 / ((betaL1_x + eps) * (betaL1_x + eps));
    double alphaL2_x = dL2 / ((betaL2_x + eps) * (betaL2_x + eps));
    double alphaL3_x = dL3 / ((betaL3_x + eps) * (betaL3_x + eps));

    double alphaR0_x = dR0 / ((betaR0_x + eps) * (betaR0_x + eps));
    double alphaR1_x = dR1 / ((betaR1_x + eps) * (betaR1_x + eps));
    double alphaR2_x = dR2 / ((betaR2_x + eps) * (betaR2_x + eps));
    double alphaR3_x = dR3 / ((betaR3_x + eps) * (betaR3_x + eps));

    // Alpha coefficients for left and right sides (y-direction)
    double alphaL0_y = dL0 / ((betaL0_y + eps) * (betaL0_y + eps));
    double alphaL1_y = dL1 / ((betaL1_y + eps) * (betaL1_y + eps));
    double alphaL2_y = dL2 / ((betaL2_y + eps) * (betaL2_y + eps));
    double alphaL3_y = dL3 / ((betaL3_y + eps) * (betaL3_y + eps));

    double alphaR0_y = dR0 / ((betaR0_y + eps) * (betaR0_y + eps));
    double alphaR1_y = dR1 / ((betaR1_y + eps) * (betaR1_y + eps));
    double alphaR2_y = dR2 / ((betaR2_y + eps) * (betaR2_y + eps));
    double alphaR3_y = dR3 / ((betaR3_y + eps) * (betaR3_y + eps));

    // Weight coefficients for left and right sides (x-direction)
    double wL0_x = alphaL0_x / (alphaL0_x + alphaL1_x + alphaL2_x + alphaL3_x);
    double wL1_x = alphaL1_x / (alphaL0_x + alphaL1_x + alphaL2_x + alphaL3_x);
    double wL2_x = alphaL2_x / (alphaL0_x + alphaL1_x + alphaL2_x + alphaL3_x);
    double wL3_x = alphaL3_x / (alphaL0_x + alphaL1_x + alphaL2_x + alphaL3_x);

    double wR0_x = alphaR0_x / (alphaR0_x + alphaR1_x + alphaR2_x + alphaR3_x);
    double wR1_x = alphaR1_x / (alphaR0_x + alphaR1_x + alphaR2_x + alphaR3_x);
    double wR2_x = alphaR2_x / (alphaR0_x + alphaR1_x + alphaR2_x + alphaR3_x);
    double wR3_x = alphaR3_x / (alphaR0_x + alphaR1_x + alphaR2_x + alphaR3_x);

    // Weight coefficients for left and right sides (y-direction)
    double wL0_y = alphaL0_y / (alphaL0_y + alphaL1_y + alphaL2_y + alphaL3_y);
    double wL1_y = alphaL1_y / (alphaL0_y + alphaL1_y + alphaL2_y + alphaL3_y);
    double wL2_y = alphaL2_y / (alphaL0_y + alphaL1_y + alphaL2_y + alphaL3_y);
    double wL3_y = alphaL3_y / (alphaL0_y + alphaL1_y + alphaL2_y + alphaL3_y);

    double wR0_y = alphaR0_y / (alphaR0_y + alphaR1_y + alphaR2_y + alphaR3_y);
    double wR1_y = alphaR1_y / (alphaR0_y + alphaR1_y + alphaR2_y + alphaR3_y);
    double wR2_y = alphaR2_y / (alphaR0_y + alphaR1_y + alphaR2_y + alphaR3_y);
    double wR3_y = alphaR3_y / (alphaR0_y + alphaR1_y + alphaR2_y + alphaR3_y);

    // Reconstructed values for left and right sides (x-direction)
    double JL_x = wL0_x * (-(1. / 4.) * J_4x + (13. / 12.) * J_3x -
                           (23. / 12.) * J_2x + (25. / 12.) * J_1x) +
                  wL1_x * ((1. / 12.) * J_3x - (5. / 12.) * J_2x +
                           (13. / 12.) * J_1x + (1. / 4.) * Jx) +
                  wL2_x * (-(1. / 12.) * J_2x + (7. / 12.) * J_1x +
                           (7. / 12.) * Jx - (1. / 12.) * J1x) +
                  wL3_x * ((1. / 4.) * J_1x + (13. / 12.) * Jx -
                           (5. / 12.) * J1x + (1. / 12.) * J2x);

    double JR_x = wR0_x * ((1. / 4.) * Jx + (13. / 12.) * J_1x -
                           (5. / 12.) * J_2x + (1. / 12.) * J_3x) +
                  wR1_x * (-(1. / 12.) * J1x + (7. / 12.) * Jx +
                           (7. / 12.) * J_1x - (1. / 12.) * J_2x) +
                  wR2_x * ((1. / 12.) * J2x - (5. / 12.) * J1x +
                           (13. / 12.) * Jx + (1. / 4.) * J_1x) +
                  wR3_x * (-(1. / 4.) * J3x + (13. / 12.) * J2x -
                           (23. / 12.) * J1x + (25. / 12.) * Jx);

    // Reconstructed values for left and right sides (y-direction)
    double JL_y = wL0_y * (-(1. / 4.) * J_4y + (13. / 12.) * J_3y -
                           (23. / 12.) * J_2y + (25. / 12.) * J_1y) +
                  wL1_y * ((1. / 12.) * J_3y - (5. / 12.) * J_2y +
                           (13. / 12.) * J_1y + (1. / 4.) * Jy) +
                  wL2_y * (-(1. / 12.) * J_2y + (7. / 12.) * J_1y +
                           (7. / 12.) * Jy - (1. / 12.) * J1y) +
                  wL3_y * ((1. / 4.) * J_1y + (13. / 12.) * Jy -
                           (5. / 12.) * J1y + (1. / 12.) * J2y);

    double JR_y = wR0_y * ((1. / 4.) * Jy + (13. / 12.) * J_1y -
                           (5. / 12.) * J_2y + (1. / 12.) * J_3y) +
                  wR1_y * (-(1. / 12.) * J1y + (7. / 12.) * Jy +
                           (7. / 12.) * J_1y - (1. / 12.) * J_2y) +
                  wR2_y * ((1. / 12.) * J2y - (5. / 12.) * J1y +
                           (13. / 12.) * Jy + (1. / 4.) * J_1y) +
                  wR3_y * (-(1. / 4.) * J3y + (13. / 12.) * J2y -
                           (23. / 12.) * J1y + (25. / 12.) * Jy);

    if (Dir == Dir::X) {
      JL = JL_x;
      JR = JR_x;
    } else if (Dir == Dir::Y) {
      JL = JL_y;
      JR = JR_y;
    }

    // Set LaplJL and LaplJR to 0.0 as in the WENO5 implementation
    LaplJL = 0.0;
    LaplJR = 0.0;
  } else {
    throw std::invalid_argument("Invalid reconstruction type");
  }

  return std::make_pair(std::make_pair(JL, LaplJL), std::make_pair(JR, LaplJR));
}

#endif
