#include "ConstrainedTransport.hpp"
#include "ComputeRiemannJ.hpp"
#include "ReconstructedValues.hpp"

static const PhysicalConstants &pc = PhysicalConstants::getInstance();

std::vector<std::vector<double>>
ConstrainedTransportAverage(const ConservativeVariables &Cn, double Dx,
                            double Dy, int nghost, Reconstruction rec, Slope sl,
                            Riemann rs, OptionalPhysics OptP) {
  // Edge-centered (reconstructions)
  std::vector<std::vector<double>> vx(
      Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
  std::vector<std::vector<double>> vy(
      Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
  std::vector<std::vector<double>> Bx(
      Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
  std::vector<std::vector<double>> By(
      Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));

  std::vector<std::vector<double>> Ez(
      Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));

  // Cell-centered
  std::vector<std::vector<double>> vxcc(Cn.ny, std::vector<double>(Cn.nx));
  std::vector<std::vector<double>> vycc(Cn.ny, std::vector<double>(Cn.nx));

  for (size_t i = 0; i < Cn.rho.size(); ++i) {
    for (size_t j = 0; j < Cn.rho[i].size(); ++j) {
      if (Cn.rho[i][j] == 0) {
        throw std::runtime_error("rho is zero");
      }
      vxcc[i][j] = Cn.rhovx[i][j] / Cn.rho[i][j];
      vycc[i][j] = Cn.rhovy[i][j] / Cn.rho[i][j];
    }
  }

  for (int j = nghost; j <= Cn.ny - nghost; ++j) {
    for (int i = nghost; i <= Cn.nx - nghost; ++i) {
      vx[j - nghost][i - nghost] = AVERAGEXY(vxcc, i - 1, j - 1);
      vy[j - nghost][i - nghost] = AVERAGEXY(vycc, i - 1, j - 1);
      Bx[j - nghost][i - nghost] = AVERAGEY(Cn.Bxf, i, j - 1);
      By[j - nghost][i - nghost] = AVERAGEX(Cn.Byf, i - 1, j);

      Ez[j - nghost][i - nghost] =
          vy[j - nghost][i - nghost] * Bx[j - nghost][i - nghost] -
          vx[j - nghost][i - nghost] * By[j - nghost][i - nghost];
    }
  }

  if (OptP == OptionalPhysics::HallResHyper) {
    // Edge-centered (reconstructions)
    std::vector<std::vector<double>> rho(
        Cn.ny + 1 - 2 * nghost,
        std::vector<double>(Cn.nx + 1 -
                            2 * nghost)); // used in hall effect. Since we only
                                          // have one mass, m=1 -> n = rho
    std::vector<std::vector<double>> jx(
        Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));
    std::vector<std::vector<double>> jy(
        Cn.ny + 1 - 2 * nghost, std::vector<double>(Cn.nx + 1 - 2 * nghost));

    for (int j = nghost; j <= Cn.ny - nghost; ++j) {
      for (int i = nghost; i <= Cn.nx - nghost; ++i) {
        int i1 = i + 1;
        int j1 = j + 1;

        rho[j - nghost][i - nghost] = AVERAGEXY(Cn.rho, i - 1, j - 1);
        jx[j - nghost][i - nghost] = AVERAGEX(Cn.Jx, i1 - 1, j1);
        jy[j - nghost][i - nghost] = AVERAGEY(Cn.Jy, i1, j1 - 1);

        Ez[j - nghost][i - nghost] +=
            (1.0 / rho[j - nghost][i - nghost]) *
            (jx[j - nghost][i - nghost] * By[j - nghost][i - nghost] -
             jy[j - nghost][i - nghost] *
                 Bx[j - nghost][i - nghost]); // Hall contribution
        // Ez[j - nghost][i - nghost] += pc.eta * (Cn.Jz[j1][i1]); //
        // Resistivity Ez[j - nghost][i - nghost] -= pc.nu * LAPLACIAN(Cn.Jz,
        // i1, j1, Dx, Dy); // Hyper resistivity
      }
    }
  }

  return Ez;
}

std::vector<std::vector<double>>
ConstrainedTransportArithmetic(const ConservativeVariables &Cn, double Dx,
                               double Dy, int nghost, Reconstruction rec,
                               Slope sl, Riemann rs, OptionalPhysics OptP) {
  PrimitiveVariables Pn(Cn);
  std::vector<std::vector<Interface>> InterfacesX(
      Pn.ny + 2 - 2 * nghost, std::vector<Interface>(Pn.nx + 1 - 2 * nghost));
  std::vector<std::vector<Interface>> InterfacesY(
      Pn.ny + 1 - 2 * nghost, std::vector<Interface>(Pn.nx + 2 - 2 * nghost));

  std::vector<std::vector<double>> Ez(
      Pn.ny + 1 - 2 * nghost, std::vector<double>(Pn.nx + 1 - 2 * nghost));

  RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);

  for (int j = nghost - 1; j < Cn.ny - nghost + 1; ++j) {
    for (int i = nghost; i <= Cn.nx - nghost; ++i) {
      InterfacesX[j - nghost + 1][i - nghost] =
          Interface(Pn, i, j, Dx, Dy, nghost, rec, sl, OptP, Dir::X);
    }
  }

  for (int j = nghost; j <= Cn.ny - nghost; ++j) {
    for (int i = nghost - 1; i < Cn.nx - nghost + 1; ++i) {
      InterfacesY[j - nghost][i - nghost + 1] =
          Interface(Pn, i, j, Dx, Dy, nghost, rec, sl, OptP, Dir::Y);
    }
  }

  for (int j = 0; j <= Cn.ny - 2 * nghost; ++j) {
    for (int i = 0; i <= Cn.nx - 2 * nghost; ++i) {
      Interface x = InterfacesX[j][i];
      Interface x1 = InterfacesX[j + 1][i];
      Interface y = InterfacesY[j][i];
      Interface y1 = InterfacesY[j][i + 1];

      double ezi = -ChosenRiemannSolver(x).By;
      double ezi1 = -ChosenRiemannSolver(x1).By;
      double ezj = ChosenRiemannSolver(y).Bx;
      double ezj1 = ChosenRiemannSolver(y1).Bx;
      Ez[j][i] = 0.25 * (ezi + ezi1 + ezj + ezj1);
    }
  }
  return Ez;
}

std::vector<std::vector<double>>
ConstrainedTransportContact(const ConservativeVariables &Cn, double Dx,
                            double Dy, int nghost, Reconstruction rec, Slope sl,
                            Riemann rs, OptionalPhysics OptP) {
  PrimitiveVariables Pn(Cn);
  std::vector<std::vector<Interface>> InterfacesX(
      Pn.ny + 2 - 2 * nghost, std::vector<Interface>(Pn.nx + 1 - 2 * nghost));
  std::vector<std::vector<Interface>> InterfacesY(
      Pn.ny + 1 - 2 * nghost, std::vector<Interface>(Pn.nx + 2 - 2 * nghost));

  std::vector<std::vector<double>> Ez(
      Pn.ny + 1 - 2 * nghost, std::vector<double>(Pn.nx + 1 - 2 * nghost));

  RiemannSolverFunction ChosenRiemannSolver = getRiemannSolver(rs);

  for (int j = nghost - 1; j < Cn.ny - nghost + 1; ++j) {
    for (int i = nghost; i <= Cn.nx - nghost; ++i) {
      InterfacesX[j - nghost + 1][i - nghost] =
          Interface(Pn, i, j, Dx, Dy, nghost, rec, sl, OptP, Dir::X);
    }
  }

  for (int j = nghost; j <= Cn.ny - nghost; ++j) {
    for (int i = nghost - 1; i < Cn.nx - nghost + 1; ++i) {
      InterfacesY[j - nghost][i - nghost + 1] =
          Interface(Pn, i, j, Dx, Dy, nghost, rec, sl, OptP, Dir::Y);
    }
  }

  for (int j = 0; j <= Cn.ny - 2 * nghost; ++j) {
    for (int i = 0; i <= Cn.nx - 2 * nghost; ++i) {
      Interface x = InterfacesX[j][i];
      Interface x1 = InterfacesX[j + 1][i];
      Interface y = InterfacesY[j][i];
      Interface y1 = InterfacesY[j][i + 1];

      // Compute contact sign
      double Frhox = ChosenRiemannSolver(x).rho;
      double Frhox1 = ChosenRiemannSolver(x1).rho;
      double Frhoy = ChosenRiemannSolver(y).rho;
      double Frhoy1 = ChosenRiemannSolver(y1).rho;

      double sx = Frhox > 0 ? 1 : Frhox < 0 ? -1 : 0;
      double sx1 = Frhox1 > 0 ? 1 : Frhox1 < 0 ? -1 : 0;
      double sy = Frhoy > 0 ? 1 : Frhoy < 0 ? -1 : 0;
      double sy1 = Frhoy1 > 0 ? 1 : Frhoy1 < 0 ? -1 : 0;

      // Cell-centered EMF
      int offset = nghost - 1;

      ReconstructedValues Ux = Pn(i + offset, j + offset);
      ReconstructedValues Ux1 = Pn(i + offset, j + 1 + offset);
      ReconstructedValues Uy = Pn(i + offset, j + offset);
      ReconstructedValues Uy1 = Pn(i + 1 + offset, j + offset);
      ReconstructedValues Ux1y1 = Pn(i + 1 + offset, j + 1 + offset);
      ReconstructedValues Uy1x1 = Pn(i + 1 + offset, j + 1 + offset);

      ReconstructedValues Fx = ComputeFluxVector(Ux, Dir::X);
      ReconstructedValues Fx1 = ComputeFluxVector(Ux1, Dir::X);
      ReconstructedValues Fy = ComputeFluxVector(Uy, Dir::Y);
      ReconstructedValues Fy1 = ComputeFluxVector(Uy1, Dir::Y);
      ReconstructedValues Fx1y1 = ComputeFluxVector(Ux1y1, Dir::X);
      ReconstructedValues Fy1x1 = ComputeFluxVector(Uy1x1, Dir::Y);

      if (OptP == OptionalPhysics::HallResHyper) {
        // Lambda function to compute and add non-ideal flux
        auto NonIdealFlux = [&](ReconstructedValues &F, ReconstructedValues &U,
                                int i_offset, int j_offset, Dir dir) {
          auto [Lx, Rx] = ComputeRiemannJ(Pn.Jx, AVERAGEY, i + i_offset,
                                          j + j_offset, Dx, Dy, rec, sl, dir);
          auto [JxL, LaplJxL] = Lx;
          auto [JxR, LaplJxR] = Rx;

          auto [Ly, Ry] = ComputeRiemannJ(Pn.Jy, AVERAGEX, i + i_offset,
                                          j + j_offset, Dx, Dy, rec, sl, dir);
          auto [JyL, LaplJyL] = Ly;
          auto [JyR, LaplJyR] = Ry;

          auto [Lz, Rz] = ComputeRiemannJ(Pn.Jz, AVERAGEXY, i + i_offset,
                                          j + j_offset, Dx, Dy, rec, sl, dir);
          auto [JzL, LaplJzL] = Lz;
          auto [JzR, LaplJzR] = Rz;

          AddNonIdealFlux(F, U, JxL, JyL, JzL, LaplJxL, LaplJyL, LaplJzL, OptP,
                          dir);
          AddNonIdealFlux(F, U, JxR, JyR, JzR, LaplJxR, LaplJyR, LaplJzR, OptP,
                          dir);
        };

        // Apply non-ideal flux calculations to each reconstructed value
        NonIdealFlux(Fx, Ux, offset, offset, Dir::X);
        NonIdealFlux(Fx1, Ux1, offset, 1 + offset, Dir::X);
        NonIdealFlux(Fy, Uy, offset, offset, Dir::Y);
        NonIdealFlux(Fy1, Uy1, 1 + offset, offset, Dir::Y);
        NonIdealFlux(Fx1y1, Ux1y1, 1 + offset, 1 + offset, Dir::X);
        NonIdealFlux(Fy1x1, Uy1x1, 1 + offset, 1 + offset, Dir::Y);
      }

      double Ezx = -Fx.By;
      double Ezx1 = -Fx1.By;
      double Ezy = Fy.Bx;
      double Ezy1 = Fy1.Bx;

      double Ezx1y1 = -Fx1y1.By;
      double Ezy1x1 = Fy1x1.Bx;

      // Face centered EMF
      double ezx = -ChosenRiemannSolver(x).By;
      double ezx1 = -ChosenRiemannSolver(x1).By;
      double ezy = ChosenRiemannSolver(y).Bx;
      double ezy1 = ChosenRiemannSolver(y1).Bx;

      // Compute Ez derivatives in each directions
      double dyEzS = (1.0 + sx) * (ezy - Ezy) + (1.0 - sx) * (ezy1 - Ezy1);
      double dyEzN = (1.0 + sx1) * (Ezy1 - ezy) + (1.0 - sx1) * (Ezy1x1 - ezy1);
      double dxEzW = (1.0 + sy) * (ezx - Ezx) + (1.0 - sy) * (ezx1 - Ezx1);
      double dxEzE = (1.0 + sy1) * (Ezx1 - ezx) + (1.0 - sy1) * (Ezx1y1 - ezx1);

      Ez[j][i] = 0.25 * (ezx + ezx1 + ezy + ezy1) + 0.125 * (dyEzS - dyEzN) +
                 0.125 * (dxEzW - dxEzE);
    }
  }
  return Ez;
}

std::vector<std::vector<double>> UCTHLL(const ConservativeVariables &Cn,
                                        double Dx, double Dy, int nghost,
                                        Reconstruction rec, Slope sl,
                                        Riemann rs, OptionalPhysics OptP) {
  PrimitiveVariables Pn(Cn);
  std::vector<std::vector<Interface>> InterfacesX(
      Pn.ny + 2 - 2 * nghost, std::vector<Interface>(Pn.nx + 1 - 2 * nghost));
  std::vector<std::vector<Interface>> InterfacesY(
      Pn.ny + 1 - 2 * nghost, std::vector<Interface>(Pn.nx + 2 - 2 * nghost));

  std::vector<std::vector<double>> Ez(
      Pn.ny + 1 - 2 * nghost, std::vector<double>(Pn.nx + 1 - 2 * nghost));

  for (int j = nghost - 1; j < Cn.ny - nghost + 1; ++j) {
    for (int i = nghost; i <= Cn.nx - nghost; ++i) {
      InterfacesX[j - nghost + 1][i - nghost] =
          Interface(Pn, i, j, Dx, Dy, nghost, rec, sl, OptP, Dir::X);
    }
  }

  for (int j = nghost; j <= Cn.ny - nghost; ++j) {
    for (int i = nghost - 1; i < Cn.nx - nghost + 1; ++i) {
      InterfacesY[j - nghost][i - nghost + 1] =
          Interface(Pn, i, j, Dx, Dy, nghost, rec, sl, OptP, Dir::Y);
    }
  }

  for (int j = 0; j <= Cn.ny - 2 * nghost; ++j) {
    for (int i = 0; i <= Cn.nx - 2 * nghost; ++i) {
      Interface x = InterfacesX[j][i];
      Interface x1 = InterfacesX[j + 1][i];
      Interface y = InterfacesY[j][i];
      Interface y1 = InterfacesY[j][i + 1];

      // Compute solver's coeficients
      double alphaxR = std::max(0.0, x.uR.vx / x.uR.rho + x.cfastxR);
      double alphaxL = -std::min(0.0, x.uL.vx / x.uL.rho - x.cfastxL);

      double alphax1R = std::max(0.0, x1.uR.vx / x1.uR.rho + x1.cfastxR);
      double alphax1L = -std::min(0.0, x1.uL.vx / x1.uL.rho - x1.cfastxL);

      double alphayR = std::max(0.0, y.uR.vy / y.uR.rho + y.cfastyR);
      double alphayL = -std::min(0.0, y.uL.vy / y.uL.rho - y.cfastyL);

      double alphay1R = std::max(0.0, y1.uR.vy / y1.uR.rho + y1.cfastyR);
      double alphay1L = -std::min(0.0, y1.uL.vy / y1.uL.rho - y1.cfastyL);

      double axR = alphaxL / (alphaxR + alphaxL);
      double axL = alphaxR / (alphaxR + alphaxL);

      double ax1R = alphax1L / (alphax1R + alphax1L);
      double ax1L = alphax1R / (alphax1R + alphax1L);

      double ayR = alphayL / (alphayR + alphayL);
      double ayL = alphayR / (alphayR + alphayL);

      double ay1R = alphay1L / (alphay1R + alphay1L);
      double ay1L = alphay1R / (alphay1R + alphay1L);

      double dxR = (alphaxR * alphaxL) / (alphaxR + alphaxL);
      double dxL = (alphaxR * alphaxL) / (alphaxR + alphaxL);

      double dx1R = (alphax1R * alphax1L) / (alphax1R + alphax1L);
      double dx1L = (alphax1R * alphax1L) / (alphax1R + alphax1L);

      double dyR = (alphayR * alphayL) / (alphayR + alphayL);
      double dyL = (alphayR * alphayL) / (alphayR + alphayL);

      double dy1R = (alphay1R * alphay1L) / (alphay1R + alphay1L);
      double dy1L = (alphay1R * alphay1L) / (alphay1R + alphay1L);

      // Averaging
      double aW = 0.5 * (axL + ax1L);
      double aE = 0.5 * (axR + ax1R);
      double aS = 0.5 * (ayL + ay1L);
      double aN = 0.5 * (ayR + ay1R);

      double dW = 0.5 * (dxL + dx1L);
      double dE = 0.5 * (dxR + dx1R);
      double dS = 0.5 * (dyL + dy1L);
      double dN = 0.5 * (dyR + dy1R);

      double vxW = 0.5 * (y.uL.vx / y.uL.rho + y.uR.vx / y.uR.rho);
      double vxE = 0.5 * (y1.uL.vx / y1.uL.rho + y1.uR.vx / y1.uR.rho);
      double vyS = 0.5 * (x.uL.vy / x.uL.rho + x.uR.vy / x.uR.rho);
      double vyN = 0.5 * (x1.uL.vy / x1.uL.rho + x1.uR.vy / x1.uR.rho);

      double ByW = Cn.Byf[j + nghost][i + nghost - 1];
      double ByE = Cn.Byf[j + nghost][i + nghost];
      double BxS = Cn.Bxf[j + nghost - 1][i + nghost];
      double BxN = Cn.Bxf[j + nghost][i + nghost];

      Ez[j][i] = -(aW * vxW * ByW + aE * vxE * ByE) +
                 (aS * vyS * BxS + aN * vyN * BxN) + (dE * ByE - dW * ByW) -
                 (dN * BxN - dS * BxS);
    }
  }

  return Ez;
}

void ApplyConstrainedTransport(ConservativeVariables &Cn1,
                               const ConservativeVariables &Cn, double Dx,
                               double Dy, double Dt, int nghost,
                               Reconstruction rec, Slope sl, Riemann rs,
                               CTMethod ct, OptionalPhysics OptP) {
  CTFunction ChosenCT = getCT(ct);
  std::vector<std::vector<double>> Ez =
      ChosenCT(Cn, Dx, Dy, nghost, rec, sl, rs, OptP);

  for (int j = nghost; j < Cn1.ny - nghost; ++j) {
    for (int i = nghost; i <= Cn1.nx - nghost; ++i) {
      Cn1.Bxf[j][i] =
          Cn.Bxf[j][i] - (Dt / Dy) * (Ez[j + 1 - nghost][i - nghost] -
                                      Ez[j - nghost][i - nghost]);
    }
  }

  for (int j = nghost; j <= Cn1.ny - nghost; ++j) {
    for (int i = nghost; i < Cn1.nx - nghost; ++i) {
      Cn1.Byf[j][i] =
          Cn.Byf[j][i] + (Dt / Dx) * (Ez[j - nghost][i + 1 - nghost] -
                                      Ez[j - nghost][i - nghost]);
    }
  }
}
