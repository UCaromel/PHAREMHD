#ifndef COMPUTE_J_HPP_
#define COMPUTE_J_HPP_

// J has 1 more ghost cells for laplacian calculations in Flux vector
template<typename Variables>
void ComputeJ(Variables& Vn, double Dx, double Dy, int nghost) {

    for (int j = nghost; j <= Vn.ny - nghost; ++j)
    {
        for (int i = nghost; i < Vn.nx - nghost; ++i)
        {
            Vn.Jx[j + 1][i + 1] = (Vn.Bz[j][i] - Vn.Bz[j - 1][i]) / (Dy);
        }
    }

    for (int j = nghost; j < Vn.ny - nghost; ++j)
    {
        for (int i = nghost; i <= Vn.nx - nghost; ++i)
        {
            Vn.Jy[j + 1][i + 1] = - (Vn.Bz[j][i] - Vn.Bz[j][i - 1]) / (Dx);
        }
    }

    for (int j = nghost; j <= Vn.ny - nghost; ++j)
    {
        for (int i = nghost; i <= Vn.nx - nghost; ++i)
        {
            // In 1D (no CT), use cell centered Bx and By
            Vn.Jz[j + 1][i + 1] =  (Vn.By[j][i] - Vn.By[j][i - 1]) / (Dx) - (Vn.Bx[j][i] - Vn.Bx[j - 1][i]) / (Dy);
            // Vn.Jz[j + 1][i + 1] = (Vn.Byf[j][i] - Vn.Byf[j][i - 1]) / (Dx) - (Vn.Bxf[j][i] - Vn.Bxf[j - 1][i]) / (Dy);
        }
    }
}

#endif // COMPUTE_J_HPP_