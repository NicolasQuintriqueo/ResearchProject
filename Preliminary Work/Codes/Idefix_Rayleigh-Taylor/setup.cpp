#include "idefix.hpp"
#include "setup.hpp"


real g_y;
real rho_heavy;
real rho_light;

//potencial phi para cada punto de la grid en Y
void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1,
               IdefixArray1D<real>& x2, IdefixArray1D<real>& x3,
               IdefixArray3D<real>& phi) {
    idefix_for("Potential", 0, data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
            phi(k, j, i) = g_y * x2(j);
        });
}

//constructor
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
    rho_heavy = input.Get<real>("Setup", "rho_heavy", 1);
    rho_light = input.Get<real>("Setup", "rho_light", 1);
    g_y = input.Get<real>("Setup", "g_y", 1);

    data.gravity->EnrollPotential(&Potential);
}


void Setup::InitFlow(DataBlock &data) {
    DataBlockHost d(data);

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {

                real yInterface = HALF_F * (1 + 0.1 * (cos(4.0 * M_PI * d.x[IDIR](i))));

                d.Vc(RHO, k, j, i) = (d.x[JDIR](j) > yInterface) ? rho_heavy : rho_light;

                EXPAND(
                    d.Vc(VX1, k, j, i) = ZERO_F; ,
                    d.Vc(VX2, k, j, i) = ZERO_F; ,
                    d.Vc(VX3, k, j, i) = ZERO_F;
                )

                #if HAVE_ENERGY
                //para un problema adiabático, se ajustaria la presion para lograr el equilibrio hidrostático. 
                    d.Vc(PRS, k, j, i) = 2.5;
                #endif
            }
        }
    }
    d.SyncToDevice();
}

void MakeAnalysis(DataBlock & data) {
}