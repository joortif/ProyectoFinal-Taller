#ifndef PROYECTO_GLOBALS_H
#define PROYECTO_GLOBALS_H

extern double **eopdata;
extern double **Cnm;
extern double **Snm;
struct AuxParam{
    double Mjd_TT;
    int n;
    int m;
    double Mjd_UTC;
};
extern AuxParam auxParam;
extern int n_eqn;


#endif //PROYECTO_GLOBALS_H
