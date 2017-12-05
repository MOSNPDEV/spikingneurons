#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _cacumm_reg();
extern void _cagk_reg();
extern void _cal2_reg();
extern void _can2_reg();
extern void _cat_reg();
extern void _kaprox_reg();
extern void _kd_reg();
extern void _kdrca1_reg();
extern void _km_reg();
extern void _na3n_reg();

modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," cacumm.mod");
fprintf(stderr," cagk.mod");
fprintf(stderr," cal2.mod");
fprintf(stderr," can2.mod");
fprintf(stderr," cat.mod");
fprintf(stderr," kaprox.mod");
fprintf(stderr," kd.mod");
fprintf(stderr," kdrca1.mod");
fprintf(stderr," km.mod");
fprintf(stderr," na3n.mod");
fprintf(stderr, "\n");
    }
_cacumm_reg();
_cagk_reg();
_cal2_reg();
_can2_reg();
_cat_reg();
_kaprox_reg();
_kd_reg();
_kdrca1_reg();
_km_reg();
_na3n_reg();
}
