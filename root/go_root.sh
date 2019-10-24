#!/bin/bash
#
#

tag="Ds"
##tag="etac"

echo "Start of job at " `date`
log=log_${tag}
echo "Output written to " ${log}

echo "Start of job at " `date` > ${log}

root -b <<EOF  >> ${log}

.L ${MYROOTLIB}/boot_library.C
.L ${MYROOTLIB}/io_library.C
.L ${MYROOTLIB}/run_alpha_s_lib.C
.L ${MYROOTLIB}/chisq_prob.C

double a_fm_cut     = 0.8 ;
double m_pi_gev_cut = 0.9 ;

char fileN[] = "loop_data_${tag}.txt"  ;

.x loop_boot_myfit_3param_2anal_B.C
EOF

echo "End of job at " `date` >> ${log}

echo "Summary"
grep SUM ${log}  | sed 's/SUM //'

exit 0
