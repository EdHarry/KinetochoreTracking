/*
 * mmfMex_data.h
 *
 * Code generation for function 'mmfMex_data'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

#ifndef __MMFMEX_DATA_H__
#define __MMFMEX_DATA_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "mmfMex_types.h"

/* Variable Declarations */
extern const volatile char_T *emlrtBreakCheckR2012bFlagVar;
extern emlrtRSInfo s_emlrtRSI;
extern emlrtRSInfo t_emlrtRSI;
extern emlrtRSInfo v_emlrtRSI;
extern emlrtRSInfo w_emlrtRSI;
extern emlrtRSInfo x_emlrtRSI;
extern emlrtRSInfo bb_emlrtRSI;
extern emlrtRSInfo cb_emlrtRSI;
extern emlrtRSInfo eb_emlrtRSI;
extern emlrtRSInfo gb_emlrtRSI;
extern emlrtRSInfo hb_emlrtRSI;
extern emlrtRSInfo ib_emlrtRSI;
extern emlrtRSInfo jb_emlrtRSI;
extern emlrtRSInfo mb_emlrtRSI;
extern emlrtRSInfo qb_emlrtRSI;
extern emlrtRSInfo rb_emlrtRSI;
extern emlrtRSInfo sb_emlrtRSI;
extern emlrtRSInfo vb_emlrtRSI;
extern emlrtRSInfo wb_emlrtRSI;
extern emlrtRSInfo xb_emlrtRSI;
extern emlrtRSInfo yb_emlrtRSI;
extern emlrtRSInfo ac_emlrtRSI;
extern emlrtRSInfo bc_emlrtRSI;
extern emlrtRSInfo cc_emlrtRSI;
extern emlrtRSInfo fc_emlrtRSI;
extern emlrtRSInfo gc_emlrtRSI;
extern emlrtRSInfo hc_emlrtRSI;
extern emlrtRSInfo ic_emlrtRSI;
extern emlrtRSInfo jc_emlrtRSI;
extern emlrtRSInfo kc_emlrtRSI;
extern emlrtRSInfo lc_emlrtRSI;
extern emlrtRSInfo mc_emlrtRSI;
extern emlrtRSInfo nc_emlrtRSI;
extern emlrtRSInfo oc_emlrtRSI;
extern emlrtRSInfo pc_emlrtRSI;
extern emlrtRSInfo qc_emlrtRSI;
extern emlrtRSInfo rc_emlrtRSI;
extern emlrtRSInfo sc_emlrtRSI;
extern emlrtRSInfo tc_emlrtRSI;
extern emlrtRSInfo uc_emlrtRSI;
extern emlrtRSInfo vc_emlrtRSI;
extern emlrtRSInfo wc_emlrtRSI;
extern emlrtRSInfo xc_emlrtRSI;
extern emlrtRSInfo yc_emlrtRSI;
extern emlrtRSInfo ad_emlrtRSI;
extern emlrtRSInfo bd_emlrtRSI;
extern emlrtRSInfo cd_emlrtRSI;
extern emlrtRSInfo dd_emlrtRSI;
extern emlrtRSInfo fd_emlrtRSI;
extern emlrtRSInfo gd_emlrtRSI;
extern emlrtRSInfo hd_emlrtRSI;
extern emlrtRSInfo id_emlrtRSI;
extern emlrtRSInfo jd_emlrtRSI;
extern emlrtRSInfo kd_emlrtRSI;
extern emlrtRSInfo ld_emlrtRSI;
extern emlrtRSInfo sd_emlrtRSI;
extern emlrtRSInfo td_emlrtRSI;
extern emlrtRSInfo de_emlrtRSI;
extern emlrtRSInfo ee_emlrtRSI;
extern emlrtRSInfo fe_emlrtRSI;
extern emlrtRSInfo mf_emlrtRSI;
extern emlrtRSInfo tf_emlrtRSI;
extern emlrtRSInfo uf_emlrtRSI;
extern emlrtRSInfo vf_emlrtRSI;
extern emlrtRSInfo wf_emlrtRSI;
extern emlrtRSInfo xf_emlrtRSI;
extern emlrtRSInfo bg_emlrtRSI;
extern emlrtRSInfo dg_emlrtRSI;
extern emlrtRSInfo eg_emlrtRSI;
extern emlrtRSInfo fg_emlrtRSI;
extern emlrtRSInfo gg_emlrtRSI;
extern emlrtRSInfo rg_emlrtRSI;
extern emlrtRSInfo sg_emlrtRSI;
extern emlrtRSInfo tg_emlrtRSI;
extern emlrtRSInfo ug_emlrtRSI;
extern emlrtRSInfo vg_emlrtRSI;
extern emlrtRSInfo wg_emlrtRSI;
extern emlrtRSInfo xg_emlrtRSI;
extern emlrtRSInfo yg_emlrtRSI;
extern emlrtRSInfo ah_emlrtRSI;
extern emlrtRSInfo bh_emlrtRSI;
extern emlrtRSInfo ch_emlrtRSI;
extern emlrtRSInfo dh_emlrtRSI;
extern emlrtRSInfo eh_emlrtRSI;
extern emlrtRSInfo fh_emlrtRSI;
extern emlrtRSInfo gh_emlrtRSI;
extern emlrtRSInfo hh_emlrtRSI;
extern emlrtRSInfo ih_emlrtRSI;
extern emlrtRSInfo jh_emlrtRSI;
extern emlrtRSInfo nh_emlrtRSI;
extern emlrtRSInfo ph_emlrtRSI;
extern emlrtRSInfo qh_emlrtRSI;
extern emlrtRSInfo rh_emlrtRSI;
extern emlrtRSInfo sh_emlrtRSI;
extern emlrtRSInfo th_emlrtRSI;
extern emlrtRSInfo uh_emlrtRSI;
extern emlrtRSInfo vh_emlrtRSI;
extern emlrtRSInfo wh_emlrtRSI;
extern emlrtRSInfo xh_emlrtRSI;
extern emlrtRSInfo yh_emlrtRSI;
extern emlrtRSInfo ei_emlrtRSI;
extern emlrtRSInfo fi_emlrtRSI;
extern emlrtRSInfo gi_emlrtRSI;
extern emlrtRSInfo hi_emlrtRSI;
extern emlrtRSInfo ii_emlrtRSI;
extern emlrtRSInfo ji_emlrtRSI;
extern emlrtRSInfo ki_emlrtRSI;
extern emlrtRSInfo li_emlrtRSI;
extern emlrtRSInfo mi_emlrtRSI;
extern emlrtRSInfo ni_emlrtRSI;
extern emlrtRSInfo oi_emlrtRSI;
extern emlrtRSInfo ik_emlrtRSI;
extern emlrtRSInfo kk_emlrtRSI;
extern emlrtRSInfo lk_emlrtRSI;
extern emlrtRSInfo mk_emlrtRSI;
extern emlrtRSInfo nk_emlrtRSI;
extern emlrtRSInfo jl_emlrtRSI;
extern emlrtRSInfo kl_emlrtRSI;
extern emlrtRSInfo ll_emlrtRSI;
extern emlrtRSInfo ml_emlrtRSI;
extern emlrtRSInfo nl_emlrtRSI;
extern emlrtRSInfo ol_emlrtRSI;
extern emlrtRSInfo pl_emlrtRSI;
extern emlrtRSInfo ql_emlrtRSI;
extern emlrtRSInfo dm_emlrtRSI;
extern emlrtRSInfo em_emlrtRSI;
extern emlrtRSInfo fm_emlrtRSI;
extern emlrtRSInfo gm_emlrtRSI;
extern emlrtRSInfo hm_emlrtRSI;
extern emlrtRSInfo im_emlrtRSI;
extern emlrtRSInfo jm_emlrtRSI;
extern emlrtRSInfo km_emlrtRSI;
extern emlrtRSInfo lm_emlrtRSI;
extern emlrtRSInfo mm_emlrtRSI;
extern emlrtRSInfo nm_emlrtRSI;
extern emlrtRSInfo om_emlrtRSI;
extern emlrtRSInfo pm_emlrtRSI;
extern emlrtRSInfo qm_emlrtRSI;
extern emlrtRSInfo rm_emlrtRSI;
extern emlrtRSInfo sm_emlrtRSI;
extern emlrtRSInfo tm_emlrtRSI;
extern emlrtRSInfo ym_emlrtRSI;
extern emlrtRSInfo an_emlrtRSI;
extern emlrtRSInfo bn_emlrtRSI;
extern emlrtRSInfo cn_emlrtRSI;
extern emlrtRSInfo dn_emlrtRSI;
extern emlrtRSInfo en_emlrtRSI;
extern emlrtRSInfo fn_emlrtRSI;
extern emlrtRSInfo gn_emlrtRSI;
extern emlrtRSInfo hn_emlrtRSI;
extern emlrtRSInfo in_emlrtRSI;
extern emlrtRSInfo jn_emlrtRSI;
extern emlrtRSInfo kn_emlrtRSI;
extern emlrtRSInfo ln_emlrtRSI;
extern emlrtRSInfo mn_emlrtRSI;
extern emlrtRSInfo nn_emlrtRSI;
extern emlrtRSInfo on_emlrtRSI;
extern emlrtRSInfo pn_emlrtRSI;
extern emlrtRSInfo qn_emlrtRSI;
extern emlrtRSInfo rn_emlrtRSI;
extern emlrtRSInfo sn_emlrtRSI;
extern emlrtRSInfo tn_emlrtRSI;
extern emlrtRSInfo vn_emlrtRSI;
extern emlrtRSInfo ao_emlrtRSI;
extern emlrtRSInfo co_emlrtRSI;
extern emlrtRSInfo do_emlrtRSI;
extern emlrtRSInfo eo_emlrtRSI;
extern emlrtRSInfo up_emlrtRSI;
extern emlrtMCInfo m_emlrtMCI;
extern emlrtMCInfo n_emlrtMCI;
extern emlrtRTEInfo b_emlrtRTEI;
extern emlrtRTEInfo c_emlrtRTEI;
extern emlrtRTEInfo t_emlrtRTEI;
extern emlrtRTEInfo ab_emlrtRTEI;
extern emlrtRTEInfo bb_emlrtRTEI;
extern emlrtRTEInfo cb_emlrtRTEI;
extern emlrtRTEInfo db_emlrtRTEI;
extern emlrtRSInfo iq_emlrtRSI;
#endif
/* End of code generation (mmfMex_data.h) */
