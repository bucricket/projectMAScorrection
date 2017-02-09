extern void rttov_load_inst_(
    int* inst_id,
    char* opts_str,
    int* nchannels,
    int  channels[],
    int l);

extern void rttov_call_direct_(
    int* err,
    int* inst_id,
    int channel_list[],
    int datetimes[],            // profile dates/times                                     [nprofiles][6]
    double angles[],            // satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
    double surfgeom[],          // lat, lon, elevation                                     [nprofiles][3]
    int surftype[],             // surftype, watertype                                     [nprofiles][2]
    double skin[],              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
    double s2m[],               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
    double simplecloud[],       // ctp, cfraction                                          [nprofiles][2]
    int icecloud[],             // ish, idg                                                [nprofiles][2]
    double zeeman[],            // Be, cosbk                                               [nprofiles][2]
    double p[],                 // pressure                                                [nprofiles][nlevels]
    double t[],                 // temperature                                             [nprofiles][nlevels]
    int* gas_units,             // units for gas profiles
    int gas_id[],               // gas ID list                                             [ngases]
    double gases[],             // gas profiles                                            [ngases][nprofiles][nlevels]
    double surfemisrefl[],      // input/output surface emissivities/BRDFs                 [2][nprofiles][nchannels]
    double btrefl[],            // output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
    double rads[],              // output radiances                                        [nprofiles][nchannels]
    int* nchannels, int* ngases, int* nlevels, int* nprofiles);

extern void rttov_call_k_(
    int* err,
    int* inst_id,
    int channel_list[],
    int datetimes[],            // profile dates/times                                     [nprofiles][6]
    double angles[],            // satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
    double surfgeom[],          // lat, lon, elevation                                     [nprofiles][3]
    int surftype[],             // surftype, watertype                                     [nprofiles][2]
    double skin[],              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
    double skin_k[],            // output skin K                                           [nprofiles][nchannels][9]
    double s2m[],               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
    double s2m_k[],             // output 2m K                                             [nprofiles][nchannels][6]
    double simplecloud[],       // ctp, cfraction                                          [nprofiles][2]
    double simplecloud_k[],     // output ctp, cfraction K                                 [nprofiles][nchannels][2]
    int icecloud[],             // ish, idg                                                [nprofiles][2]
    double zeeman[],            // Be, cosbk                                               [nprofiles][2]
    double p[],                 // pressure                                                [nprofiles][nlevels]
    double p_k[],               // output pressure K                                       [nprofiles][nchannels][nlevels]
    double t[],                 // temperature                                             [nprofiles][nlevels]
    double t_k[],               // output temperature K                                    [nprofiles][nchannels][nlevels]
    int* gas_units,             // units for gas profiles
    int gas_id[],               // gas ID list                                             [ngases]
    double gases[],             // gas profiles                                            [ngases][nprofiles][nlevels]
    double gases_k[],           // output gas profiles K                                   [ngases][nprofiles][nchannels][nlevels]
    double surfemisrefl[],      // input/output surface emissivities/BRDFs                 [2][nprofiles][nchannels]
    double surfemisrefl_k[],    // input/output surface emissivities/BRDFs K               [2][nprofiles][nchannels]
    double btrefl[],            // output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
    double rads[],              // output radiances                                        [nprofiles][nchannels]
    double btrefl_k[],          // input BT perturbations                                  [nprofiles][nchannels]
    double rads_k[],            // input radiance perturbations                            [nprofiles][nchannels]
    int* nchannels, int* ngases, int* nlevels, int* nprofiles);

extern void rttov_drop_inst_(int* err, int* inst_id);

extern void rttov_drop_all_(int* err);

extern void rttov_set_options_(int* err, int* inst_id, char* opts_str, int l);

extern void rttov_print_options_(int* err, int* inst_id);


extern void rttov_ir_emis_atlas_setup_(
    int* err,
    char* path,
    int* month,
    int* version,
    int* inst_id,
    int* ang_corr,
    int l);

extern void rttov_mw_emis_atlas_setup_(
    int* err,
    char* path,
    int* month,
    int* version,
    int* inst_id,
    int l);

extern void rttov_brdf_atlas_setup_(
    int* err,
    char* path,
    int* month,
    int* version,
    int* inst_id,
    int l);

extern void rttov_ir_emis_atlas_dealloc_();

extern void rttov_mw_emis_atlas_dealloc_();

extern void rttov_brdf_atlas_dealloc_();


extern void rttov_get_rad_clear_(int* err, int* inst_id, double rad_clear[], int* nchanprof);

extern void rttov_get_rad_total_(int* err, int* inst_id, double rad_total[], int* nchanprof);

extern void rttov_get_bt_clear_(int* err, int* inst_id, double bt_clear[], int* nchanprof);

extern void rttov_get_bt_(int* err, int* inst_id, double bt[], int* nchanprof);

extern void rttov_get_refl_clear_(int* err, int* inst_id, double refl_clear[], int* nchanprof);

extern void rttov_get_refl_(int* err, int* inst_id, double refl[], int* nchanprof);

extern void rttov_get_rad_cloudy_(int* err, int* inst_id, double rad_cloudy[], int* nchanprof);

extern void rttov_get_overcast_(int* err, int* inst_id, double overcast[], int* nchanprof, int* nlayers);


extern void rttov_get_rad2_upclear_(int* err, int* inst_id, double rad2_upclear[], int* nchanprof);

extern void rttov_get_rad2_dnclear_(int* err, int* inst_id, double rad2_dnclear[], int* nchanprof);

extern void rttov_get_rad2_refldnclear_(int* err, int* inst_id, double rad2_refldnclear[], int* nchanprof);

extern void rttov_get_rad2_up_(int* err, int* inst_id, double rad2_up[], int* nchanprof, int* nlayers);

extern void rttov_get_rad2_down_(int* err, int* inst_id, double rad2_down[], int* nchanprof, int* nlayers);

extern void rttov_get_rad2_surf_(int* err, int* inst_id, double rad2_surf[], int* nchanprof, int* nlayers);


extern void rttov_get_tau_total_(int* err, int* inst_id, double tau_total[], int* nchanprof);

extern void rttov_get_tau_levels_(int* err, int* inst_id, double tau_levels[], int* nchanprof, int* nlevels);

extern void rttov_get_tausun_total_path2_(int* err, int* inst_id, double tau_total_path2[], int* nchanprof);

extern void rttov_get_tausun_levels_path2_(int* err, int* inst_id, double tau_levels_path2[], int* nchanprof, int* nlevels);

extern void rttov_get_tausun_total_path1_(int* err, int* inst_id, double tau_total_path1[], int* nchanprof);

extern void rttov_get_tausun_levels_path1_(int* err, int* inst_id, double tau_levels_path1[], int* nchanprof, int* nlevels);


// The following are used by the Rttov and RttovSafe classes to interrogate the RTTOV coefficients structure
extern void rttov_get_coef_val_i0_(int* err, int* inst_id, char* varch, int* i0, int l);

extern void rttov_get_coef_val_r0_(int* err, int* inst_id, char* varch, double* r0, int l);

extern void rttov_get_coef_val_c0_(int* err, int* inst_id, char* varch, char* c0, int l);

extern void rttov_get_coef_val_i1_(int* err, int* inst_id, char* varch, int* m, int* i1, int l);

extern void rttov_get_coef_val_i2_(int* err, int* inst_id, char* varch, int* m, int* n, int* i2, int l);

extern void rttov_get_coef_val_r1_(int* err, int* inst_id, char* varch, int* m, double* r1, int l);




