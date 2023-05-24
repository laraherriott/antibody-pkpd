function dy = ODEs(t, y, k_in, k_olig_inc, k_olig_sep, k_clear_Abeta, ...
    k_onPP, k_off_ma0, k_off_ma1, k_plaque_inc, k_plaque_sep, k_clear_olig, ...
    k_clear_P, k_onPD, k_off_ma2, k_synth_FcR, k_clear_FcR, k_onPF, k_offPF, ...
    k_ADCP, clearance, k_mAb_transport_back, k_mAb_transport, k_mAbcomplex_clear, ...
    dose_list)

    % Variables
    ABeta = y(1);
    Oligomer = y(2);
    Plaque = y(3);
    FcR = y(4);
    mAb_plasma = y(5);
    mAb = y(6);
    ABeta_mAb = y(7);
    Oligomer_mAb = y(8);
    Plaque_mAb = y(9);
    Oligomer_mAb_FcR = y(10);
    Plaque_mAb_FcR = y(11);

    % Differentials
    dABeta = k_in - k_olig_inc * ABeta + k_olig_sep * Oligomer - k_clear_Abeta * ABeta - k_onPP * ABeta * mAb + k_off_ma0 * ABeta_mAb;
    dOligomer = k_olig_inc * ABeta - k_olig_sep * Oligomer - k_plaque_inc * Oligomer + k_plaque_sep * Plaque - k_clear_olig * Oligomer - k_onPP * Oligomer * mAb + k_off_ma1 * Oligomer_mAb;
    dPlaque = k_plaque_inc * Oligomer - k_plaque_sep * Plaque - k_clear_P * Plaque - k_onPD * Plaque * mAb + k_off_ma2 * Plaque_mAb;
    dFcR = k_synth_FcR - k_clear_FcR * FcR - k_onPF * Oligomer_mAb * FcR + k_offPF * Oligomer_mAb_FcR - k_onPF * Plaque_mAb * FcR + k_offPF * Plaque_mAb_FcR + k_ADCP * Oligomer_mAb_FcR + k_ADCP * Plaque_mAb_FcR;
    dmAb_plasma = dosefn(dose_list, t) - clearance * mAb_plasma + k_mAb_transport_back * mAb - k_mAb_transport * mAb_plasma;
    dmAb = -k_mAb_transport_back * mAb + k_mAb_transport * mAb_plasma - k_onPP * ABeta * mAb + k_off_ma0 * ABeta_mAb - k_onPP * Oligomer * mAb + k_off_ma1 * Oligomer_mAb - k_onPD * Plaque * mAb + k_off_ma2 * Plaque_mAb - k_mAbcomplex_clear * mAb;
    dABeta_mAb = k_onPP * ABeta * mAb - k_off_ma0 * ABeta_mAb - k_mAbcomplex_clear * ABeta_mAb;
    dOligomer_mAb = k_onPP * Oligomer * mAb - k_off_ma1 * Oligomer_mAb - k_mAbcomplex_clear * Oligomer_mAb + k_offPF * Oligomer_mAb_FcR - k_onPF * Oligomer_mAb * FcR;
    dPlaque_mAb = k_onPD * Plaque * mAb - k_off_ma2 * Plaque_mAb - k_onPF * Plaque_mAb * FcR + k_offPF * Plaque_mAb_FcR;
    dOligomer_mAb_FcR = k_onPF * Oligomer_mAb * FcR - k_offPF * Oligomer_mAb_FcR - k_ADCP * Oligomer_mAb_FcR;
    dPlaque_mAb_FcR = k_onPF * Plaque_mAb * FcR - k_offPF * Plaque_mAb_FcR - k_ADCP * Plaque_mAb_FcR;

    dy = [dABeta; dOligomer; dPlaque; dFcR; dmAb_plasma; dmAb; dABeta_mAb; dOligomer_mAb; dPlaque_mAb; dOligomer_mAb_FcR; dPlaque_mAb_FcR];
end