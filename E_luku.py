# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 22:36:00 2022

@author: laukkara
"""

import numpy as np


A_netto = 1.0

f_sahko = 1.2
f_kaukolampo = 0.5
f_kaukojaahdytys = 0.28
f_polttoaine = 1.0
f_uusiutuvat = 0.5

dt = 24.0 * np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
rho_i = 1.2
c_pi = 1000.0

###########




T_s = np.array([21.0])

T_u = np.array([-5, -4, -1, 3, 5, 10, 12, 15, 10, 3, -1, -2.0])







# 3.2 Rakennusvaipan johtumislämpöhäviöt

U_us = 0.17
A_us = 200.0
H_ulkoseina = U_us * A_us

U_yp = 0.09
A_yp = 120.0
H_ylapohja = U_yp * A_yp

U_ap = 0.19
A_ap = A_yp
H_alapohja = U_ap * A_ap



ikk = [{'A': 10.0,
     'U': 1.0,
     'ilmansuunta': 0.0,
     'phi': 30.0,
     'alpha': 10.0,
     'beta': 10.0},
    {'A': 10.0,
     'U': 1.0,
     'ilmansuunta': 90.0,
     'phi': 30.0,
     'alpha': 10.0,
     'beta': 10.0},
    {'A': 10.0,
     'U': 1.0,
     'ilmansuunta': 180.0,
     'phi': 30.0,
     'alpha': 10.0,
     'beta': 10.0},
    {'A': 10.0,
     'U': 1.0,
     'ilmansuunta': 270.0,
     'phi': 30.0,
     'alpha': 10.0,
     'beta': 10.0}]

A_ikk = np.sum([i['A'] for i in ikk])

H_ikkuna = np.sum([i['U']*i['A'] for i in ikk])

U_ovi = 1.0
A_ovi = 3.0
H_ovi = U_ovi * A_ovi

U_muu = 0.0
A_muu = 0.0
H_muu = U_muu * A_muu

psi_ks = 0.04
L_ks = 100.0
H_kylmasilta = psi_ks * L_ks


T_u_vuosi = np.mean(T_u)

dT_maa_vuosi = 5.0
T_maa_vuosi = T_u_vuosi + dT_maa_vuosi

dT_maa_kuukausi = np.array([0.0, -1, -2, -3, -3, -2, 0, 1, 2, 3, 3, 2])
T_maa_kuukausi = T_maa_vuosi + dT_maa_kuukausi
T_ap = T_maa_kuukausi


Q_ulkoseina = H_ulkoseina * (T_s - T_u) * dt / 1000
Q_ylapohja = H_ylapohja * (T_s - T_u) * dt / 1000
Q_alapohja = H_alapohja * (T_s - T_ap) * dt / 1000
Q_ikkuna = H_ikkuna * (T_s - T_u) * dt / 1000
Q_ovi = H_ovi * (T_s - T_u) * dt / 1000
Q_muu = H_muu * (T_s - T_u) * dt / 1000
Q_kylmasilta = H_kylmasilta * (T_s - T_u) * dt / 1000



Q_joht = Q_ulkoseina \
        + Q_ylapohja \
        + Q_alapohja \
        + Q_ikkuna \
        + Q_ovi \
        + Q_muu \
        + Q_kylmasilta



# 3.3 Vuotoilman lämpenemisen lämpöenergian tarve

q_50 = 4.0
A_vaippa = np.sum([A_us, A_yp, A_ap, A_ikk, A_ovi, A_muu])

x = 35.0


q_v_vuotoilma = (q_50/(3600.0 * x)) * A_vaippa

Q_vuotoilma = rho_i * c_pi * q_v_vuotoilma * (T_s - T_u) * dt / 1000




# 3.4 Ilmanvaihdon lämmitysenergian nettotarve
A_netto = A_ap

t_d = 1.0
t_v = 1.0

q_v_tulo_omin = 0.4/1000
q_v_tulo = q_v_tulo_omin * A_netto

T_sp = 18.0
dT_puhallin = 0.5

idx_lto = np.ones(12)
idx_lto[5:8] = 0.0
eta_a_ivkone = 0.55 * idx_lto
q_v_poisto = q_v_tulo


phi_lto = eta_a_ivkone*t_d*t_v*rho_i*c_pi*q_v_poisto*(T_s - T_u)
T_lto = T_u + phi_lto / (t_d*t_v*rho_i*c_pi*q_v_tulo)

Q_iv = t_d * t_v * rho_i * c_pi * q_v_tulo * ((T_sp-dT_puhallin)-T_lto) * dt / 1000
Q_iv = np.maximum(Q_iv, 0.0)




# 3.5 Tuloilman ja korvausilman lämmitysenergian tarve

Q_iv_tuloilma = t_d*t_v*rho_i*c_pi*q_v_tulo*(T_s-T_u)*dt/1000

q_v_korvausilma = t_d*t_v*q_v_poisto - t_d*t_v*q_v_tulo
Q_iv_korvausilma = rho_i*c_pi*q_v_korvausilma*(T_s-T_u)*dt/1000



# 3.6 Ilmanvaihdosta talteenotettu energia

Q_lto = t_d*t_v*rho_i*c_pi*q_v_tulo*(T_lto-T_u)*dt/1000


# 3.7 Lämpimän käyttöveden lämmitysenergian nettotarve

rho_v = 1000.0
c_pv = 4200.0

T_lkv = 55.0
T_kv = 5.0

V_lkv = 0.6*A_netto

Q_lkv_netto = rho_v * c_pv * V_lkv * (T_lkv - T_kv)/(1000*3600)
Q_lkv_netto = Q_lkv_netto * (dt/np.sum(dt))



# 4.2 Valaistuksen sähköenergian kulutus
tau_d = 24
tau_w = 7
k_valaistus = 0.1
P_valaistus = 6.0
W_valaistus = k_valaistus * P_valaistus * (tau_d/24.0)*(tau_w/7)*(8760/1000)
W_valaistus = W_valaistus * (dt/np.sum(dt))



# 5.1 Lämpökuorma henkilöistä
k_henk = 0.6
P_henk = 2.0
Q_henk = k_henk * P_henk * (tau_d/24.0)*(tau_w/7)*(8760/1000)
Q_henk = Q_henk * (dt/np.sum(dt))



# 5.2 Lämpökuorma valaistuksesta ja sähkölaitteista

k_kl = 0.6
P_kl = 3.0
W_kuluttajalaitteet = k_kl * P_kl * (tau_d/24.0)*(tau_w/7)*(8760/1000)
W_kuluttajalaitteet = W_kuluttajalaitteet * (dt/np.sum(dt))

Q_säh = W_valaistus + W_kuluttajalaitteet




# 5.3 Ikkunoiden kautta rakennukseen tuleva auringon säteilyenergia

g_kohtisuora = 0.7
g = 0.9 * g_kohtisuora

F_kehä = 0.75
F_verho = 0.75

F_varjostus = F_ymparisto * F_ylavarjostus * F_sivuvarjostus
F_läpäisy = F_kehä * F_verho * F_varjostus



















# 3.1 Tilojen lämmitysenergian nettotarve


Q_tila = Q_joht \
        + Q_vuotoilma \
        + Q_iv_tuloilma \
        + Q_iv_korvausilma

Q_lammitys_tilat_netto = Q_tila - Q_sis_lampo



# 



E = (f_kaukolampo*Q_kaukolampo \
     + f_kaukojaahdytys*Q_kaukojaahdytys \
     + f_polttoaine*Q_polttoaine \
     + f_sahko*W_sahko) / A_netto

    
RAK_ek = (Q_lammitys_tilat \
        + Q_lammitys_iv \
        + Q_lammitys_lkv \
        + Q_jk \
        + W_tilat \
        + W_ilmanvaihto \
        + W_lkv_pumppu \
        + W_jaahd_apu \
        + W_kuluttajalaitteet \
        + W_valaistus) / A_netto