# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import numpy as np
import matplotlib.pyplot as plt

#%%

# Constante gravitationnelle (CODATA recommandée)
G = 6.67408e-11  # m^3 kg^-1 s^-2

def demi_grand_axe_1(r_a, r_p):
    return (r_a+r_p)/2

def demi_grand_axe_2(r_a, e):
    return r_a/(1+e)

def demi_grand_axe_3(r_p, e):
    return r_p/(1-e)

def alt_apogee(a,e, r_T = 6378):
    return (a*(1+e)) -  r_T

def demi_grand_axe_periode(T, M, G=G):
    # Entrée: T en jours, M en kg. Retour: demi-grand axe en km.
    T_s = T * 86400
    a_m = ((T_s**2 * G * M)/(4 * (np.pi**2)))**(1/3)
    return a_m / 1e3

def periode_demi_grand_axe(a, M, G=G):
    # Entrée: a en km, M en kg. Retour: période en jours.
    a_m = a * 1e3
    T_s = 2*np.pi*np.sqrt((a_m**3)/(G*M))
    return T_s / 86400

def mouvement_moyen(T):
    return (2*np.pi)/T
    
#%%
r_T = 6378*10**3
r_p = 620*10**3 + r_T
e = 0.04

a = demi_grand_axe_3(r_p, e)/1000
print(a, "km")
print(alt_apogee(a, e), "km")

#%%

M = 2*10**30 #masse du soleil
M_T = 6*10**24 #masse terre

print(demi_grand_axe_periode(87.969, M))
print(demi_grand_axe_periode(686.886, M))
print(periode_demi_grand_axe(1426700000, M))
print(periode_demi_grand_axe(384399, M_T))

#%%
r_T = 6378 #km
r_p = 620 + r_T #km
e = 0.04
a = demi_grand_axe_3(r_p, e)
print(periode_demi_grand_axe(a, M_T)*24*60, "min")
print(mouvement_moyen(periode_demi_grand_axe(a, M_T)*24*3600), "rad/s")

#%%

def vitesse_orbite_circulaire(r, mu= 398600.4418, G = 6.67408e-11):
    return np.sqrt((mu)/r)

def duree_orbite_circulaire(v, rayon):
    return v * (2*np.pi*rayon)

R_T = 6378

#arevoir

#utiliser le temps sidéral /!\

print(vitesse_orbite_circulaire(200+r_T), "km/s", duree_orbite_circulaire(vitesse_orbite_circulaire(200+r_T), 200+r_T)/(3600*24), "secondes", (60*24)/(duree_orbite_circulaire(vitesse_orbite_circulaire(200+r_T), 200+r_T)/3600))
print(vitesse_orbite_circulaire(1000+r_T), "km/s",duree_orbite_circulaire(vitesse_orbite_circulaire(1000+r_T), 1000+r_T)/(3600*24), "secondes", (60*24)/(duree_orbite_circulaire(vitesse_orbite_circulaire(1000+r_T), 200+r_T)/3600))
print(vitesse_orbite_circulaire(10000+r_T), "km/s",duree_orbite_circulaire(vitesse_orbite_circulaire(10000+r_T), 10000+r_T)/(3600*24), "secondes", (60*24)/(duree_orbite_circulaire(vitesse_orbite_circulaire(10000+r_T), 200+r_T)/3600))
print(vitesse_orbite_circulaire(35786+r_T), "km/s",duree_orbite_circulaire(vitesse_orbite_circulaire(35786+r_T), 35786+r_T)/(3600*24), "secondes", (60*24)/(duree_orbite_circulaire(vitesse_orbite_circulaire(35786+r_T), 200+r_T)/3600))

#%%
r_T = 6378 #km
alt_perigee = 300 #km
alt_apogee = 35786 #km
M_T = 6*10**24 #masse terre #kg

r_p = alt_perigee + r_T
r_a = alt_apogee + r_T

demi_grand_axe = demi_grand_axe_1(r_a, r_p) #km
excentricite = 1 - (r_p/demi_grand_axe)
periode = periode_demi_grand_axe(demi_grand_axe, M=M_T) #jours

print(demi_grand_axe, "km")
print(excentricite)
print(periode, "jours")

def vitesse_selon_rayon(r, a, mu_terre = 398600.4418):
    return np.sqrt(mu_terre*((2/r)-(1/a)))

vitesse_perigee = vitesse_selon_rayon(r_p, demi_grand_axe) #km/s
vitesse_apogee = vitesse_selon_rayon(r_a, demi_grand_axe) #km/s

print("en ellipse au périgée", vitesse_perigee, "km/s")
print("en ellipse à l'apogée", vitesse_apogee, "km/s")

vitesse_p = vitesse_orbite_circulaire(r_p)
vitesse_a = vitesse_orbite_circulaire(r_a)

print("en circulaire au périgée", vitesse_p, "km/s")
print("en circulaire à l'apogée", vitesse_a, "km/s")

delta_vitesse_p = vitesse_perigee - vitesse_p
delta_vitesse_a = vitesse_apogee - vitesse_a

print(delta_vitesse_p)
print(delta_vitesse_a)


#%%

def vitesse_liberation(r, mu_terre = 398600.4418, r_T = 6378):
    return np.sqrt(2*mu_terre/(r + r_T))

print(vitesse_liberation(0), "km/s")
print(vitesse_liberation(200), "km/s")
print(vitesse_liberation(10000), "km/s")
print(vitesse_liberation(35786), "km/s")


