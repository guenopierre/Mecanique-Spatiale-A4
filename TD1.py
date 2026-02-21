# -*- coding: utf-8 -*-
"""
Script de calculs de mécanique orbitale.

Calculs des paramètres orbitaux, vitesses, périodes et altitudes
pour des orbites terrestres et héliocentriques.
"""

import numpy as np
import matplotlib.pyplot as plt


# ============================================================================
# CONSTANTES PHYSIQUES
# ============================================================================

G = 6.67408e-11  # Constante gravitationnelle universelle (m³ kg⁻¹ s⁻²)

# Masses célestes
M_SOLEIL = 2e30  # kg
M_TERRE = 6e24  # kg

# Paramètres de la Terre
R_TERRE = 6378  # km (rayon équatorial)
MU_TERRE = 398600.4418  # km³ s⁻² (paramètre gravitationnel)

# Conversions de temps
SECONDES_PAR_JOUR = 86400


# ============================================================================
# FONCTIONS DE CALCUL - DEMI-GRAND AXE
# ============================================================================

def demi_grand_axe_from_ray(r_apogee, r_perigee):
    """
    Calcul du demi-grand axe à partir des rayons d'apogée et de périgée.
    
    Args:
        r_apogee: Rayon à l'apogée (km)
        r_perigee: Rayon au périgée (km)
    
    Returns:
        Demi-grand axe (km)
    """
    return (r_apogee + r_perigee) / 2


def demi_grand_axe_from_apogee(r_apogee, excentricite):
    """
    Calcul du demi-grand axe à partir du rayon d'apogée et l'excentricité.
    
    Args:
        r_apogee: Rayon à l'apogée (km)
        excentricite: Excentricité de l'orbite
    
    Returns:
        Demi-grand axe (km)
    """
    return r_apogee / (1 + excentricite)


def demi_grand_axe_from_perigee(r_perigee, excentricite):
    """
    Calcul du demi-grand axe à partir du rayon de périgée et l'excentricité.
    
    Args:
        r_perigee: Rayon au périgée (km)
        excentricite: Excentricité de l'orbite
    
    Returns:
        Demi-grand axe (km)
    """
    return r_perigee / (1 - excentricite)


def demi_grand_axe_from_periode(periode, masse, G=G):
    """
    Calcul du demi-grand axe à partir de la période orbitale.
    
    Args:
        periode: Période orbitale (jours)
        masse: Masse du corps central (kg)
        G: Constante gravitationnelle (m³ kg⁻¹ s⁻²)
    
    Returns:
        Demi-grand axe (km)
    """
    periode_s = periode * SECONDES_PAR_JOUR
    a_m = ((periode_s**2 * G * masse) / (4 * np.pi**2))**(1/3)
    return a_m / 1e3


# ============================================================================
# FONCTIONS DE CALCUL - PÉRIODE ORBITALE
# ============================================================================

def periode_from_demi_grand_axe(demi_grand_axe, masse, G=G):
    """
    Calcul de la période orbitale à partir du demi-grand axe.
    
    Args:
        demi_grand_axe: Demi-grand axe (km)
        masse: Masse du corps central (kg)
        G: Constante gravitationnelle (m³ kg⁻¹ s⁻²)
    
    Returns:
        Période orbitale (jours)
    """
    a_m = demi_grand_axe * 1e3
    periode_s = 2 * np.pi * np.sqrt(a_m**3 / (G * masse))
    return periode_s / SECONDES_PAR_JOUR


def mouvement_moyen(periode):
    """
    Calcul du mouvement moyen d'une orbite.
    
    Args:
        periode: Période orbitale (secondes)
    
    Returns:
        Mouvement moyen (rad/s)
    """
    return (2 * np.pi) / periode


# ============================================================================
# FONCTIONS DE CALCUL - VITESSES
# ============================================================================

def vitesse_orbite_circulaire(rayon, mu=MU_TERRE):
    """
    Calcul de la vitesse pour une orbite circulaire.
    
    Args:
        rayon: Rayon orbital (km)
        mu: Paramètre gravitationnel (km³ s⁻²)
    
    Returns:
        Vitesse orbitale (km/s)
    """
    return np.sqrt(mu / rayon)


def vitesse_selon_rayon(rayon, demi_grand_axe, mu=MU_TERRE):
    """
    Calcul de la vitesse à un rayon donné pour une orbite elliptique.
    
    Args:
        rayon: Rayon au point considéré (km)
        demi_grand_axe: Demi-grand axe de l'orbite (km)
        mu: Paramètre gravitationnel (km³ s⁻²)
    
    Returns:
        Vitesse (km/s)
    """
    return np.sqrt(mu * ((2 / rayon) - (1 / demi_grand_axe)))


def vitesse_liberation(altitude, mu=MU_TERRE, rayon_terre=R_TERRE):
    """
    Calcul de la vitesse de libération à une altitude donnée.
    
    Args:
        altitude: Altitude au-dessus de la surface (km)
        mu: Paramètre gravitationnel (km³ s⁻²)
        rayon_terre: Rayon de la Terre (km)
    
    Returns:
        Vitesse de libération (km/s)
    """
    return np.sqrt(2 * mu / (altitude + rayon_terre))


# ============================================================================
# FONCTIONS DE CALCUL - ALTITUDES
# ============================================================================

def altitude_apogee(demi_grand_axe, excentricite, rayon_terre=R_TERRE):
    """
    Calcul de l'altitude à l'apogée.
    
    Args:
        demi_grand_axe: Demi-grand axe (km)
        excentricite: Excentricité de l'orbite
        rayon_terre: Rayon de la Terre (km)
    
    Returns:
        Altitude à l'apogée (km)
    """
    return (demi_grand_axe * (1 + excentricite)) - rayon_terre


def duree_orbite_circulaire(vitesse, rayon):
    """
    Calcul de la durée d'une orbite circulaire.
    
    Args:
        vitesse: Vitesse orbitale (km/s)
        rayon: Rayon orbital (km)
    
    Returns:
        Durée de l'orbite (secondes)
    """
    return 2 * np.pi * rayon / vitesse


# ============================================================================
# EXEMPLES ET TESTS
# ============================================================================

if __name__ == "__main__":
    
    # Exemple 1: Orbite basse terrestre
    print("\n--- Exemple 1: Orbite basse terrestre ---")
    r_t = R_TERRE * 1e3  # en mètres
    r_p = 620e3 + r_t
    e = 0.04
    a = demi_grand_axe_from_perigee(r_p / 1e3, e)
    print(f"Demi-grand axe: {a:.2f} km")
    print(f"Altitude à l'apogée: {altitude_apogee(a, e):.2f} km")
    
    # Exemple 2: Orbites héliocentriques
    print("\n--- Exemple 2: Orbites héliocentriques ---")
    a_mercure = demi_grand_axe_from_periode(87.969, M_SOLEIL)
    a_mars = demi_grand_axe_from_periode(686.886, M_SOLEIL)
    periode_asteroide = periode_from_demi_grand_axe(1426700000, M_SOLEIL)
    periode_lune = periode_from_demi_grand_axe(384399, M_TERRE)
    print(f"Demi-grand axe Mercure: {a_mercure:.0f} km")
    print(f"Demi-grand axe Mars: {a_mars:.0f} km")
    print(f"Période astéroïde: {periode_asteroide:.2f} jours")
    print(f"Période Lune: {periode_lune:.2f} jours")
    
    # Exemple 3: Orbite spécifique
    print("\n--- Exemple 3: Orbite à 620 km d'altitude ---")
    r_t = R_TERRE
    r_p = 620 + r_t
    e = 0.04
    a = demi_grand_axe_from_perigee(r_p, e)
    periode = periode_from_demi_grand_axe(a, M_TERRE)
    mm = mouvement_moyen(periode * SECONDES_PAR_JOUR)
    print(f"Période: {periode * 24 * 60:.2f} min")
    print(f"Mouvement moyen: {mm:.6f} rad/s")
    
    # Exemple 4: Vitesses en orbite circulaire
    print("\n--- Exemple 4: Vitesses en orbite circulaire ---")
    altitudes = [200, 1000, 10000, 35786]
    for alt in altitudes:
        v = vitesse_orbite_circulaire(alt + R_TERRE)
        durée = duree_orbite_circulaire(v, alt + R_TERRE) / SECONDES_PAR_JOUR
        print(f"À {alt} km: v = {v:.3f} km/s, durée = {durée:.4f} jours")
    
    # Exemple 5: Transfert d'orbite (Hohmann)
    print("\n--- Exemple 5: Transfert d'orbite Terre-GEO ---")
    r_t = R_TERRE
    alt_perigee = 300
    alt_apogee = 35786
    
    r_p = alt_perigee + r_t
    r_a = alt_apogee + r_t
    
    a_transfert = demi_grand_axe_from_ray(r_a, r_p)
    e_transfert = 1 - (r_p / a_transfert)
    
    print(f"Demi-grand axe de transfert: {a_transfert:.2f} km")
    print(f"Excentricité de transfert: {e_transfert:.6f}")
    
    v_p_ellipse = vitesse_selon_rayon(r_p, a_transfert)
    v_a_ellipse = vitesse_selon_rayon(r_a, a_transfert)
    v_p_circulaire = vitesse_orbite_circulaire(r_p)
    v_a_circulaire = vitesse_orbite_circulaire(r_a)
    
    delta_v_p = v_p_ellipse - v_p_circulaire
    delta_v_a = v_a_ellipse - v_a_circulaire
    
    print(f"Vitesse ellipse au périgée: {v_p_ellipse:.3f} km/s")
    print(f"Vitesse ellipse à l'apogée: {v_a_ellipse:.3f} km/s")
    print(f"Vitesse circulaire au périgée: {v_p_circulaire:.3f} km/s")
    print(f"Vitesse circulaire à l'apogée: {v_a_circulaire:.3f} km/s")
    print(f"Delta-v au périgée: {delta_v_p:.3f} km/s")
    print(f"Delta-v à l'apogée: {delta_v_a:.3f} km/s")
    
    # Exemple 6: Vitesses de libération
    print("\n--- Exemple 6: Vitesses de libération ---")
    altitudes = [0, 200, 10000, 35786]
    for alt in altitudes:
        v_lib = vitesse_liberation(alt)
        print(f"À {alt} km d'altitude: {v_lib:.3f} km/s")


