import numpy as np

def get_num_atoms_in_formula_unit(configuration):
    """Return the number ot atoms per formula unit.

    For unaries, it's the number in the primitive cell: 1 for SC, BCC, FCC; 2 for diamond.
    """
    mapping = {
        'XO': 2,
        'XO2': 3,
        'X2O': 3,
        'X2O3': 5,
        'XO3': 4,
        'X2O5': 7,
        'X/Diamond': 2,
        'X/SC': 1,
        'X/FCC': 1,
        'X/BCC': 1
    }

    try:
        return mapping[configuration]
    except KeyError:
        raise ValueError(f"Unknown number of atoms in formula unit for configuration '{configuration}'")


def get_volume_scaling_to_formula_unit(num_atoms_in_cell, element, configuration):
    """Return the scaling factor for extensive quantities (energies, volumes, ...).

    The volume per formula unit is the volume in the simulation cell (with `num_atoms_in_cell`)
    DIVIDED by the value returned by this function.
    """
    num_atoms_in_formula_unit = get_num_atoms_in_formula_unit(configuration)

    if element != "O":
        # For oxygen I might get a smaller cell;
        # the remaining logic is still correct (I divide by 0.5 => multiply by 2
        # if there is e.g. 1 atom in XO but the number of atoms in the formula unit is 2)
        assert num_atoms_in_cell % num_atoms_in_formula_unit == 0, (
            f"Weird! {num_atoms_in_cell} atoms in cell but "
            f"{num_atoms_in_formula_unit} in formula unit for {element}-{configuration}!")

    scaling = num_atoms_in_cell / num_atoms_in_formula_unit
    return scaling

def birch_murnaghan(V,E0,V0,B0,B01):
    """
    Return the energy for given volume (V - it can be a vector) according to
    the Birch Murnaghan function with parameters E0,V0,B0,B01.
    """
    r = (V0/V)**(2./3.)
    return (E0 +
            9./16. * B0 * V0 * (
            (r-1.)**3 * B01 + 
            (r-1.)**2 * (6. - 4.* r)))


def intE12sq(v0w,b0w,b1w,v0f,b0f,b1f,V1,V2):
    """
    Integral of (E1(V) - E2(V))**2 in dV evaluated between volume V1 and volume V2
    """
    F1 = antiderE12sq(v0w,b0w,b1w,v0f,b0f,b1f,V1)
    F2 = antiderE12sq(v0w,b0w,b1w,v0f,b0f,b1f,V2)
    integral = F2 - F1

    return integral

def antiderE12sq(v0w,b0w,b1w,v0f,b0f,b1f,V):
    """
    Antiderivative of (E1(V) - E2(V))**2 where E1(V) and E2(V) are birch murnaghan
    functions with different parameters
    """
    antider = (81*(\
            6*b0w*b0f*(-16 + 3*b1w)*(-16 + 3*b1f)*V*v0w*(v0w/V)**(2/3)*v0f*(v0f/V)**(2/3) - \
            2*b0w*b0f*(-14 + 3*b1w)*(-16 + 3*b1f)*V*v0w*(v0w/V)**(4/3)*v0f*(v0f/V)**(2/3) - \
            2*b0w*b0f*(-16 + 3*b1w)*(-14 + 3*b1f)*V*v0w*(v0w/V)**(2/3)*v0f*(v0f/V)**(4/3) + \
            (6*b0w*b0f*(-14 + 3*b1w)*(-14 + 3*b1f)*V*v0w*(v0w/V)**(4/3)*v0f*(v0f/V)**(4/3))/5. + \
            V*(b0w*(-6 + b1w)*v0w - b0f*(-6 + b1f)*v0f)**2 - \
            (b0w*(-4 + b1w)*v0w**3 - b0f*(-4 + b1f)*v0f**3)**2/(3.*V**3) + \
            (3*b0f*(v0f/V)**(7/3)*(-2*b0w*(-14 + 3*b1f)*v0w*(-7*(-6 + b1w)*V**2 + \
            (-4 + b1w)*v0w**2) \
            - 7*b0f*(424 + 5*b1f*(-32 + 3*b1f))*V**2*\
            v0f + 2*b0f*(-4 + b1f)*(-14 + 3*b1f)*\
            v0f**3))/7. - \
            (3*b0f*(v0f/V)**(5/3)*\
            (-2*b0w*(-16 + 3*b1f)*v0w*\
            (5*(-6 + b1w)*V**2 + (-4 + b1w)*v0w**2) \
            + 10*b0f*(-6 + b1f)*(-16 + 3*b1f)*V**2*\
            v0f + b0f*(324 + 5*b1f*(-28 + 3*b1f))*\
            v0f**3))/5. +\
            (4*b0w**2*(124 + 5*(-10 + b1w)*b1w)*v0w**4 - \
            2*b0w*b0f*(-4 + b1w)*(-6 + b1f)*v0w**3*\
            v0f - 2*b0w*b0f*(-6 + b1w)*(-4 + b1f)*v0w*\
            v0f**3 + 4*b0f**2*(124 + 5*(-10 + b1f)*b1f)*v0f**4)/V + \
            (3*b0w*(v0w/V)**(7/3)*\
            (-7*b0w*(424 + 5*b1w*(-32 + 3*b1w))*V**2*\
            v0w + 2*b0w*(-4 + b1w)*(-14 + 3*b1w)*v0w**3 + \
            2*b0f*(-14 + 3*b1w)*v0f*\
            (7*(-6 + b1f)*V**2 - (-4 + b1f)*v0f**2)))/7. - \
            (3*b0w*(v0w/V)**(5/3)*\
            (10*b0w*(-6 + b1w)*(-16 + 3*b1w)*V**2*v0w + \
            b0w*(324 + 5*b1w*(-28 + 3*b1w))*v0w**3 - \
            2*b0f*(-16 + 3*b1w)*v0f*\
            (5*(-6 + b1f)*V**2 + (-4 + b1f)*v0f**2)))/5.))/256.

    return antider

def intEdV(V0,B0,B0pr,V1,V2):
    """
    integral of E(V) in dV evaluated between volumes V1 and V2
    """
    F1 = antiderE(V0,B0,B0pr,V1)
    F2 = antiderE(V0,B0,B0pr,V2)
    integral = F2 - F1

    return integral

def antiderE(V0,B0,B0pr,V):
    """
    antiderivative of the Birch Murnaghan E(V)
    """
    antider = (9*B0*V0*(-((-6 + B0pr)*V) - ((-4 + B0pr)*V0**2)/V + \
            3*(-14 + 3*B0pr)*V0*(V0/V)**(1/3) + \
            3*(-16 + 3*B0pr)*V*(V0/V)**(2/3)))/16

    return antider

def intE2dV(V0,B0,B0pr,V1,V2):
    """
    Integral of E**2(V) in dV evaluated between volume V1 and volume V2
    """
    F1 = antiderE2(V0,B0,B0pr,V1)
    F2 = antiderE2(V0,B0,B0pr,V2)
    integral = F2 - F1

    return integral

def antiderE2(V0,B0,B0pr,V):
    """
    Antiderivative of the Birch Murnaghan squared (E**2(V))
    """
    antider = (81*B0**2*V0**2*((-6 + B0pr)**2*V + \
            (4*(124 + 5*(-10 + B0pr)*B0pr)*V0**2)/V - \
            ((-4 + B0pr)**2*V0**4)/(3.*V**3) - \
            (3*(V0/V)**(2/3)* \
            (10*(-6 + B0pr)*(-16 + 3*B0pr)*V**2 + \
            (324 + 5*B0pr*(-28 + 3*B0pr))*V0**2))/(5.*V) \
            + (V0/V)**(1/3)* \
            (-3*(424 + 5*B0pr*(-32 + 3*B0pr))*V0 + \
            (6*(-4 + B0pr)*(-14 + 3*B0pr)*V0**3)/(7.*V**2))) \
            )/256.

    return antider

def delta(v0w, b0w, b1w, v0f, b0f, b1f, prefact, weight_b0, weight_b1):
    """
    Calculate the Delta value, function copied from the official DeltaTest repository.
    I don't understand what it does, but it works. Result is in meV
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    """

    Vi = 0.94 * (v0w + v0f) / 2.
    Vf = 1.06 * (v0w + v0f) / 2.

    a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
    a2f = 9. * v0f**(7. / 3.) * b0f / 16. * (14. - 3. * b1f)
    a1f = 9. * v0f**(5. / 3.) * b0f / 16. * (3. * b1f - 16.)
    a0f = 9. * v0f * b0f / 16. * (6. - b1f)

    a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
    a2w = 9. * v0w**(7. / 3.) * b0w / 16. * (14. - 3. * b1w)
    a1w = 9. * v0w**(5. / 3.) * b0w / 16. * (3. * b1w - 16.)
    a0w = 9. * v0w * b0w / 16. * (6. - b1w)

    x = [0, 0, 0, 0, 0, 0, 0]

    x[0] = (a0f - a0w)**2
    x[1] = 6. * (a1f - a1w) * (a0f - a0w)
    x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
    x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
    x[4] = -3. / 5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
    x[5] = -6. / 7. * (a3f - a3w) * (a2f - a2w)
    x[6] = -1. / 3. * (a3f - a3w)**2.

    y = [0, 0, 0, 0, 0, 0, 0]

    y[0] = (a0f + a0w)**2 / 4.
    y[1] = 3. * (a1f + a1w) * (a0f + a0w) / 2.
    y[2] = -3. * (2. * (a2f + a2w) * (a0f + a0w) + (a1f + a1w)**2.) / 4.
    y[3] = -(a3f + a3w) * (a0f + a0w) / 2. - (a2f + a2w) * (a1f + a1w) / 2.
    y[4] = -3. / 20. * (2. * (a3f + a3w) * (a1f + a1w) + (a2f + a2w)**2.)
    y[5] = -3. / 14. * (a3f + a3w) * (a2f + a2w)
    y[6] = -1. / 12. * (a3f + a3w)**2.

    Fi = np.zeros_like(Vi)
    Ff = np.zeros_like(Vf)

    Gi = np.zeros_like(Vi)
    Gf = np.zeros_like(Vf)

    for n in range(7):
        Fi = Fi + x[n] * Vi**(-(2. * n - 3.) / 3.)
        Ff = Ff + x[n] * Vf**(-(2. * n - 3.) / 3.)

        Gi = Gi + y[n] * Vi**(-(2. * n - 3.) / 3.)
        Gf = Gf + y[n] * Vf**(-(2. * n - 3.) / 3.)

    Delta = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi))
    #Deltarel = 100. * np.sqrt((Ff - Fi) / (Gf - Gi))
    #vref = 30.
    #bref = 100. * 10.**9. / 1.602176565e-19 / 10.**30. #100 GPa in ev_ang3
    #Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
    #    / (v0w + v0f) / (b0w + b0f) * 4. * vref * bref

    return Delta  #, Deltarel, Delta1



def epsilon(v0w, b0w, b1w, v0f, b0f, b1f, prefact, weight_b0, weight_b1):
    """
    Calculate alternative Delta2 based on 2 EOS fits
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    """

    # volume range
    Vi = 0.94 * (v0w + v0f) / 2.
    Vf = 1.06 * (v0w + v0f) / 2.
    deltaV = Vf - Vi

    intdiff2 = intE12sq(v0w,b0w,b1w,v0f,b0f,b1f,Vi,Vf)
    Eavg1 = intEdV(v0w,b0w,b1w,Vi,Vf)/deltaV
    Eavg2 = intEdV(v0w,b0w,b1w,Vi,Vf)/deltaV
    int3 = intE2dV(v0w,b0w,b1w,Vi,Vf) - \
            2*Eavg1*intEdV(v0w,b0w,b1w,Vi,Vf) + \
            deltaV*Eavg1**2 # integrate (ene - mean(ene))**2
    int4 = intE2dV(v0f,b0f,b1f,Vi,Vf) - \
            2*Eavg2*intEdV(v0f,b0f,b1f,Vi,Vf) + \
            deltaV*Eavg2**2
    eps2 = intdiff2/np.sqrt(int3*int4)

    #We saw a case when, for numerical error, intdiff2 was negative
    #(about -1*10^{-13}). For this reason, we add a safty check.
    if eps2 < 0.0:
        print(f"eps2 = {eps2}, negative due to numerical error probably, we take the absolute value")
        eps2 = abs(eps2)
    
    return np.sqrt(eps2)*prefact


def V0_rel_diff(v0w, b0w, b1w, v0f, b0f, b1f, prefact, weight_b0, weight_b1):
    """
    Returns the relative difference in the volumes.
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though several inputs are useless here.
    """
    return prefact*2*(v0w-v0f)/(v0w+v0f)


def B0_rel_diff(v0w, b0w, b1w, v0f, b0f, b1f, prefact, weight_b0, weight_b1):
    """
    Returns the relative difference in the bulk modulus.
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though several inputs are useless here.
    """
    return prefact*2*(b0w-b0f)/(b0w+b0f)

def B1_rel_diff(v0w, b0w, b1w, v0f, b0f, b1f, prefact, weight_b0, weight_b1):
    """
    Returns the reletive difference in the derivative of the bulk modulus.
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    Even though several inputs are useless here.
    """
    return prefact*2*(b1w-b1f)/(b1w+b1f)

def rel_errors_vec_length_unsquared(v0w, b0w, b1w, v0f, b0f, b1f, prefact, weight_b0, weight_b1):
    """
    Returns the length of the vector formed by the relative error of V0, B0, B1.
    Note that this is NOT the nu, but a version in which the weights are outside the square.
    This version should *not* be used when using the weights discussed in the paper for nu.
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    """
    V0err =  2*(v0w-v0f)/(v0w+v0f)
    B0err =  2*(b0w-b0f)/(b0w+b0f)
    B1err =  2*(b1w-b1f)/(b1w+b1f)
    leng = np.sqrt(V0err**2+weight_b0*(B0err)**2+weight_b1*(B1err)**2)
    return leng*prefact

def nu(v0w, b0w, b1w, v0f, b0f, b1f, prefact, weight_b0, weight_b1):
    """
    Returns the "nu" measure defined with appropriate weights applied to the relative error of V0, B0, B1.
    (Note that the weights are directy multiplied to the relative errors, i.e., they are inside the square).
    THE SIGNATURE OF THIS FUNCTION HAS BEEN CHOSEN TO MATCH THE ONE OF ALL THE OTHER FUNCTIONS
    RETURNING A QUANTITY THAT IS USEFUL FOR COMPARISON, THIS SIMPLIFIES THE CODE LATER.
    """
    V0err =  2*(v0w-v0f)/(v0w+v0f)
    B0err =  2*(b0w-b0f)/(b0w+b0f)
    B1err =  2*(b1w-b1f)/(b1w+b1f)
    leng = np.sqrt(V0err**2+(weight_b0*B0err)**2+(weight_b1*B1err)**2)
    return leng*prefact
