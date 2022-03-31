last_updated = "220331"

# load libraries
import sys, os
import re
from datetime import datetime
import pandas as pd
import pkg_resources

rt = 3 # round output to how many decimal places?


# GB: function to write output to a csv
def write_output(filename, CATION, LIGAND, nth_complex, G, H, S, CP, V, a1, a2, a3, a4, c1, c2, wcon, Z):
    
    cat_nocharge = re.sub("\+.*$", "", CATION)
    lig_nocharge = re.sub("\-.*$", "", LIGAND)
    
    if Z > 1:
        this_charge = "+" + str(Z)
    elif Z == 1:
        this_charge = "+"
    elif Z == 0:
        this_charge = ""
    elif Z == -1:
        this_charge = "-"
    else:
        this_charge = str(Z)
    
    if nth_complex == 1:
        ligand_subscript = ""
    else:
        ligand_subscript = str(nth_complex)
    
    this_date = datetime.today().strftime('%Y%m%d') 
    
    data={
        "name": [cat_nocharge + "(" + lig_nocharge + ")" + ligand_subscript + this_charge],
        "abbrv": [cat_nocharge + "(" + lig_nocharge + ")" + ligand_subscript + this_charge],
        "formula": [cat_nocharge + "(" + lig_nocharge + ")" + ligand_subscript + this_charge],
        "state": ["aq"],
        "ref1": ["autocomplicator"],
        "ref2": ["1" + " " + CATION + "," + str(nth_complex) + " " + LIGAND],  # future: "1" can be replaced with number of cations
        "date": [this_date],
        "E_units": ["cal"],
        "G": [G],
        "H": [H],
        "S": [S],
        "Cp": [CP],
        "V": [V],
        "a1.a": [a1],
        "a2.b": [a2],
        "a3.c": [a3],
        "a4.d": [a4],
        "c1.e": [c1],
        "c2.f": [c2],
        "omega.lambda": [wcon],
        "z.T":[Z],
    }
    df = pd.DataFrame(data)
    
    try:
        file_exists = os.path.isfile(filename+'.csv')
        if file_exists:
            with open(filename+'.csv', 'a') as f:
                df.to_csv(f, header=False, index=False)
        else:
            df.to_csv(filename+'.csv', index=False)
            
    except:
        import warnings
        warnings.warn("The file '" + filename + ".csv' could not be written! Make sure this is a valid filename.")
    
    return df


# function to count significant digits
# adapted from https://stackoverflow.com/questions/8142676/python-counting-significant-digits
# function to count significant digits
# adapted from https://stackoverflow.com/questions/8142676/python-counting-significant-digits
def find_sigfigs(x):
    '''Returns the number of significant digits in a number. This takes into account
       strings formatted in 1.23e+3 format and even strings such as 123.450'''
    # change all the 'E' to 'e'
    x = x.lower()
    x = re.sub("^-", "", x)# remove negative sign
    if ('e' in x):
        # return the length of the numbers before the 'e'
        myStr = x.split('e')
        if ("." in myStr[0]):
            return len( myStr[0] ) - 1 # to compenstate for the decimal point
        else:
            return len( myStr[0] )
    else:
        # put it in e format and return the result of that
        ### NOTE: because of the 8 below, it may do crazy things when it parses 9 sigfigs
        n = ('%.*e' %(8, float(x))).split('e')
        # remove and count the number of removed user added zeroes. (these are sig figs)
        if '.' in x:
            s = x.replace('.', '')
            #number of zeroes to add back in
            l = len(s) - len(s.rstrip('0'))
            #strip off the python added zeroes and add back in the ones the user added
            n[0] = n[0].rstrip('0') + ''.join(['0' for num in range(l)])
        else:
            #the user had no trailing zeroes so just strip them all
            n[0] = n[0].rstrip('0')
        #pass it back to the beginning to be parsed
    return find_sigfigs('e'.join(n))


def complicate(cation_override, ligand_override, beta, sest, out_name, 
               sigfigs=False, metal_ref_path=None, ligand_ref_path=None):
    
    BETA1, BETA2, BETA3, BETA4 = beta
    sest1, sest2, sest3, sest4 = sest
    
    # look up metal cation values from reference sheet
    # open cation reference sheet
    try:
        if metal_ref_path == None:
            metal_ref_path = "Metal_HKF_Parameters.txt"
            metal_ref = pkg_resources.resource_string(__name__, metal_ref_path).decode("utf-8")
            d = "\n"
            metal_ref_split = metal_ref.split(d)
            metal_ref = [line+d for line in metal_ref_split]
        else:
            metal_ref = open(metal_ref_path, "r").readlines()
    except:
        print("Error in cation override. Could not find or open", metal_ref_path)
        sys.exit(0)

    # look up user-supplied cation in reference sheet
    try:
        to_exec = metal_ref[metal_ref.index('CATION = "' + cation_override + '"\n') + 1]
        to_exec = to_exec.replace("\n", "")
    except:
        print("Error in cation override. Could not find entry for", cation_override, "in reference sheet.")
        sys.exit(0)

    # get significant digits for each value
    try:
        cation_vals = re.sub(" #.*$", "", to_exec) # removes comment
        cation_vals = cation_vals.split(" = ")[1]
        cation_vals = cation_vals.split(", ")
        cation_vals = [re.sub("\.$", "", val) for val in cation_vals] # removes terminating periods TODO: fix this in cation database
        GC_sf = find_sigfigs(cation_vals[1])
        HC_sf = find_sigfigs(cation_vals[2])
        SC_sf = find_sigfigs(cation_vals[3])
        CPC_sf = find_sigfigs(cation_vals[4])
        VC_sf = find_sigfigs(cation_vals[5])
    except:
        print("Error in cation override. Could not count significant digits for values supplied for", cation_override)
        print("Cation values:", cation_vals)
        sys.exit(0)

    # evaluate entry as code
    try:
        exec('CATION = "' + cation_override + '"')
        exec(to_exec)
    except:
        print("Error in cation override. Could not evaluate", to_exec, "as Python code.")
        sys.exit(0)

        
        
    # look up ligand values from reference sheet
    # open ligand reference sheet
    try:
        if ligand_ref_path == None:
            ligand_ref_path = "Monovalent_ligand_HKF_Parameters.txt"
            ligand_ref = pkg_resources.resource_string(__name__, ligand_ref_path).decode("utf-8")
            d = "\n"
            ligand_ref_split = ligand_ref.split(d)
            ligand_ref = [line+d for line in ligand_ref_split]
        else:
            ligand_ref = open(ligand_ref_path, "r").readlines()
    except:
        print("Error in ligand override. Could not find or open", ligand_ref_path)
        sys.exit(0)

    # look up user-supplied ligand in reference sheet
    try:
        to_exec = ligand_ref[ligand_ref.index('LIGAND = "' + ligand_override + '"\n') + 1]
        to_exec = to_exec.replace("\n", "")
    except:
        print("Error in ligand override. Could not find entry for", ligand_override, "in reference sheet.")
        sys.exit(0)

    # get significant digits for each value
    try:
        ligand_vals = re.sub(" #.*$", "", to_exec) # removes comment
        ligand_vals = ligand_vals.split(" = ")[1]
        ligand_vals = ligand_vals.split(", ")
        ligand_vals = [re.sub("\.$", "", val) for val in ligand_vals] # removes terminating periods TODO: fix this in cation database
        GL_sf = find_sigfigs(ligand_vals[1])
        HL_sf = find_sigfigs(ligand_vals[2])
        SL_sf = find_sigfigs(ligand_vals[3])
        CPL_sf = find_sigfigs(ligand_vals[4])
        VL_sf = find_sigfigs(ligand_vals[5])
    except:
        print("Error in ligand override. Could not count significant digits for values supplied for", ligand_override)
        print("Ligand values:", ligand_vals)
        sys.exit(0)

    # evaluate entry as code
    try:
        exec('LIGAND = "' + ligand_override + '"')
        exec(to_exec)
    except:
        print("Error in ligand override. Could not evaluate", to_exec, "as Python code.")
        sys.exit(0)

    CATION = cation_override
    ZC, GC, HC, SC, CPC, VC = cation_vals[0], cation_vals[1], cation_vals[2], cation_vals[3], cation_vals[4], cation_vals[5]
    ZC, GC, HC, SC, CPC, VC = int(float(ZC)), float(GC), float(HC), float(SC), float(CPC), float(VC)       
    
    LIGAND = ligand_override
    ZL, GL, HL, SL, CPL, VL = ligand_vals[0], ligand_vals[1], ligand_vals[2], ligand_vals[3], ligand_vals[4], ligand_vals[5]
    ZL, GL, HL, SL, CPL, VL = int(float(ZL)), float(GL), float(HL), float(SL), float(CPL), float(VL)

    # Calculations for the first complex
    Z = ZC + ZL
    DELGR1 = (2.30259)*(1.98719)*(298.15)*BETA1
    G1 = DELGR1 + GC + GL

    # Now the entropy predictor starts for the first complex

    # if sest1 is available, use it. If not, predict it.
    if sest1 != 0:
        S1 = sest1
        DELSR1 = S1 - SC - SL
        DELS1 = DELSR1

    else:
        AZ= 0.016241*Z - 0.000479
        AZP= -0.360972*Z + 0.3209
        BZ= 0.32102*Z  - 0.05996
        BZP= 8.2198*Z - 1.557

        ALPHA = AZ*(SL+(-5.0*ZL)) + AZP
        BETA = BZ*(SL+(-5.0*ZL)) + BZP

        DELS1 = ALPHA*(SC+(-5.0*ZC)) + BETA
        DELSR1= DELS1

        S1 = DELSR1 + SC + SL

    Z1 = Z

    DELHR1 = DELGR1 + (298.15*DELSR1)
    H1 = DELHR1 + HC + (1.0*HL)

    # Here the Cp predictor starts for the first complex. Modified 8 Aug 1992 after latest revision Modified 17 March 1992 to incorporate changes from last July 1991
    dz = 0.856*CPL - 2.1 + 45.3*ZC
    DELCP1 = 1.25*CPC + dz
    CP1 = DELCP1 + CPC + CPL

    # Here the Volume predictor starts for the first complex
    DELVR1  = 0.11419*VC + 8.9432
    V1 = DELVR1 + VC + VL
    Z1 = Z

    # Calculations for the second complex
    Z = Z + ZL

    DELGR2 = 2.30259*1.98719*298.15*BETA2
    G2 = DELGR2 + GC + (2.0*GL)

    # if sest2 is available, use it. If not, predict it.
    if sest2 != 0:
        S2=sest2
        DELSR2 = S2 - SC - (2.0*SL)
        DELS2 = DELSR2 - DELS1
    else:
        AZ= 0.016241*Z - 0.000479
        AZP= -0.360972*Z + 0.3209
        BZ= 0.32102*Z  - 0.05996
        BZP= 8.2198*Z - 1.557

        ALPHA = AZ*(SL+(-5.0*ZL)) + AZP
        BETA = BZ*(SL+(-5.0*ZL)) + BZP

        DELS2 = ALPHA*(S1+(-5.0*(ZC+ZL))) + BETA
        DELSR2= DELS1 + DELS2

        S2 = DELSR2 + SC + (2.0*SL)

    Z2 = Z

    DELHR2 = DELGR2 + (298.15*DELSR2)
    H2 = DELHR2 + HC + (2.0*HL)

    # CHANGED 17 MARCH 1992!!!

    gz = 0.89*CPC + 0.72*CPL + 16.3

    DELCP2 = DELCP1 + gz

    CP2 = DELCP2 + CP1 + CPL

    Z2 = Z

    DELVR2  = 0.11419*V1 + 8.9432
    V2 = DELVR2 + DELVR1 + VC + (2.0*VL)
    Z2 = Z

    # Calculations for the third complex

    Z = Z + ZL  
    DELGR3 = 2.30259*1.98719*298.15*BETA3
    G3 = DELGR3 + GC + (3.0*GL)

    # if sest3 is available, use it. If not, predict it.
    if sest3 != 0:
        S3=sest3
        DELSR3 = S3 - SC - (3.0*SL)
        DELS3 = DELSR3 - DELS1 - DELS2
    else:
        AZ= 0.016241*Z - 0.000479
        AZP= -0.360972*Z + 0.3209
        BZ= 0.32102*Z  - 0.05996
        BZP= 8.2198*Z - 1.557

        ALPHA = AZ*(SL+(-5.0*ZL)) + AZP
        BETA = BZ*(SL+(-5.0*ZL)) + BZP

        DELS3 = ALPHA*(S2+(-5.0*(ZC+(2*ZL)))) + BETA
        DELSR3= DELS1 + DELS2 + DELS3

        S3 = DELSR3 + SC + (3.0*SL)

    Z3 = Z

    DELHR3 = DELGR3 + (298.15*DELSR3)
    H3 = DELHR3 + HC + (3.0*HL)

    # CHANGED 17 MARCH 1992!!!

    gz = 0.89*CPC + 0.72*CPL + 16.3

    DELCP3 = DELCP2 + gz

    CP3 = DELCP3 + CP2 + CPL

    Z3 = Z

    DELVR3  = 0.11419*V2 + 8.9432
    V3 = DELVR3 + DELVR2 + DELVR1 + VC + (3.0*VL)
    Z3 = Z

    # Calculations for the fourth complex

    Z = Z + ZL

    DELGR4 = 2.30259*1.98719*298.15*BETA4
    G4 = DELGR4 + GC + (4.0*GL)

    # if sest4 is available, use it. If not, predict it.
    if sest4 != 0:
        S4=sest4
        DELSR4 = S4 - SC - (4.0*SL)
        DELS4 = DELSR4 - DELS1 - DELS2 - DELS3
        Z4 = Z
    else:
        AZ= 0.016241*Z - 0.000479
        AZP= -0.360972*Z + 0.3209
        BZ= 0.32102*Z  - 0.05996
        BZP= 8.2198*Z - 1.557

        ALPHA = AZ*(SL+(-5.0*ZL)) + AZP
        BETA = BZ*(SL+(-5.0*ZL)) + BZP

        DELS4 = ALPHA*(S3+(-5.0*(ZC+(3*ZL)))) + BETA
        DELSR4= DELS1 + DELS2 + DELS3 + DELS4

        S4 = DELSR4 + SC + (4.0*SL)

    Z4 = Z
    DELHR4 = DELGR4 + (298.15*DELSR4)
    H4 = DELHR4 + HC + (4.0*HL)

    # CHANGED 17 MARCH 1992!!!

    gz = 0.89*CPC + 0.72*CPL + 16.3

    DELCP4 = DELCP3 + gz

    CP4 = DELCP4 + CP3 + CPL

    Z4 = Z

    DELVR4  = 0.11419*V3 + 8.9432
    V4 = DELVR4 + DELVR3 + DELVR2 + DELVR1 + VC + (4.0*VL)
    Z4 = Z

    # Now the partial molal properties are used in PARCOR to estimate equation of state coefficients for each complex.

    def calc_params(Z, G, H, S, CP, V):
        # The following are the Born functions at 25 C, and various other constants.

        AK = 0.0
        RX = 0.0
        sigma = 0.0
        POLAR = 2.0

        Q= 5.903E-07
        X= -3.090E-07
        Y= -5.802E-05
        xN= -1.820E-10
        conv=41.8393
        tr=298.15
        theta = 228.0
        pfunk=2601.0
        eta=1.66027E+05
        wabsh = 0.5387E+05

        # Calculation of GAMMA term (Helgeson & Kirkham 1976, page 153)

        if Z > 0:
            gamma=0.94
        elif Z < 0:
            gamma=0.0

        if Z != 0:
            if abs(Z) == 1:
                alphaz=72
            elif abs(Z) == 2:
                alphaz = 141
            elif abs(Z) == 3:
                alphaz = 211
            elif abs(Z) == 4:
                alphaz = 286
            else:
                print("alphaz cannot be calculated!")
                sys.exit(0)

            if S != 0.0:
                if RX != 0.0:
                    RE = RX+R*gamma
                else:
                    ra=(Z**2 *(eta * Y - 100.)/(S - alphaz))
                    ire=int(100.*ra +.5)
                    re=ire/100.
            else:
                if RX != 0.0:
                    re = rx+z*gamma
                    S =(z**2)*(eta*Y - 100.) / re + alphaz
                else:
                    print('IF Z NE 0.0, EITHER S OR RX MUST BE GIVEN!')
                    sys.exit(0)


            wabs = (eta*Z**2)/re
            wcon = wabs-Z*wabsh
        else:
            if POLAR == 2:
                wcon= -0.03*100000
                iwcon = int(wcon + .5)
                wcon = iwcon

        # Calculation of solvation contributions to the standard partial molal properties
        vs=-wcon*Q*conv
        cps=wcon*tr*X
        ss=wcon*Y

        # Calculation of the non-solvation contributions to the standard partial molal properties
        vn=V-vs
        cpn=CP-cps
        sn=S-ss

        # THE CORRELATIONS FOR E.O.S. PARAMETERS BEGIN !
        if AK != 0.0:
            aks=wcon*xN*conv
            akn=AK-aks
            if sigma != 0.0:
                a2=17.19E+4*akn+421.1
                # rounding
                ia2 = int(a2*100 + .5)
                a2 = ia2/100

                a1=sigma/conv-(a2/pfunk)
                # rounding
                ia1 = int(a1*100000 + .5)
                a1 = ia1/100000

            else:
                a2=17.19E+4*akn+421.1
                # rounding
                ia2 = int(a2*100 + .5)
                a2 = ia2/100

                sigma=1.11*vn+1.8
                # rounding
                isigma = int(sigma*100 + .5)
                sigma = isigma/100

                a1=sigma/conv-(a2/pfunk)
                # rounding
                ia1 = int(a1*100000 + .5)
                a1 = ia1/100000
        else:
            if sigma != 0.0:
                a1=1.3684E-2*vn+0.1765
                # rounding
                ia1 = int(a1*100000. + .5)
                a1 = ia1/100000.

                a2=((sigma/conv)-a1)*pfunk
                # rounding
                ia2 = int(a2*100. + .5)
                a2 = ia2/100.

            else:
                sigma=1.11*vn+1.8
                # rounding
                isigma = int(sigma*100. + .5)
                sigma = isigma/100.

                a1=1.3684E-2*vn+0.1765
                # rounding
                ia1 = int(a1*100000 + .5)
                a1 = ia1/100000

                a2=((sigma/conv)-a1)*pfunk
                # rounding
                ia2 = int(a2*100 + .5)
                a2 = ia2/100

        a4 = -4.134*a2-27790.0
        # rounding
        ia4 = int(a4 +.5)
        a4 = ia4

        a3=((vn/conv)-a1-a2/pfunk)*(tr-theta)-(a4/pfunk)
        # rounding
        ia3 = int(a3*10000 + .5)
        a3 = ia3/10000

        c2= .2037*CP-3.0346
        # rounding
        ic2=int(c2*10000 + .5)
        c2 = ic2

        c1=cpn-c2/(tr-theta)**2
        # rounding
        ic1 = int(c1*10000 + .5)
        c1= ic1/10000

        akn=(a2+a4/(tr-theta))*(1/pfunk**2)
        aks=wcon*conv*xN
        ak=akn+aks

        return G, H, S, a1, a2, a3, a4, c1, c2, wcon

        # THAT'S ALL FOLKS !

    # The following section has been altered at GEOPIG so that there will be no output if Beta = 0.0

    G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params(Z1, G1, H1, S1, CP1, V1)

    def sf(val, sf):
        val = "{0:.{prec}g}".format(val, prec=sf)
        if "e" in val and abs(float(val)) >= 1:
            val = str(int(float(val)))
        elif "e" not in val and abs(float(val)) < 1:
            val = val + "0"*(sf - find_sigfigs(val))
        elif "e" in val and abs(float(val)) < 1:
            val = "{0:.{prec}f}".format(float(val), prec=20) # add arbitrarily high trailing zeros
            val = re.sub("0+$", "", val) # trim off all trailing zeros
            val = val + "0"*(sf - find_sigfigs(val)) # add back correct number of trailing zeros
        return val
    
    if sigfigs:
        G_out = sf(G, min(GC_sf, GL_sf))
        H_out = sf(H, min(HC_sf, HL_sf))
        S_out = sf(S, min(SC_sf, SL_sf))
        CP_out = sf(CP1, min(CPC_sf, CPL_sf))
        V_out = sf(V1, min(VC_sf, VL_sf))
        a1_out = sf(a1*10, min(VC_sf, VL_sf))
        a2_out = sf(a2/100, min(VC_sf, VL_sf))
        a3_out = sf(a3, min(VC_sf, VL_sf))
        a4_out = sf(a4/10000, min(VC_sf, VL_sf))
        c1_out = sf(c1, min(CPC_sf, CPL_sf))
        c2_out = sf(c2/10000, min(CPC_sf, CPL_sf))
        wcon_out = sf(wcon/100000, min(SC_sf, SL_sf))
    else:
        G_out = round(G, rt)
        H_out = round(H, rt)
        S_out = round(S, rt)
        CP_out = round(CP1, rt)
        V_out = round(V1, rt)
        a1_out = round(a1*10, rt)
        a2_out = round(a2/100, rt)
        a3_out = round(a3, rt)
        a4_out = round(a4/10000, rt)
        c1_out = round(c1, rt)
        c2_out = round(c2/10000, rt)
        wcon_out = round(wcon/100000, rt)
    Z_out = Z1


    write_output(filename=out_name, CATION=CATION, LIGAND=LIGAND, nth_complex=1, G=G_out, H=H_out, S=S_out,
                 CP=CP_out, V=V_out, a1=a1_out, a2=a2_out, a3=a3_out, a4=a4_out, c1=c1_out, c2=c2_out,
                 wcon=wcon_out, Z=Z_out)

    if BETA2 != 0:

        G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params(Z2, G2, H2, S2, CP2, V2)
        
        if sigfigs:
            G_out = sf(G, min(GC_sf, GL_sf))
            H_out = sf(H, min(HC_sf, HL_sf))
            S_out = sf(S, min(SC_sf, SL_sf))
            CP_out = sf(CP1, min(CPC_sf, CPL_sf))
            V_out = sf(V1, min(VC_sf, VL_sf))
            a1_out = sf(a1*10, min(VC_sf, VL_sf))
            a2_out = sf(a2/100, min(VC_sf, VL_sf))
            a3_out = sf(a3, min(VC_sf, VL_sf))
            a4_out = sf(a4/10000, min(VC_sf, VL_sf))
            c1_out = sf(c1, min(CPC_sf, CPL_sf))
            c2_out = sf(c2/10000, min(CPC_sf, CPL_sf))
            wcon_out = sf(wcon/100000, min(SC_sf, SL_sf))
        else:
            G_out = round(G, rt)
            H_out = round(H, rt)
            S_out = round(S, rt)
            CP_out = round(CP2, rt)
            V_out = round(V2, rt)
            a1_out = round(a1*10, rt)
            a2_out = round(a2/100, rt)
            a3_out = round(a3, rt)
            a4_out = round(a4/10000, rt)
            c1_out = round(c1, rt)
            c2_out = round(c2/10000, rt)
            wcon_out = round(wcon/100000, rt)
        Z_out = Z2

        write_output(filename=out_name, CATION=CATION,
                     LIGAND=LIGAND, nth_complex=2,
                     G=G_out, H=H_out, S=S_out,
                     CP=CP_out, V=V_out, a1=a1_out,
                     a2=a2_out, a3=a3_out, a4=a4_out,
                     c1=c1_out, c2=c2_out,
                     wcon=wcon_out, Z=Z_out)

    if BETA3 != 0:

        G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params(Z3, G3, H3, S3, CP3, V3)
    
        if sigfigs:
            G_out = sf(G, min(GC_sf, GL_sf))
            H_out = sf(H, min(HC_sf, HL_sf))
            S_out = sf(S, min(SC_sf, SL_sf))
            CP_out = sf(CP1, min(CPC_sf, CPL_sf))
            V_out = sf(V1, min(VC_sf, VL_sf))
            a1_out = sf(a1*10, min(VC_sf, VL_sf))
            a2_out = sf(a2/100, min(VC_sf, VL_sf))
            a3_out = sf(a3, min(VC_sf, VL_sf))
            a4_out = sf(a4/10000, min(VC_sf, VL_sf))
            c1_out = sf(c1, min(CPC_sf, CPL_sf))
            c2_out = sf(c2/10000, min(CPC_sf, CPL_sf))
            wcon_out = sf(wcon/100000, min(SC_sf, SL_sf))
        else:
            G_out = round(G, rt)
            H_out = round(H, rt)
            S_out = round(S, rt)
            CP_out = round(CP3, rt)
            V_out = round(V3, rt)
            a1_out = round(a1*10, rt)
            a2_out = round(a2/100, rt)
            a3_out = round(a3, rt)
            a4_out = round(a4/10000, rt)
            c1_out = round(c1, rt)
            c2_out = round(c2/10000, rt)
            wcon_out = round(wcon/100000, rt)
        Z_out = Z3


        write_output(filename=out_name, CATION=CATION,
                     LIGAND=LIGAND, nth_complex=3,
                     G=G_out, H=H_out, S=S_out,
                     CP=CP_out, V=V_out, a1=a1_out,
                     a2=a2_out, a3=a3_out, a4=a4_out,
                     c1=c1_out, c2=c2_out,
                     wcon=wcon_out, Z=Z_out)

    if BETA4 != 0:

        G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params(Z4, G4, H4, S4, CP4, V4)
        
        if sigfigs:
            G_out = sf(G, min(GC_sf, GL_sf))
            H_out = sf(H, min(HC_sf, HL_sf))
            S_out = sf(S, min(SC_sf, SL_sf))
            CP_out = sf(CP1, min(CPC_sf, CPL_sf))
            V_out = sf(V1, min(VC_sf, VL_sf))
            a1_out = sf(a1*10, min(VC_sf, VL_sf))
            a2_out = sf(a2/100, min(VC_sf, VL_sf))
            a3_out = sf(a3, min(VC_sf, VL_sf))
            a4_out = sf(a4/10000, min(VC_sf, VL_sf))
            c1_out = sf(c1, min(CPC_sf, CPL_sf))
            c2_out = sf(c2/10000, min(CPC_sf, CPL_sf))
            wcon_out = sf(wcon/100000, min(SC_sf, SL_sf))
        else:
            G_out = round(G, rt)
            H_out = round(H, rt)
            S_out = round(S, rt)
            CP_out = round(CP4, rt)
            V_out = round(V4, rt)
            a1_out = round(a1*10, rt)
            a2_out = round(a2/100, rt)
            a3_out = round(a3, rt)
            a4_out = round(a4/10000, rt)
            c1_out = round(c1, rt)
            c2_out = round(c2/10000, rt)
            wcon_out = round(wcon/100000, rt)
        Z_out = Z4

        df = write_output(filename=out_name, CATION=CATION,
                     LIGAND=LIGAND, nth_complex=4,
                     G=G_out, H=H_out, S=S_out,
                     CP=CP_out, V=V_out, a1=a1_out,
                     a2=a2_out, a3=a3_out, a4=a4_out,
                     c1=c1_out, c2=c2_out,
                     wcon=wcon_out, Z=Z_out)
        
        return df