# load libraries
import sys, os
import re
from datetime import datetime
import pandas as pd
import pkg_resources
import decimal
import math
import chemparse
import WORMutils


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


def __get_formula_ox_dict(name, df):
    
    entry = df[df["name"]==name]
    formula_ox = list(entry["formula_ox"])[0]
    formula_ox_split = formula_ox.split(" ")

    # remove any "" (empty strings) from list
    formula_ox_split = list(filter(("").__ne__, formula_ox_split))
    
    formula_ox_dict = {}
    for item in formula_ox_split:
        coeff = re.search(r'^\D*(\d+(?:\.\d+)?)', item)
        if coeff != None:
            coeff = coeff.group()
        if coeff == item or coeff == None:
            coeff = 1
            elem = item
        else:
            elem = item.split(coeff)[1]
        formula_ox_dict[elem] = float(coeff)
        
    return formula_ox_dict


def __write_output(filename, cation, ligand, nth_complex, G, H, S, CP, V,
                   a1, a2, a3, a4, c1, c2, wcon, Z, azero,
                   cation_dissrxn_dict, ligand_dissrxn_dict, thermo_data,
                   cation_formula_ox_dict, ligand_formula_ox_dict, replace,
                   skip_duplicates, ligand_name_abbrv_pairs):
    """
    Write output to a CSV.
    """

    duplicate_list = []
    cat_nocharge = re.sub("\+.*$", "", cation)
    lig_nocharge = re.sub("\-.*$", "", ligand)

    cat_formula = thermo_data[thermo_data["name"]==cation]["formula"].values[0]
    lig_formula = thermo_data[thermo_data["name"]==ligand]["formula"].values[0]
    cat_formula_nocharge = re.sub("\+.*$", "", cat_formula)
    lig_formula_nocharge = re.sub("\-.*$", "", lig_formula)
    
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
    
    azero_val = 4
    if Z == 2:
        azero_val = 6
    elif Z == 3:
        azero_val = 9
    elif Z == 4:
        azero_val = 11
    
    if nth_complex == 1:
        ligand_subscript = ""
    else:
        ligand_subscript = str(nth_complex)
    
    this_date = datetime.today().strftime('%Y%m%d') 
    
    complex_name = cat_nocharge + "(" + lig_nocharge + ")" + ligand_subscript + this_charge
    
    cation_dissrxn_dict_coeff = {name:coeff*1 for name, coeff in cation_dissrxn_dict.items()} # future: "1" can be replaced with number of cations
    ligand_dissrxn_dict_coeff = {name:coeff*nth_complex for name, coeff in ligand_dissrxn_dict.items()}
    
    cation_ligand_dissrxn_dict = {k: cation_dissrxn_dict_coeff.get(k, 0) + ligand_dissrxn_dict_coeff.get(k, 0) for k in set(cation_dissrxn_dict_coeff) | set(ligand_dissrxn_dict_coeff)}
    cation_ligand_dissrxn_formatted = ["{:.4f}".format(coeff) + " " + name for name, coeff in cation_ligand_dissrxn_dict.items()]
    
    cation_ligand_dissrxn_formatted = " ".join(cation_ligand_dissrxn_formatted)
    
    complex_dissrxn = "-1.0000 " + complex_name + " " + cation_ligand_dissrxn_formatted
    
    cation_formula_ox_dict_coeff = {name:coeff*1 for name, coeff in cation_formula_ox_dict.items()} # future: "1" can be replaced with number of cations
    ligand_formula_ox_dict_coeff = {name:coeff*nth_complex for name, coeff in ligand_formula_ox_dict.items()}
    complex_formula_ox_dict = {k: cation_formula_ox_dict_coeff.get(k, 0) + ligand_formula_ox_dict_coeff.get(k, 0) for k in set(cation_formula_ox_dict_coeff) | set(ligand_formula_ox_dict_coeff)}
    
    complex_formula_ox = []
    for elem,coeff in complex_formula_ox_dict.items():
        if coeff == 1:
            complex_formula_ox.append(str(elem))
        elif coeff == int(coeff):
            complex_formula_ox.append(str(int(coeff))+str(elem))
        else:
            complex_formula_ox.append(str(coeff)+str(elem))
            
    complex_formula_ox = " ".join(complex_formula_ox)

    data={
        "name": [complex_name],
        "abbrv": [cat_nocharge + "(" + lig_nocharge + ")" + ligand_subscript + this_charge],
        "formula": [cat_formula_nocharge + "(" + lig_formula_nocharge + ")" + ligand_subscript + this_charge],
        "state": ["aq"],
        "ref1": ["autocomplicator"],
        "ref2": ["NA"],
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
        "azero":[azero_val],
        "neutral_ion_type":[0],
        "dissrxn":[complex_dissrxn],
        "tag":[""],
        "formula_ox":[complex_formula_ox],
        "category_1":["inorganic_aq"],
        "category_2":[""],
    }
    df = pd.DataFrame(data)

    if isinstance(filename, str):
        
        file_exists = os.path.isfile(filename+'.csv')
        
        if file_exists:
            with open(filename+'.csv', 'a') as f:
    
                f_df = pd.read_csv(filename+'.csv')
    
                # check that a monoligand complex without parentheses isn't already
                # in the thermodynamic database. Warn or skip.
                keep = True
                complex_name_no_parentheses = complex_name.replace("(", "").replace(")", "")
    
                if complex_name_no_parentheses in list(f_df["name"]):
                    if skip_duplicates:
                        keep = False
                    else:
                        duplicate_list.append([complex_name, complex_name_no_parentheses])
    
                complex_name_organic_abbv = ""
    
                for name,abbv in ligand_name_abbrv_pairs:
                    # replace "Mg(acetate)+" with "Mg(Ac)+", etc. for the sake of
                    # checking for duplicate entries
                    complex_name_organic_abbv = complex_name.replace(name, abbv)
    
                if complex_name_organic_abbv in list(f_df["name"]):
                    if skip_duplicates:
                        keep = False
                    else:
                        duplicate_list.append([complex_name, complex_name_organic_abbv])
                
                if keep and replace and complex_name in list(f_df["name"]):
                    cols = list(f_df.columns) 
                    f_df.iloc[f_df["name"]==complex_name, :] = df
    
                    f_df.to_csv(filename+'.csv', header=True, index=False)
    
                elif keep and complex_name not in list(f_df["name"]):
                    df.to_csv(f, header=False, index=False)
                    
                else:
                    pass
    
        else:
            df.to_csv(filename+'.csv', index=False)
            
    return df, duplicate_list


# # function to count significant digits
# # adapted from https://stackoverflow.com/questions/8142676/python-counting-significant-digits
# def find_sigfigs(x):
    
#     '''
#     Get the number of significant digits in a string representing a number up to
#     eight digits long.

#     Parameters
#     ----------
#     x : str
#         A string denoting a number. This can include scientific notation.
    
#     Parameters
#     ----------
#     int
#         The number of significant digits.
    
#     Examples
#     --------
#     >>> find_sigfigs("5.220")
#     4
    
#     This also takes into account scientific notation.
    
#     >>> find_sigfigs("1.23e+3")
#     3
    
#     Insignificant zeros are ignored.
    
#     >>> find_sigfigs("4000")
#     1
    
#     A decimal point denotes that zeros are significant.
    
#     >>> find_sigfigs("4000.")
#     4
#     '''
    
#     x = str(x)
    
#     # change all the 'E' to 'e'
#     x = x.lower()
#     if ('-' == x[0]):
#         x = x[1:]
#     if ('e' in x):
#         # return the length of the numbers before the 'e'
#         myStr = x.split('e')
#         return len(myStr[0]) - 1  # to compenstate for the decimal point
#     else:
#         # put it in e format and return the result of that
#         ### NOTE: because of the 8 below, it may do crazy things when it parses 9 sigfigs
#         n = ('%.*e' % (8, float(x))).split('e')
#         # remove and count the number of removed user added zeroes. (these are sig figs)
#         if '.' in x:
#             s = x.replace('.', '')
#             #number of zeroes to add back in
#             l = len(s) - len(s.rstrip('0'))
#             #strip off the python added zeroes and add back in the ones the user added
#             n[0] = n[0].rstrip('0') + ''.join(['0' for num in range(l)])
#         else:
#             #the user had no trailing zeroes so just strip them all
#             n[0] = n[0].rstrip('0')
#         #pass it back to the beginning to be parsed
#     return find_sigfigs('e'.join(n))



def find_sigfigs(number_string):
    """
    Count significant figures (cannot handle scientific notation)
    """
    
    if isinstance(number_string, float):
        return ""

    if "-" in number_string:
        # delete negative symbol
        number_string = number_string.replace("-", "")
    
    if "." in number_string:
        # delete decimal point
        number_string = number_string.replace(".", "")
        # delete leading zeros
        number_string = number_string.lstrip("0")
    else:
        # delete leading and trailing zeros
        number_string = number_string.strip("0")

    return len(number_string)



def calc_params_WORMutils(Z, G, H, S, CP, V,
                          a1, a2, a3, a4, c1, c2, omega,
                          water_model="SUPCRT92"):

    kwargs = {"V": V, "Cp": CP, "Gf": G, "Hf": H, "Saq": S,
              "a1":a1, "a2":a2, "a3":a3, "a4":a4, "c1":c1, "c2":c2,
              "aq_complex": True, "HKF_scale":False, "Z":Z}
    
    if water_model == "DEW":
        kwargs["DEW"] = True

    if isinstance(omega, float):
        kwargs["omega"] = omega
    
    kwargs["phase_TrPr"] = "cr"

    out_dict, eqs = WORMutils.find_HKF(**kwargs)
    
    G = out_dict["G"]
    H = out_dict["H"]
    S = out_dict["S"]
    a1 = out_dict["a1"]
    a2 = out_dict["a2"]
    a3 = out_dict["a3"]
    a4 = out_dict["a4"]
    c1 = out_dict["c1"]
    c2 = out_dict["c2"]
    wcon = out_dict["omega"]
    
    return G, H, S, a1, a2, a3, a4, c1, c2, wcon


def calc_params(Z, G, H, S, CP, V, a1, a2, a3, a4, c1, c2, omega, water_model="SUPCRT92"):
    # The following are the Born functions at 25 C, and various other constants.

    AK = float('nan')
    RX = float('nan')
    sigma = float('nan')
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
    else:
        gamma=float('NaN')
    
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
            return
        
        if not math.isnan(S):
            # GB: RX is always 'NaN' because it's set at the start of the function. Why is this here?
            if not math.isnan(RX):
                re = RX+R*gamma
            else:
                ra=(Z**2 *(eta * Y - 100.)/(S - alphaz))
                ire=int(100.*ra +.5)
                re=ire/100.
        else:
            # GB: RX is always 'NaN' because it's set at the start of the function. Why is this here?
            if not math.isnan(RX):
                re = RX+Z*gamma
                S =(Z**2)*(eta*Y - 100.) / re + alphaz
            else:
                print('IF Z IS 0.0, EITHER S OR RX MUST BE GIVEN!')
                return

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

    # AK is always NaN because it is set at the beginning of the function.
    # Is there a reason it is this way in the original code?
    # This section will always be skipped!
    if not math.isnan(AK):
        aks=wcon*xN*conv
        akn=AK-aks

        # sigma is always 0 because it is set at the beginning of the function.
        # Is there a reason it is this way in the original code?
        if not math.isnan(sigma):
            a2=17.19E+4*akn+421.1
            # rounding
            ia2 = int(a2*100 + .5)
            a2 = ia2/100

            if water_model == "DEW":
                if Z == 0:
                    a1 = (0.1942*vn + 1.52)/10
                else:
                    a1 = ((0.1304*abs(Z) - 0.0217)*vn  + (1.4567*abs(Z)+0.6187))/10
            else:
                a1 = sigma/conv-(a2/pfunk)
                
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

            if water_model == "DEW":
                if Z == 0:
                    a1 = (0.1942*vn + 1.52)/10
                else:
                    a1 = ((0.1304*abs(Z) - 0.0217)*vn  + (1.4567*abs(Z)+0.6187))/10
            else:
                a1=sigma/conv-(a2/pfunk)
                
            # rounding
            ia1 = int(a1*100000 + .5)
            a1 = ia1/100000
    else:

        # sigma is always NaN because it is set at the beginning of the function.
        # Is there a reason it is this way in the original code?
        # This section will always be skipped!
        if not math.isnan(sigma):
    
            if water_model == "DEW":
                if Z == 0:
                    a1 = (0.1942*vn + 1.52)/10
                else:
                    a1 = ((0.1304*abs(Z) - 0.0217)*vn  + (1.4567*abs(Z)+0.6187))/10
            else:
                a1 = 1.3684E-2*vn+0.1765
            
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
            
            if water_model == "DEW":
                if Z == 0:
                    a1 = (0.1942*vn + 1.52)/10
                else:
                    a1 = ((0.1304*abs(Z) - 0.0217)*vn  + (1.4567*abs(Z)+0.6187))/10
            else:
                a1 = 1.3684E-2*vn+0.1765
                
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


def complicate(cation=None, ligand=None, beta=None, sass=None, out_name=None,
               G=None, H=None, S=None, V=None, Cp=None, a1=None, a2=None, a3=None, a4=None,
               c1=None, c2=None, omega=None, azero=4, rt=3, sigfigs=False, df=None,
               data_path=None, ligand_abbrv_data_path=None, water_model="SUPCRT92",
               correct_basis=True, replace=False, skip_duplicates=True, print_warnings=True,
               legacy_HKF_method=False):
    
    """
    Estimate the thermodynamic properties and Helgeson-Kirkham-Flowers (HKF)
    equation of state parameters for aqueous inorganic complexes.
    
    Parameters
    ----------
    cation : str
        Name of the cation in the aqueous complex.
        
    ligand : str
        Name of the ligand in the aqueous complex.
    
    beta : list of float
        List of four association equilibrium constants for the first, second,
        third, and fourth association, respectively. The first value should
        correspond to the equilibrium constant of a reaction like this:
        Mg+2 + Cl- = MgCl+
        While the second should represent a reaction written like this:
        Mg+2 + 2Cl- = MgCl2
        And so on, up to the fourth association.
       
    sass : list of float, optional
        List of four entropies of association (in cal mol-1 K-1) for the first,
        second, third, and fourth aqueous complex, respectively. If one or more
        values are not provided, they will be estimated from the properties of
        the cation, ligand, and association constants. The first entropy of
        association would be for a reaction like this:
        Mg+2 + Cl- = MgCl+
        While the second would be for a reaction written like this:
        Mg+2 + 2Cl- = MgCl2
        And so on, up to the fourth association.

    G : list of numeric, optional
        List of four partial molal Gibbs free energies of formation (in cal
        mol-1) corresponding to the aqueous complex formed by the first, second,
        third, and fourth association. If values are not provided, they will be
        estimated from the properties of the cation, ligand, and association
        constants.

    H : list of numeric, optional
        List of four partial molal enthalpies of formation (in cal
        mol-1) corresponding to the aqueous complex formed by the first, second,
        third, and fourth association. If values are not provided, they will be
        estimated from the properties of the cation, ligand, and association
        constants.

    Cp : list of float, optional
        List of four partial molal isobaric heat capacities (in cal mol-1 K-1)
        corresponding to the aqueous complex formed by the first, second, third,
        and fourth association. If values are not provided, they will be
        estimated from the properties of the cation and ligand.

    V : list of float, optional
        List of four partial molal volumes (in cm3 mol-1) corresponding to the
        aqueous complex formed by the first, second, third, and fourth
        association. If values are not provided, they will be estimated from the
        properties of the cation and ligand.

    a1, a2, a3, a4, c1, c2, omega : list of float, optional
        Parameters for the revised Helgeson Kirkam Flowers
        equation of state. The user may provide a list of four values
        corresponding to the aqueous complex formed by the first, second, third,
        and fourth association. If values are not provided, they will be
        estimated using published correlation methods.
    
    azero : float, default 4
        Aqueous species hard core diameter, a parameter in the Debeye-Hückel
        B-dot equation (Helgeson 1969) used to calculate activity coefficients.
        Units are in angstroms.
    
    out_name : str
        Name of the CSV output file to write.
    
    rt : int, default 3
        Round output values to how many decimal places? Ignored if `sigfigs` is
        True.
    
    sigfigs : bool, default False
        Experimental. Round output values to a number of significant figures
        determined by number of significant figures in values used to make the
        estimation?

    df : pandas.DataFrame, optional
        A dataframe containing columns for "Metal", "Ligand", "BETA_1", "BETA_2",
        "BETA_3", "BETA_4", "S_1", "S_2", "S_3", "S_4". If a dataframe is
        provided to this parameter, thermodynamic estimates will be created for
        the metal-ligand complex given on each row.
    
    data_path : str, optional
        File path and name of a custom thermodynamic database CSV. If undefined,
        the default WORM database will be used.

    water_model : str, default "SUPCRT92"
        Can be either "SUPCRT92" or "DEW" (Deep Earth Water) model. Influences
        the compatibility of estimated aqueous complex equation of state
        parameters for the chosen water model. If "DEW" is selected, the a1
        parameter is estimated according to Eq 16 in Sverjensky et al 2014
        (https://doi.org/10.1016/j.gca.2013.12.019)
    
    correct_basis : bool, default True
        If a cation or ligand is not a strict or auxiliary basis species in the
        thermodynamic database, ensure the dissociation reaction of the
        resulting complex includes only basis species? The default is True
        because complexes that are sent to aqueous equilibration package AqEquil
        require dissociation reactions into strict or auxiliary basis species.
    
    skip_duplicates : bool, default True
        Detect whether a complex with a monovalent ligand is already in the
        thermodynamic database without parentheses. E.g., if "Ca(HCO3)+" is being
        added, check that a species named "CaHCO3+" isn't already in the database.
        If set to False, the species will be added anyway, but the user will
        be warned. If set to True, the species will not be added.
    
    Returns
    -------
    df_out : pandas dataframe
        A dataframe with estimated properties and HKF parameters of the aqueous
        inorganic complex.
    """

    if isinstance(df, pd.DataFrame):

        def complicator_check_assign(prop, row):
            if prop in row.index:
                if not math.isnan(row[prop]):
                    return row[prop]
                else:
                    return float('NaN')
            else:
                return float('NaN')

        
        warnings = []
        df_out_list = []
        for index, row in df.iterrows():
            if 'BETA_1' in row.index:
                BETA1 = row['BETA_1']
            else:
                BETA1 = float('NaN')
            if 'BETA_2' in row.index:
                BETA2 = row['BETA_2']
            else:
                BETA2 = float('NaN')
            if 'BETA_3' in row.index:
                BETA3 = row['BETA_3']
            else:
                BETA3 = float('NaN')
            if 'BETA_4' in row.index:
                BETA4 = row['BETA_4']
            else:
                BETA4 = float('NaN')
            if math.isnan(BETA1) and math.isnan(BETA2) and math.isnan(BETA3) and math.isnan(BETA4):
                continue

            if 'Sass_1' in row.index:
                sass1 = row['Sass_1']
            else:
                sass1 = float('NaN')
            if 'Sass_2' in row.index:
                sass2 = row['Sass_2']
            else:
                sass2 = float('NaN')
            if 'Sass_3' in row.index:
                sass3 = row['Sass_3']
            else:
                sass3 = float('NaN')
            if 'Sass_4' in row.index:
                sass4 = row['Sass_4']
            else:
                sass4 = float('NaN')

            G_1 = complicator_check_assign("G_1", row)
            H_1 = complicator_check_assign("H_1", row)
            S_1 = complicator_check_assign("S_1", row)
            V_1 = complicator_check_assign("V_1", row)
            Cp_1 = complicator_check_assign("Cp_1", row)
            a1_1 = complicator_check_assign("a1_1", row)
            a2_1 = complicator_check_assign("a2_1", row)
            a3_1 = complicator_check_assign("a3_1", row)
            a4_1 = complicator_check_assign("a4_1", row)
            c1_1 = complicator_check_assign("c1_1", row)
            c2_1 = complicator_check_assign("c2_1", row)
            omega_1 = complicator_check_assign("omega_1", row)
            G_2 = complicator_check_assign("G_2", row)
            H_2 = complicator_check_assign("H_2", row)
            S_2 = complicator_check_assign("S_2", row)
            V_2 = complicator_check_assign("V_2", row)
            Cp_2 = complicator_check_assign("Cp_2", row)
            a1_2 = complicator_check_assign("a1_2", row)
            a2_2 = complicator_check_assign("a2_2", row)
            a3_2 = complicator_check_assign("a3_2", row)
            a4_2 = complicator_check_assign("a4_2", row)
            c1_2 = complicator_check_assign("c1_2", row)
            c2_2 = complicator_check_assign("c2_2", row)
            omega_2 = complicator_check_assign("omega_2", row)
            G_3 = complicator_check_assign("G_3", row)
            H_3 = complicator_check_assign("H_3", row)
            S_3 = complicator_check_assign("S_3", row)
            V_3 = complicator_check_assign("V_3", row)
            Cp_3 = complicator_check_assign("Cp_3", row)
            a1_3 = complicator_check_assign("a1_3", row)
            a2_3 = complicator_check_assign("a2_3", row)
            a3_3 = complicator_check_assign("a3_3", row)
            a4_3 = complicator_check_assign("a4_3", row)
            c1_3 = complicator_check_assign("c1_3", row)
            c2_3 = complicator_check_assign("c2_3", row)
            omega_3 = complicator_check_assign("omega_3", row)
            G_4 = complicator_check_assign("G_4", row)
            H_4 = complicator_check_assign("H_4", row)
            S_4 = complicator_check_assign("S_4", row)
            V_4 = complicator_check_assign("V_4", row)
            Cp_4 = complicator_check_assign("Cp_4", row)
            a1_4 = complicator_check_assign("a1_4", row)
            a2_4 = complicator_check_assign("a2_4", row)
            a3_4 = complicator_check_assign("a3_4", row)
            a4_4 = complicator_check_assign("a4_4", row)
            c1_4 = complicator_check_assign("c1_4", row)
            c2_4 = complicator_check_assign("c2_4", row)
            omega_4 = complicator_check_assign("omega_4", row)

            a1_1 = a1_1/10
            a2_1 = a2_1/10**-2
            a3_1 = a3_1
            a4_1 = a4_1/10**-4
            c1_1 = c1_1
            c2_1 = c2_1/10**-4
            omega_1 = omega_1/10**-5
            a1_2 = a1_2/10
            a2_2 = a2_2/10**-2
            a3_2 = a3_2
            a4_2 = a4_2/10**-4
            c1_2 = c1_2
            c2_2 = c2_2/10**-4
            omega_2 = omega_2/10**-5
            a1_3 = a1_3/10
            a2_3 = a2_3/10**-2
            a3_3 = a3_3
            a4_3 = a4_3/10**-4
            c1_3 = c1_3
            c2_3 = c2_3/10**-4
            omega_3 = omega_3/10**-5
            a1_4 = a1_4/10
            a2_4 = a2_4/10**-2
            a3_4 = a3_4
            a4_4 = a4_4/10**-4
            c1_4 = c1_4
            c2_4 = c2_4/10**-4
            omega_4 = omega_4/10**-5

            cation = str(row['Metal'])
            ligand = str(row['Ligand'])
            beta = [BETA1, BETA2, BETA3, BETA4]
            G = [G_1, G_2, G_3, G_4]
            H = [H_1, H_2, H_3, H_4]
            S = [S_1, S_2, S_3, S_4]
            sass = [sass1, sass2, sass3, sass4]
            V = [V_1, V_2, V_3, V_4]
            Cp = [Cp_1, Cp_2, Cp_3, Cp_4]
            a1 = [a1_1, a1_2, a1_3, a1_4]
            a2 = [a2_1, a2_2, a2_3, a2_4]
            a3 = [a3_1, a3_2, a3_3, a3_4]
            a4 = [a4_1, a4_2, a4_3, a4_4]
            c1 = [c1_1, c1_2, c1_3, c1_4]
            c2 = [c2_1, c2_2, c2_3, c2_4]
            omega = [omega_1, omega_2, omega_3, omega_4]
            
            df_out_row, warn_row = complicate(cation, ligand, beta, sass,
                    out_name, df=None, sigfigs=sigfigs, data_path=data_path,
                    replace=replace, skip_duplicates=skip_duplicates,
                    water_model=water_model, legacy_HKF_method=legacy_HKF_method,
                    G=G, H=H, S=S, V=V, Cp=Cp, a1=a1, a2=a2, a3=a3, a4=a4, c1=c1, c2=c2,
                    omega=omega, ligand_abbrv_data_path=ligand_abbrv_data_path,
                    print_warnings=False)
            
            df_out_list.append(df_out_row)
            warnings += warn_row

        df_out = pd.concat(df_out_list, ignore_index=True)
        
        if print_warnings:
            warnings = list(dict.fromkeys(warnings))
            for warning in warnings:
                print(warning)
        
        return df_out, warnings

    warning_list = []
    BETA1, BETA2, BETA3, BETA4 = beta

    # convert association constants (input) into dissociation constants used by
    # the autopp93s4 code.
    BETA1 = -BETA1
    BETA2 = -BETA2
    BETA3 = -BETA3
    BETA4 = -BETA4
    
    sass1, sass2, sass3, sass4 = sass

    G_1, G_2, G_3, G_4 = G
    H_1, H_2, H_3, H_4 = H
    S_1, S_2, S_3, S_4 = S
    Cp_1, Cp_2, Cp_3, Cp_4 = Cp
    V_1, V_2, V_3, V_4 = V
    a1_1, a1_2, a1_3, a1_4 = a1
    a2_1, a2_2, a2_3, a2_4 = a2
    a3_1, a3_2, a3_3, a3_4 = a3
    a4_1, a4_2, a4_3, a4_4 = a4
    c1_1, c1_2, c1_3, c1_4 = c1
    c2_1, c2_2, c2_3, c2_4 = c2
    omega_1, omega_2, omega_3, omega_4 = omega 

    eqs_used = [] # store list of strings describing which equations were used in the estimation
    
    if data_path == None:
        thermo_data = pkg_resources.resource_stream(__name__, "wrm_data.csv")
        thermo_data = pd.read_csv(thermo_data, dtype=object)
    else:
        thermo_data = pd.read_csv(data_path)

    if ligand_abbrv_data_path == None:
        ligand_abbrv = pkg_resources.resource_stream(__name__, "wrm_ligand_abbrv.csv")
        ligand_abbrv = pd.read_csv(ligand_abbrv)
    elif isinstance(ligand_abbrv_data_path, str):
        ligand_abbrv = pd.read_csv(ligand_abbrv_data_path)
        
    # a list of ligand names and abbreviations,
    # e.g., [["formate", "For"], ["acetate", "Ac"]], for the sake of
    # checking for duplicate entries in the thermodynamic database
    ligand_name_abbrv_pairs = [[name,abbv] for name,abbv in zip(ligand_abbrv["ligand"], ligand_abbrv["abbrv"])]
    
    cation_entry = thermo_data[thermo_data["name"]==cation]
    ligand_entry = thermo_data[thermo_data["name"]==ligand]
    
    if cation_entry.empty:
        msg = "Could not find an entry for " + cation + " in reference sheet."
        warning_list.append(msg)
        
    if ligand_entry.empty:
        msg = "Could not find an entry for " + ligand + " in reference sheet."
        warning_list.append(msg)

    if cation_entry.shape[0] == 0 or ligand_entry.shape[0] == 0:
        return pd.DataFrame(), warning_list
    
    if cation_entry["tag"].values[0] != "basis" and cation_entry["tag"].values[0] != "aux" and correct_basis:
        cation_dissrxn = cation_entry["dissrxn"].values[0]
        cation_dissrxn_split = cation_dissrxn.split(" ")
        cation_dissrxn_names = cation_dissrxn_split[1::2]
        cation_dissrxn_coeffs = cation_dissrxn_split[::2]
        
        cation_dissrxn_dict_prediv = {name:float(coeff) for name,coeff in zip(cation_dissrxn_names, cation_dissrxn_coeffs)}
        div_val = abs(cation_dissrxn_dict_prediv[cation])
        cation_dissrxn_dict = {name:coeff/div_val for name,coeff in cation_dissrxn_dict_prediv.items()}
        del cation_dissrxn_dict[cation]
        
        cation_dissrxn_dict = cation_dissrxn_dict/abs(cation_dissrxn_dict[cation])

        msg = (cation + " is not a strict or auxiliary basis species. "
              "Dissociation reactions of the complex will assume the "
              "cation dissociates according to: " + cation_dissrxn)
        warning_list.append(msg)
    else:
        cation_dissrxn_dict = {cation:1}

    if ligand_entry["tag"].values[0] != "basis" and ligand_entry["tag"].values[0] != "aux" and correct_basis:
        ligand_dissrxn = ligand_entry["dissrxn"].values[0]
        ligand_dissrxn_split = ligand_dissrxn.split(" ")
        ligand_dissrxn_names = ligand_dissrxn_split[1::2]
        ligand_dissrxn_coeffs = ligand_dissrxn_split[::2]
        
        ligand_dissrxn_dict_prediv = {name:float(coeff) for name,coeff in zip(ligand_dissrxn_names, ligand_dissrxn_coeffs)}
        div_val = abs(ligand_dissrxn_dict_prediv[ligand])
        ligand_dissrxn_dict = {name:coeff/div_val for name,coeff in ligand_dissrxn_dict_prediv.items()}
        del ligand_dissrxn_dict[ligand]

        msg = (ligand + " is not a strict or auxiliary basis species. "
              "Dissociation reactions of the complex will assume the "
              "ligand dissociates according to: " + ligand_dissrxn)
        warning_list.append(msg)
    else:
        ligand_dissrxn_dict = {ligand:1}
    
    # handle complex formula_ox
    cation_formula_ox_dict = __get_formula_ox_dict(cation, df=thermo_data)
    ligand_formula_ox_dict = __get_formula_ox_dict(ligand, df=thermo_data)
    
        
    ZC = str(cation_entry["z.T"].values[0])
    GC = str(cation_entry["G"].values[0])
    HC = str(cation_entry["H"].values[0])
    SC = str(cation_entry["S"].values[0])
    CPC = str(cation_entry["Cp"].values[0])
    VC = str(cation_entry["V"].values[0])
        
    if sigfigs:
        GC_sf = find_sigfigs(GC)
        HC_sf = find_sigfigs(HC)
        SC_sf = find_sigfigs(SC)
        CPC_sf = find_sigfigs(CPC)
        VC_sf = find_sigfigs(VC)
    else:
        GC_sf = GC
        HC_sf = HC
        SC_sf = SC
        CPC_sf = CPC
        VC_sf = VC
    
    ZL = str(ligand_entry["z.T"].values[0])
    GL = str(ligand_entry["G"].values[0])
    HL = str(ligand_entry["H"].values[0])
    SL = str(ligand_entry["S"].values[0])
    CPL = str(ligand_entry["Cp"].values[0])
    VL = str(ligand_entry["V"].values[0])
        
    GL_sf = find_sigfigs(GL)
    HL_sf = find_sigfigs(HL)
    SL_sf = find_sigfigs(SL)
    CPL_sf = find_sigfigs(CPL)
    VL_sf = find_sigfigs(VL)

    ZC, GC, HC, SC, CPC, VC = int(float(ZC)), float(GC), float(HC), float(SC), float(CPC), float(VC)
    ZL, GL, HL, SL, CPL, VL = int(float(ZL)), float(GL), float(HL), float(SL), float(CPL), float(VL)

    # Calculations for the first complex
    Z = ZC + ZL
    DELGR1 = (2.30259)*(1.98719)*(298.15)*BETA1
    G1 = DELGR1 + GC + GL
    eqs_used.append("DELGR1 = {} cal/mol = (2.30259)*(1.98719)*(298.15)*BETA1, Gibbs free energy of association of the first complex\n".format("{0:.5g}".format(DELGR1)))
    eqs_used.append("G1 = {} cal/mol = DELGR1 + GC + GL, standard state partial molal Gibbs free energy of formation of the first complex\n".format("{0:.5g}".format(G1)))

    # Now the entropy predictor starts for the first complex

    # if sass1 is available, use it. If not, predict it.
    if not math.isnan(sass1):
        # S1 = sass1
        # DELSR1 = S1 - SC - SL
        # DELS1 = DELSR1
        DELSR1 = sass1
        DELS1 = DELSR1
        S1 = SC + SL + DELS1
        eqs_used.append("Entropy of association Sass_1 = {}, third law entropy of the cation SC = {}, and ligand SL = {}, cal/mol/K\n".format(sass1, SC, SL))
        eqs_used.append("S1 = {} cal/mol/K = S1 = SC + SL + Sass_1, standard state partial molal third law entropy of the first complex\n".format("{0:.5g}".format(S1)))
    
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

    # Here the Cp predictor starts for the first complex.
    # Modified 8 Aug 1992 after latest revision Modified 17 March 1992 to
    # incorporate changes from last July 1991
    dz = 0.856*CPL - 2.1 + 45.3*ZC
    DELCP1 = 1.25*CPC + dz
    CP1 = DELCP1 + CPC + CPL

    # Here the Volume predictor starts for the first complex
    DELVR1  = 0.11419*VC + 8.9432 # Eq. 62 in SSH97
    V1 = DELVR1 + VC + VL
    Z1 = Z
    
    # Calculations for the second complex
    Z = Z + ZL

    DELGR2 = 2.30259*1.98719*298.15*BETA2
    G2 = DELGR2 + GC + (2.0*GL)

    # if sass2 is available, use it. If not, predict it.
    if not math.isnan(sass2):
        # S2=sass2
        # DELSR2 = S2 - SC - (2.0*SL)
        # DELS2 = DELSR2 - DELS1
        DELSR2 = sass2
        DELS2 = DELSR2 - DELS1
        S2 = SC + 2.0*SL + DELS1 + DELS2
    
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

    DELVR2  = 0.11419*V1 + 8.9432 # Eq. 62 in SSH97
    V2 = DELVR2 + DELVR1 + VC + (2.0*VL)
    Z2 = Z
    
    # Calculations for the third complex

    Z = Z + ZL  
    DELGR3 = 2.30259*1.98719*298.15*BETA3
    G3 = DELGR3 + GC + (3.0*GL)

    # if sass3 is available, use it. If not, predict it.
    if not math.isnan(sass3):
        # S3=sass3
        # DELSR3 = S3 - SC - (3.0*SL)
        # DELS3 = DELSR3 - DELS1 - DELS2
        DELSR3 = sass3
        DELS3 = DELSR3 - DELS1 - DELS2
        S3 = SC + 3.0*SL + DELS1 + DELS2 + DELS3
    
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

    DELVR3  = 0.11419*V2 + 8.9432 # Eq. 62 in SSH97
    V3 = DELVR3 + DELVR2 + DELVR1 + VC + (3.0*VL)
    Z3 = Z
    
    # Calculations for the fourth complex
    Z = Z + ZL

    DELGR4 = 2.30259*1.98719*298.15*BETA4
    G4 = DELGR4 + GC + (4.0*GL)
    
    # if sass4 is available, use it. If not, predict it.
    if not math.isnan(sass4):
        # S4=sass4
        # DELSR4 = S4 - SC - (4.0*SL)
        # DELS4 = DELSR4 - DELS1 - DELS2 - DELS3
        DELSR4 = sass4
        DELS4 = DELSR4 - DELS1 - DELS2 - DELS3
        S4 = SC + 4.0*SL + DELS1 + DELS2 + DELS3 + DELS4
        
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

    DELVR4  = 0.11419*V3 + 8.9432 # Eq. 62 in SSH97
    V4 = DELVR4 + DELVR3 + DELVR2 + DELVR1 + VC + (4.0*VL)
    Z4 = Z


    # override properties and parameters if they are supplied by the user
    if not math.isnan(G_1):
        G1 = G_1
    if not math.isnan(H_1):
        H1 = H_1
    if not math.isnan(S_1):
        S1 = S_1
    if not math.isnan(Cp_1):
        CP1 = Cp_1
    if not math.isnan(V_1):
        V1 = V_1
    
    # Now the partial molal properties are used in PARCOR to estimate equation of state coefficients for each complex.

    df_out_1 = pd.DataFrame()
    df_out_2 = pd.DataFrame()
    df_out_3 = pd.DataFrame()
    df_out_4 = pd.DataFrame()
    duplicate_list_1 = []
    duplicate_list_2 = []
    duplicate_list_3 = []
    duplicate_list_4 = []
    
    # The following section has been altered at GEOPIG so that there will be no output if Beta = 0.0
    if not math.isnan(BETA1):
        try:
            if legacy_HKF_method:
                G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params(Z1, G1, H1, S1, CP1, V1, a1_1, a2_1, a3_1, a4_1, c1_1, c2_1, omega_1, water_model=water_model)
            else:
                G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params_WORMutils(Z1, G1, H1, S1, CP1, V1, a1_1, a2_1, a3_1, a4_1, c1_1, c2_1, omega_1, water_model=water_model)
        except:
            return pd.DataFrame(), warning_list

        if not math.isnan(G):
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

            df_out_1, duplicate_list_1 = __write_output(filename=out_name,
                         cation=cation, ligand=ligand, nth_complex=1, G=G_out, H=H_out, S=S_out,
                         CP=CP_out, V=V_out, a1=a1_out, a2=a2_out, a3=a3_out, a4=a4_out, c1=c1_out, c2=c2_out,
                         wcon=wcon_out, Z=Z_out, azero=azero,
                         cation_dissrxn_dict=cation_dissrxn_dict,
                         ligand_dissrxn_dict=ligand_dissrxn_dict,
                         thermo_data=thermo_data,
                         cation_formula_ox_dict=cation_formula_ox_dict,
                         ligand_formula_ox_dict=ligand_formula_ox_dict,
                         replace=replace,
                         skip_duplicates=skip_duplicates,
                         ligand_name_abbrv_pairs=ligand_name_abbrv_pairs)

    # override properties and parameters if they are supplied by the user
    if not math.isnan(G_2):
        G2 = G_2
    if not math.isnan(H_2):
        H2 = H_2
    if not math.isnan(S_2):
        S2 = S_2
    if not math.isnan(Cp_2):
        CP2 = Cp_2
    if not math.isnan(V_2):
        V2 = V_2
    
    if not math.isnan(BETA2):
        try:
            if legacy_HKF_method:
                G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params(Z2, G2, H2, S2, CP2, V2, a1_2, a2_2, a3_2, a4_2, c1_2, c2_2, omega_2, water_model=water_model)
            else:
                G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params_WORMutils(Z2, G2, H2, S2, CP2, V2, a1_2, a2_2, a3_2, a4_2, c1_2, c2_2, omega_2, water_model=water_model)
        except:
            return pd.DataFrame(), warning_list

        if not math.isnan(G):
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
        
            df_out_2, duplicate_list_2 = __write_output(filename=out_name,
                         cation=cation,
                         ligand=ligand, nth_complex=2,
                         G=G_out, H=H_out, S=S_out,
                         CP=CP_out, V=V_out, a1=a1_out,
                         a2=a2_out, a3=a3_out, a4=a4_out,
                         c1=c1_out, c2=c2_out,
                         wcon=wcon_out, Z=Z_out, azero=azero,
                         cation_dissrxn_dict=cation_dissrxn_dict,
                         ligand_dissrxn_dict=ligand_dissrxn_dict,
                         thermo_data=thermo_data,
                         cation_formula_ox_dict=cation_formula_ox_dict,
                         ligand_formula_ox_dict=ligand_formula_ox_dict,
                         replace=replace,
                         skip_duplicates=skip_duplicates,
                         ligand_name_abbrv_pairs=ligand_name_abbrv_pairs)

    # override properties and parameters if they are supplied by the user
    if not math.isnan(G_3):
        G3 = G_3
    if not math.isnan(H_3):
        H3 = H_3
    if not math.isnan(S_3):
        S3 = S_3
    if not math.isnan(Cp_3):
        CP3 = Cp_3
    if not math.isnan(V_3):
        V3 = V_3
    
    if not math.isnan(BETA3):
        try:
            if legacy_HKF_method:
                G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params(Z3, G3, H3, S3, CP3, V3, a1_3, a2_3, a3_3, a4_3, c1_3, c2_3, omega_3, water_model=water_model)
            else:
                G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params_WORMutils(Z3, G3, H3, S3, CP3, V3, a1_3, a2_3, a3_3, a4_3, c1_3, c2_3, omega_3, water_model=water_model)
        except:
            return pd.DataFrame(), warning_list

        if not math.isnan(G):
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

            df_out_3, duplicate_list_3 = __write_output(filename=out_name,
                         cation=cation,
                         ligand=ligand, nth_complex=3,
                         G=G_out, H=H_out, S=S_out,
                         CP=CP_out, V=V_out, a1=a1_out,
                         a2=a2_out, a3=a3_out, a4=a4_out,
                         c1=c1_out, c2=c2_out,
                         wcon=wcon_out, Z=Z_out, azero=azero,
                         cation_dissrxn_dict=cation_dissrxn_dict,
                         ligand_dissrxn_dict=ligand_dissrxn_dict,
                         thermo_data=thermo_data,
                         cation_formula_ox_dict=cation_formula_ox_dict,
                         ligand_formula_ox_dict=ligand_formula_ox_dict,
                         replace=replace,
                         skip_duplicates=skip_duplicates,
                         ligand_name_abbrv_pairs=ligand_name_abbrv_pairs)

    # override properties and parameters if they are supplied by the user
    if not math.isnan(G_4):
        G4 = G_4
    if not math.isnan(H_4):
        H4 = H_4
    if not math.isnan(S_4):
        S4 = S_4
    if not math.isnan(Cp_4):
        CP4 = Cp_4
    if not math.isnan(V_4):
        V4 = V_4

    if not math.isnan(BETA4):
        try:
            if legacy_HKF_method:
                G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params(Z4, G4, H4, S4, CP4, V4, a1_4, a2_4, a3_4, a4_4, c1_4, c2_4, omega_4, water_model=water_model)
            else:
                G, H, S, a1, a2, a3, a4, c1, c2, wcon = calc_params_WORMutils(Z4, G4, H4, S4, CP4, V4, a1_4, a2_4, a3_4, a4_4, c1_4, c2_4, omega, water_model=water_model)
        except:
            return pd.DataFrame(), warning_list

        if not math.isnan(G):
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
            
            df_out_4, duplicate_list_4 = __write_output(filename=out_name, cation=cation,
                         ligand=ligand, nth_complex=4,
                         G=G_out, H=H_out, S=S_out,
                         CP=CP_out, V=V_out, a1=a1_out,
                         a2=a2_out, a3=a3_out, a4=a4_out,
                         c1=c1_out, c2=c2_out,
                         wcon=wcon_out, Z=Z_out, azero=azero,
                         cation_dissrxn_dict=cation_dissrxn_dict,
                         ligand_dissrxn_dict=ligand_dissrxn_dict,
                         thermo_data=thermo_data,
                         cation_formula_ox_dict=cation_formula_ox_dict,
                         ligand_formula_ox_dict=ligand_formula_ox_dict,
                         replace=replace,
                         skip_duplicates=skip_duplicates,
                         ligand_name_abbrv_pairs=ligand_name_abbrv_pairs)

    df_out = pd.concat([df_out_1, df_out_2, df_out_3, df_out_4], ignore_index=True)
    
    duplicate_list = duplicate_list_1 + duplicate_list_2 + duplicate_list_3 + duplicate_list_4
    if len(duplicate_list) > 0:
        
        msg = ("Certain complexes that have just been estimated with the "
               "Complicator are likely duplicates of complexes that are already "
               "in the thermodynamic database. Set skip_duplicates=True "
               "and rerun to skip and exclude these complexes.")
        warning_list.append(msg)

        for dup_pair in duplicate_list:
            complex_name = dup_pair[0]
            complex_name_no_parentheses = dup_pair[1]
            
            warning_list.append(complex_name + " appears to be a duplicate of " + complex_name_no_parentheses + " already present in the thermodynamic database.")

    
    return df_out, warning_list