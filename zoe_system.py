### GET DATA AND PASS TO sep_mol
# ONLY WANT TO DO ONCE PER XYZ
def systemData(path, File, check):
    from qcp.general import xyzPull

    # CUTOFF BETWEEN FRAGS
    distca = 1.7

    # PULL COORDS FROM XYZ
    coords = xyzPull(path, File)

    # CHECK CORRECT FRAGMENTS IF FMO
    if check:
        # UNLESS CHANGED
        comp   = 'nn'
        fragList, atmList, totChrg, totMult = sep_mol(coords, distca)
        if totMult == '?':
            while comp != 'n' and comp != 'y':
                comp = input('Does your cluster have '+ str(len(fragList)) +' frags? (y/n) ')
            while comp == 'n':
                distca = float(input("Covalent bond cutoff distance: "))
                fragList, atmList, totChrg, totMult = sep_mol(coords, distca)
                #print(fragList, totChrg)
                if totMult == '?':
                    comp = input('Does your cluster have '+ str(len(fragList)) +' frags? (y/n) ')
                else:
                    comp = 'y'


    # DOES NOT MATTER WITH OTHER CALC TYPES
    else:
        fragList, atmList, totChrg, totMult = sep_mol(coords, distca)

    return [fragList, atmList, totChrg, totMult]

### SEPARATE MOLECULES INTO IONS/MOLECULES
# RETURN totChrg, mult, chrgList, frgIndx, xyzData
# frgList["sym"], nu, x, y, z, ifrag, chrg : LIST OF DICTS
def sep_mol(coords, distca):

    import re, math
    import numpy       as np
    import collections as col
    from   qcp.pprint      import detec
    from   qcp.chemData    import pTable

    atmList      = []

    for ID, line in enumerate(coords):
        #print(line)
        atmDict        = {}
        atmDict["id"]  = ID
        atmDict["sym"] = line[0]
        atmDict["x"]   = float(line[1])
        atmDict["y"]   = float(line[2])
        atmDict["z"]   = float(line[3])
        atmDict["grp"] = False
        atmDict["con"] = []
        for sym, data in pTable.items():
            if atmDict["sym"] == sym:
                atmDict["nu"]  = data[0]
                atmDict["vdw"] = data[2]
        atmList.append(atmDict)

    natoms = len(atmList)
    dist   = np.zeros((natoms,natoms))
    # FIND DISTS BETWEEN ALL ATOMS i & j
    group = 0
    for val, i in enumerate(atmList):
        con = False
        # SECOND ATOM < DISTCA
        for vals, j in enumerate(atmList):
            if i != j:
                a = (i["x"]-j["x"])**2 + (i["y"]-j["y"])**2\
                + (i["z"]-j["z"])**2
                dist[val,vals] = math.sqrt(a)
                # IF PART OF THE SAME MOLECULE
                if dist[val, vals] < distca:
                    con = True
                    # CHECK IF CONNECTED USING VDW DATA
                    if dist[val, vals] < i["vdw"] + j["vdw"]:
                        if not j["id"] in i["con"]:
                            i["con"].append(j["id"])
                            j["con"].append(i["id"])
                    # IF NEITHER I NOR J PART OF A GROUP
                    # ADD THEM TO A NEW GROUP
                    if i["grp"] is False and j["grp"] is False:
                        i["grp"], j["grp"] = group, group
                        group += 1
                    # IF BOTH HAVE BEEN ASSIGNED TO A DIFF GROUP
                    elif not i["grp"] is False and not j["grp"] is False:
                        if i["grp"] != j["grp"]:
                            grp_chng = j["grp"]
                            for k in atmList:
                                if k["grp"] is grp_chng:
                                    k["grp"] = i["grp"]
                    # IF j NOT ASSIGNED
                    elif not i["grp"] is False and j["grp"] is False:
                        j["grp"] = i["grp"]
                    # IF i NOT ASSIGNED
                    elif not j["grp"] is False and i["grp"] is False:
                        i["grp"] = j["grp"]
        if not con:
            i["grp"] = group
            group += 1

    # REORDER GROUPS STARTING FROM 0
    # CREATE fragList FROM fragDicts
    print("-"*40)
    print("SYSTEM:")
    #
    totChrg  = 0
    totMult  = 1
    #
    nfrags = 0
    fragList = []
    for grp in range(group):
        fragDict = {}
        fragDict['ids'] = []
        fragDict['syms'] = []
        foundnew = False
        for atm in atmList:
            if atm["grp"] == grp:
                foundnew = True
                atm["grp"] = nfrags
                fragDict['syms'].append(atm['sym'])
                fragDict['ids'].append(atm['id'])

        if foundnew:
            fragDict['grp'] = nfrags
            # FIND CHARGE/MULT
            if isCation(fragDict['syms'], fragDict['ids']):
                fragDict["chrg"] = 1
                fragDict["mult"] = 1
                if totChrg != '?':
                    totChrg += 1

            elif isAnion(fragDict['syms'], fragDict['ids']):
                fragDict["chrg"] = -1
                fragDict["mult"] = 1
                if totChrg != '?':
                    totChrg += -1

            elif isNeutral(fragDict['syms'], fragDict['ids']):
                fragDict["chrg"] = 0
                fragDict["mult"] = 1

            elif isRadical(fragDict['syms'], fragDict['ids']):
                fragDict["chrg"] = 0
                fragDict["mult"] = 2
                totMult = 2

            else:
                # FIND CHEMICAL FORMULA
                syms = []
                chemForm = ''
                for sym in fragDict['syms']:
                    if sym not in syms:
                        numTimes = fragDict['syms'].count(sym)
                        chemForm += sym + str(numTimes)
                        syms.append(sym)
                detec("unknown", chemForm, " ".join(str(x) for x in fragDict['ids']))
                fragDict["chrg"] = '?'
                fragDict["mult"] = '?'
                totChrg = '?'
                totMult = '?'


            # FOR NEXT FRAGMENT
            fragList.append(fragDict)
            nfrags += 1

    print('-'*40)
    ### fragList = {'ids': [21], 'syms': ['Cl'], 'grp': 1, 'chrg': -1, 'mult': 1}
    #for i in fragList:
    #    print(i)

    ### atmList = {'id': 710, 'sym': 'H', 'x': 11.5282, 'y': 7.0276, 'z': -18.5563, 'grp': 80, 'nu': 1.0}
    #for i in atmList:
    #    print(i)
    #print(fragList, atmList, totChrg, totMult)
    return fragList, atmList, totChrg, totMult

#
#
# ----------- FUNCTIONS ADAPTED FROM PPQC Samual Tan


def isCation(a, b, q = False):

    import collections as col
    from qcp.chemData import CationDB
    from qcp.pprint   import detec

    #a = mol.atomListAsElem_Sym()
    isCat = False
    for key, cation in CationDB.items():
        if col.Counter(a) == col.Counter(cation):
            if not q:
                atmList = ''
                for i in b:
                    # START FROM 1 INSTEAD OF ZERO
                    atmList += str(i+1) + ' '
                detec("Cation detected", key, atmList)
            #return True
            isCat = True
            break
    return isCat

def isAnion(a, b, q = False):

    import collections as col
    from qcp.chemData import AnionDB
    from qcp.pprint   import detec

    # q for quiet
    #a = mol.atomListAsElem_Sym()
    # note the next only returns the first value--no duplicates allowed, or detected
    # checking done via Counter() from collections, no sorting required, duplicates included
    #return next((key for key, anion in AnionDB.items()
    #             if col.Counter(a) == col.Counter(anion)), None)
    isAni = False
    for key, anion in AnionDB.items():
        if col.Counter(a) == col.Counter(anion):
            if not q:
                atmList = ''
                for i in b:
                    # START FROM 1 INSTEAD OF ZERO
                    atmList += str(i+1) + ' '
                detec("Anion detected", key, atmList)
            isAni = True
            break
            #return True
        #else:
        #    print("no match ", col.Counter(a), col.Counter(anion))
        #    return False
    return isAni


def isNeutral(a, b, q = False):
    import collections as col
    from qcp.chemData import NeutralDB
    from qcp.pprint   import detec

    isNeu = False
    for key, molec in NeutralDB.items():
        if col.Counter(a) == col.Counter(molec):
            if not q:
                atmList = ''
                for i in b:
                    # START FROM 1 INSTEAD OF ZERO
                    atmList += str(i+1) + ' '
                detec("Neutral detected", key, atmList)
            isNeu = True
            break
    return isNeu


def isRadical(a, b, q = False):
    import collections as col
    from qcp.chemData import RadicalDB
    from qcp.pprint   import detec

    isRad = False
    for key, molec in RadicalDB.items():
        if col.Counter(a) == col.Counter(molec):
            if not q:
                atmList = ''
                for i in b:
                    # START FROM 1 INSTEAD OF ZERO
                    atmList += str(i+1) + ' '
                detec("Radical detected", key)
            isRad = True
            break
    return isRad
