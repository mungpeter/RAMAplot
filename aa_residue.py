#!/usr/bin/env python3

## 20.02.26 - duplicate those that have 1- or 2-letter PDB_ID code

##########################################################################
# Natural amino acid 3-letter <-> 1-letter conversion
aa_dict = {
  'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q',
  'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'MET':'M',
  'PRO':'P', 'PHE':'F', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y',
  'LEU':'L', 'VAL':'V', 'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP',
  'C':'CYS', 'Q':'GLN', 'E':'GLU', 'G':'GLY', 'H':'HIS', 'I':'ILE',
  'K':'LYS', 'M':'MET', 'P':'PRO', 'F':'PHE', 'S':'SER', 'T':'THR',
  'W':'TRP', 'Y':'TYR', 'L':'LEU', 'V':'VAL',
}

# Natural amino acid 3-letter <-> 1-letter conversion
def AA(resname):
  if resname in aa_dict:
    return aa_dict[resname]
  else:
    return ''

##########################################################################
# Unnatural amino acid <-> corresponding natural amino acid conversion
ua_dict = {
  'CSO':'CYS',  # S-hydroxycysteine
  'OCS':'CYS',  # cysteine sulfonic acid
  'CSX':'CYS',  # S-oxy cysteine
  'OCY':'CYS',  # hydroxyethyl cysteine
  'CSE':'CYS',  # selenocysteine (obsolete)
  'CSS':'CYS',  # S-mercapto cysteine
  'CSW':'CYS',  # cystein-s-dioxide (obsolete)
  'CME':'CYS',  # S,S-(2-hydroxyethyl)thiocysteine
  'CMT':'CYS',  # O-methylcysteine
  'CY0':'CYS',  # S-{3-[(4-ANILINOQUINAZOLIN-6-YL)AMINO]-3-OXOPROPYL}-L-CYS
  'CAF':'CYS',  # S-dimethylarsinoyl cysteine
  'CAS':'CYS',  # S-(dimethylarsenic) cysteine
  '2CO':'CYS',  # S-hydroperoxycysteine

  'AME':'MET',  # N-acetylmethionine
  'MHO':'MET',  # S-oxymethionine
  'MSE':'MET',  # selenomethionine

  'KCX':'LYS',  # lysine NZ-carboxylic acid
  'LGY':'LYS',  # (E)-N-6-(4-oxobytylidene)-L-lysine
  'TRG':'LYS',  # L-(N,N-dimethyl)lysine (obsolete)
  'ALY':'LYS',  # N(6)-acetyl lysine
  'CIR':'LYS',  # citrulline
  'MLY':'LYS',  # N6,N6-dimethyl lysine

  'NMM':'ARG',  # N6-methyl argenine

  'FCL':'PHE',  # 3-chloro-L-phenylalanine
  'PFF':'PHE',  # 4-fluoro-L-phenylalanine
  'MEA':'PHE',  # N-methylphenylalanine

  'B3L':'LEU',  # (3S)-3-amino-5-methylhexanoic acid
  'L3O':'LEU',  # (2S,3S)-3-amino-2-hydroxy-5-methylhexanoic acid
  'MK8':'LEU',  # 2-methyl-L-norleucine

  'CSD':'ALA',  # 3-sulfino alnine
  'SCS':'ALA',  # 3-(ethyldisulfanyl)-L-alanine

  'TPO':'THR',  # phospho-threoine

  'SEP':'SER',  # phospho-serine
  '2RX':'SER',  # O-thiophosphono-L-serine

  'PTR':'TYR',  # phospho-tyrosine
  'TYI':'TYR',  # 3,5-diiodotyrosine

  'NEP':'HIS',  # N1-phosphohistidine
}

# Convert un-natural amino acid 3-letter -> natural amino acid 3-letter
def UnnaturalAA(name):
  if name in ua_dict:
    return ua_dict[name]
  else:
    return False

# Check if amino acid is un-natural amino acid
def CheckUnnaturalAA(name):
  if name in ua_dict:
    return True
  else:
    return False


##########################################################################
# Known salt ions and detergents/lipids in cyrstallography
sh_dict = {
'OSV':'',  # ruthenium octasporine 4
'DWC':'',  # pyridocarbazole cyclopentadienyl osimum
'RPS':'',  # phtalimide ruthenium complex
'RUI':'',  # cyclopentadienyl(carbox monoxide) .... ruthenium
'DVT':'',   # DECAVANADATE
'EMC':'',   # ethyl mercury ion

'0OA':'',   # phosphate lipid + sugar + cyclohexane
'0O9':'',   # phosphate lipid + sugar
'LDA':'',   # lauryl dimethylamine-N-oxide
'HEX':'',   # hexane
'UND':'',   # undecane
'FAR':'',   # frnesyl
'MYR':'',   # myristic acid
'0O8':'',   # (1,1-dimethylpiperidin-1-ium-4-yl)octadecyl hydrogen phosphate

'DTD':'',   # dithiane diol -   (4R,5R)-1,2-dithiane-4,5-diol
'BME':'',   # beta-mercaptoethanol
'DTT':'',   # 2,3-dihydroxy-1,4-dithiobutane
'DTV':'',   # 2S,3S-1,4-dimercaptobutane-2,3-diol
'DHL':'',   # 2-aminoethanethol
'MPT':'',   # beta-mercaptopropionic acid
'SGM':'',   # (2R)-3-sulfanylpropane-1,2-diol
'COM':'',   # 2-sulfanylethanesulfonic acid
'DTU':'',   # (2S,3R)-1,4-bis-sulfanylbutane-2,3-diol

'DMS':'',   # DMSO
'ODO':'',   # 4-[(E)-C-methyl-N-oxidanyl-carbonimidoyl]benzene-1,3-diol
'MES':'',   # 2-morpholin-4-ium-4-ylethanesulfonate
'IMD':'',   # imidazole
'EPE':'',   # 4-(2-hydroxyethyl)-1-piperazine ethanesulfonic acid
'NHE':'',   # 2-(N-cyclohexylamino)ethane sulfonic acid

'FMT':'',   # fomic acid
'FUM':'',   # (E)-but-2-enedioic acid
'ACE':'',   # acetyl group
'ACT':'',   # acetate ion
'ACY':'',   # acetic acid
'TFA':'',   # trifluoroacetic acid
'AZ1':'',   # azelaic acid
'GOA':'',   # glycolic acid
'AKG':'',   # 2-oxoglutaric acid
'CIT':'',   # citrate acid
'FLC':'',   # citrate ion
'OEG':'',   # 2,2'-oxydiacetic acid
'MLA':'',   # malonic acid
'MLI':'',   # malonate ion
'TLA':'',   # L-tartaric acid
'TAR':'',   # D-tartaric acid
'SRT':'',   # S,R meso-tartaric acid
'SIN':'',   # succinic acid
'HCA':'',   # 3-hydroxy-3-carboxy-hexan-1,6-oic acid
'ABA':'',   # (2S)-2-azanylbutanoic acid - alpha-aminobutyric acid
'ACA':'',   # 6-aminohexanoic acid
'REL':'',   # D-glucuronic acid
'CAC':'',   # gamma-carboxy-glutamic acid
'66N':'',   # L-alaninamide
'OCE':'',   # octanedioic acid

'GLA':'',   # alpha-D-galactose
'GAL':'',   # beta-D-galactose
'GLC':'',   # alpha-D-glucose
'BGC':'',   # beta-D-glucose
'BG6':'',   # beta-D-glucose-6-phosphate
'BOG':'',   # B-octylglucoside
'BDP':'',   # beta-D-glycopyranuronic acid
'BMA':'',   # beta-D-mannose
'MAN':'',   # alpha-D-mannose
'FUL':'',   # beta-L-fuctose
'FUC':'',   # alpha-L-fuctose
'XYS':'',   # xylopyranose
'MAL':'',   # maltose
'SUC':'',   # sucrose
'NOJ':'',   # 1-deoxynojirimycin
'LMT':'',   # dodecyl beta-D-maltoside
'HTG':'',   # heptyl 1-thiohexopyranoside
'NAG':'',   # N-acetyl-D-glucosamine
'NDG':'',   # 2-(acetylamino)-2-deoxy-Alpha-D-glycopyranose
'NGZ':'',   # 2-(acetylamino)-2-deoxyl-alpha-L-glucopyranose
'MG8':'',   # N-octanoyl-N-methylglucamine
'BOG':'',   # B-octylglucoside
'BNG':'',   # B-nonylglucoside
'BCD':'',   # poly cyclic sugar, BETA-CYCLODEXTRIN
'HSJ':'',   # octyl beta-L-talopyranoside
'I6P':'',   # INOSITOL 1,2,3,4,5,6-HEXAKISPHOSPHATE

'2HT':'',   # 3-methylbenzonitrile
'2CH':'',   # 2-chlorophenol
'IPH':'',   # phenol
'YTP':'',   # 1-(4-hydroxy-3-methylphenyl)ethanone
'CAQ':'',   # benzene-1,2-diol

'ETA':'',   # ethanolamine
'EDO':'',   # 1,2-ethanediol
'BU3':'',   # R,R-2,3-butanediol
'3ZQ':'',   # 1S,2S-cyclohexane-1,2-diol
'BUD':'',   # (2R,3R)-butane-2,3-diol
'GOL':'',   # glycerol
'PEG':'',   # DI(HYDROXYETHYL)ETHER
'ETX':'',   # 2-ethoxyethanol
'PGE':'',   # triethylene glycol
'PGO':'',   # S-1,2-propanediol
'P4G':'',   # 1-ethoxy-2-(2-ethoxyethoxy)ethane
'PG4':'',   # tetraethyene glycol
'1PE':'',   # pentaethylene glcol
'P6G':'',   # hexaethylene glycol
'2PE':'',   # nonaethylene glycol
'PE5':'',   # 3,6,9,12,15,18,21,24-octaoxahexacosane-1ol
'P33':'',   # 3,6,9,12,15,18-hexaoxaicosane-1,20-diol
'PE4':'',   # 3,6,9,12,15,18,21-heptaoxapentacosane-1,2-diol
'7PE':'',   # 3,6,9,12,15,18,21-heptaoxapentacosane-1-ol
'BTB':'',   # 2-(bis(2-hydroxyethyl)amino)-2-(hydroxymethyl)propane-1,3-diol
'TRS':'',   # 2-amino-2-hydroxymethyl-propane-1,3-diol
'MRD':'',   # 4R-2-methylpentane-2,4-diol
'MPD':'',   # 4S-2-methyl-2,4-pentanediol
'TAM':'',   # tris(hydroxyethyl)aminomethane
'MOH':'',   # methanol
'EOH':'',   # ethanol
'IPA':'',   # isopropyl alcohol
'SBT':'',   # 2-butanol
'OCT':'',   # octanol
'F09':'',   # nonan-1-ol
'DIO':'',   # 1,4-diethylene dioxide
'ETF':'',   # trifluoroethanol
'PSE':'',   # O-phospho ethanolamine
'PG0':'',   # 2-(2-methoxyethoxy)ethanol
'MXE':'',   # 2-methoxyethanol
'PEU':'',   # poly PEG
'P4C':'',   # poly PEG ethanal
'1KA':'',   # 2-(2-hydroxyethyloxy)ethanal

'GLZ':'',   # amino-acetaldehyde
'PTL':'',   # pentanal
'GBL':'',   # gamma-butyrolactone

'16D':'',   # hexane-1,6-diamine
'SPD':'',   # spermidine
'PUT':'',   # 1,4-diaminobutane
'SPE':'',   # thermine
'TMA':'',   # tetramethylammonium ion
'PHU':'',   # 1-phenylurea
'UMS':'',   # tetrabutylammonium ion
'NH2':'',   # amino group
'NH3':'',   # ammonia
'NH4':'',   # ammonium ion
'T1A':'',   # tetraethylarsonium ion

' LI':'',   # lithium ion
' NA':'',   # sodium ion
'  K':'',   # potassium ion
' RB':'',   # rubidium ion
' CS':'',   # cesium ion
' MG':'',   # magnesium ion
' CA':'',   # calcium ion
' SR':'',   # strontium ion
' BA':'',   # barium ion
'LI':'',   # lithium ion
'NA':'',   # sodium ion
'K':'',   # potassium ion
'RB':'',   # rubidium ion
'CS':'',   # cesium ion
'MG':'',   # magnesium ion
'CA':'',   # calcium ion
'SR':'',   # strontium ion
'BA':'',   # barium ion

' MN':'',   # managese ion
' ZN':'',   # zinc ion
' CO':'',   # cobalt2 ion
' MO':'',   # molybdenum 
' CD':'',   # cadium ion
' NI':'',   # nickel ion
' PB':'',   # lead2 ion
' HG':'',   # mercery ion
' AU':'',   # gold ion
' FE':'',   # iron3 ion
'FE2':'',   # iron2 ion
' CU':'',   # copper2 ion
'CUA':'',   # Cu-Cu ion
'  W':'',   # tungstein ion
' TL':'',   # thallium ion
'YT3':'',   # yttrium(+3) cation
'MN':'',   # managese ion
'ZN':'',   # zinc ion
'CO':'',   # cobalt2 ion
'MO':'',   # molybdenum 
'CD':'',   # cadium ion
'NI':'',   # nickel ion
'PB':'',   # lead2 ion
'HG':'',   # mercery ion
'AU':'',   # gold ion
'FE':'',   # iron3 ion
'CU':'',   # copper2 ion
'W':'',   # tungstein ion
'TL':'',   # thallium ion


' CL':'',
'CL':'',
' BR':'',
'BR':'',
'IOD':'',   # iodide ion
' OH':'',   # hydroxide ion
'OH':'',   # hydroxide ion
'SCN':'',   # thiocyanate ion
'CYN':'',   # cyanide ion
'AZI':'',   # azide ion
'CO3':'',   # carbonate ion
'NO3':'',   # nitrate ion
'SO3':'',   # sulfite ion
'SO4':'',   # sulfate ion
'PO4':'',
'DPO':'',   # diphosphate
'AF3':'',   # aluminum floride
'ALF':'',   # tetrafluoroalumium ion
'SF4':'',   
'OHX':'',   # osmium3 hexamine
'BEF':'',   # beryllium trifluoride ion
'BO3':'',   # boric acid
'ARS':'',   # arsenic
'MGF':'',   # trifluoromagnesium
'PO2':'',   # HYPOPHOSPHITE
'WO4':'',   # TUNGSTATE(VI)ION

'HOH':'',
'UNK':'',   # unknown
'UNX':'',   # unknown
' UN':'',
'UN':'',
}

# Check if residue is salt or un-natural amino acid
def SaltAdditive(name):
  if UnnaturalAA(name):
    return UnnaturalAA(name)
  if name in sh_dict:
    return True
  else:
    return False

##########################################################################
# Nucleotides

nc_dict = {
  'ATP':'', # ADENOSINE-5'-TRIPHOSPHATE	
  'ADP':'', 
  'AMP':'',
  '01G':'', #
  '08T':'', #
  '0DC':'', #
  '0DG':'', #
  '0G4':'', #
  '0G8':'', #
  '0KX':'', # 2'-deoxy-5'-O-[(R)-hydroxy{[(R)-hydroxy(phosphonooxy)phosphoryl]amino}phosphoryl]cytidine	
  '0O2':'', # guanosine 5'-(tetrahydrogen triphosphate) 3'-(trihydrogen diphosphate)
  'ANP':'', # PHOSPHOAMINOPHOSPHONIC ACID-ADENYLATE ESTER	
  '  A':'', # ADENOSINE-5'-MONOPHOSPHATE	
  'A':'', # ADENOSINE-5'-MONOPHOSPHATE	
  'A23':'', # ADENOSINE-5'-PHOSPHATE-2',3'-CYCLIC PHOSPHATE	
  'A2M':'', # 2'-O-methyladenosine 5'-(dihydrogen phosphate)	
  'A2P':'', # ADENOSINE-2'-5'-DIPHOSPHATE	
  'A3A':'', # 2'DEOXY-ALPHA-ANOMERIC-ADENOSINE-5'-PHOSPHATE	
  'A3P':'', # ADENOSINE-3'-5'-DIPHOSPHATE	
  'AAM':'', # ALPHA-ADENOSINE MONOPHOSPHATE	
  'AD3':'', # 3-DEAZA-ADENOSINE	
  'ADN':'', # ADENOSINE
  'ADS':'', # ADENOSINE-5'-(DITHIO)PHOSPHATE	
  'ADW':'', # ADENOSINE-5'-DITUNGSTATE	
  'AP2':'', # PHOSPHOMETHYLPHOSPHONIC ACID ADENOSYL ESTER	
  'APC':'', # DIPHOSPHOMETHYLPHOSPHONIC ACID ADENOSYL ESTER	
  ' AS':'', # 2-DEOXY-ADENOSINE -5'-THIO-MONOPHOSPHATE	
  'AS':'', # 2-DEOXY-ADENOSINE -5'-THIO-MONOPHOSPHATE	
  'ABP':'', # 8-BROMOADENOSINE-5'-DIPHOSPHATE	
  'ACP':'', # PHOSPHOMETHYLPHOSPHONIC ACID ADENYLATE ESTER	
  'ACQ':'', # DIPHOSPHOMETHYLPHOSPHONIC ACID ADENYLATE ESTER	
  'AD9':'', # ADP METAVANADATE	
  'ADX':'', # ADENOSINE-5'-DIPHOSPHATE-GLUCOSE	
  'AGS':'', # PHOSPHOTHIOPHOSPHORIC ACID-ADENYLATE ESTER	
  'AN2':'', # AMP PHOSPHORAMIDATE	
  'AOV':'', # ADP ORTHOVANADATE	
  'AP7':'', # N1-PROTONATED ADENOSINE-5'-MONOPHOSPHATE	
  'APR':'', # ADENOSINE-5-DIPHOSPHORIBOSE	
  'APW':'', # {5'-O-[(R)-{[(S)-AMINO(HYDROXY-KAPPAO)PHOSPHORYL]OXY}(HYDROXY-KAPPAO)PHOSPHORYL]ADENOSINATO(2-)}MAGNESIUM	
  'AQH':'', # [(2R,3S,4R,5R)-5-(6-amino-9H-purin-9-yl)-3,4-dihydroxytetrahydrofuran-2-yl]methyl (2R,3R,4R,5R,6S)-6-[(1R)-1,2-dihydroxyethyl]-3,4,5-trihydroxytetrahydro-2H-pyran-2-yl dihydrogen diphosphate	
  'AQP':'', # ADENOSINE-5'-TETRAPHOSPHATE	
  'AT4':'', # 5'-O-[(R)-HYDROXY(THIOPHOSPHONOOXY)PHOSPHORYL]ADENOSINE	
  'ATF':'', # PHOSPHODIFLUOROMETHYLPHOSPHONIC ACID-ADENYLATE ESTER	
  'ATR':'', # 2'-MONOPHOSPHOADENOSINE-5'-DIPHOSPHATE	
  'ATS':'', # GAMMA-ARSONO-BETA, GAMMA-METHYLENEADENOSINE-5'-DIPHOSPHATE	
  'AV2':'', # ADENOSINE-5'-DIPHOSPHATE-2',3'-VANADATE	
  'AVC':'', # ADENOSINE-5'-MONOPHOSPHATE-2',3'-VANADATE	
  'A12':'', # PHOSPHOMETHYLPHOSPHONIC ACID ADENOSYL ESTER	
  'AU1':'', # 5'-O-[(R)-hydroxy(phosphonoamino)phosphoryl]adenosine	
  'DAT':'', # 2'-DEOXYADENOSINE-5'-DIPHOSPHATE	
  ' DA':'', # 2'-DEOXYADENOSINE-5'-MONOPHOSPHATE	
  'DA':'', # 2'-DEOXYADENOSINE-5'-MONOPHOSPHATE	
  'CTP':'', # CYTIDINE-5'-TRIPHOSPHATE	
  'JZU':'', # 5'-deoxy-5'-(sulfamoylamino)adenosine	
  'LA8':'', # L-ADENOSINE-5'-DIPHOSPHATE	
  'M33':'', # 5'-O-[(S)-hydroxy{[(S)-hydroxy(methyl)phosphoryl]oxy}phosphoryl]adenosine	
}

def Nucleotides(name):
  if name in nc_dict:
    return True
  else:
    return False