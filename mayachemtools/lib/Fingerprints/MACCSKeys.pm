package Fingerprints::MACCSKeys;
#
# File: MACCSKeys.pm
# Author: Manish Sud <msud@san.rr.com>
#
# Copyright (C) 2018 Manish Sud. All rights reserved.
#
# This file is part of MayaChemTools.
#
# MayaChemTools is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# MayaChemTools is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability of fitness
# for a particular purpose.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MayaChemTools; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation Inc., 59 Temple Place, Suite 330,
# Boston, MA, 02111-1307, USA.
#

use strict;
use Carp;
use Exporter;
use Fingerprints::Fingerprints;
use TextUtil ();
use Molecule;
use PeriodicTable;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Fingerprints::Fingerprints Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

%EXPORT_TAGS = (all  => [@EXPORT, @EXPORT_OK]);

# Setup class variables...
my($ClassName);
_InitializeClass();

# Overload Perl functions...
use overload '""' => 'StringifyMACCSKeys';

# Class constructor...
sub new {
  my($Class, %NamesAndValues) = @_;

  # Initialize object...
  my $This = $Class->SUPER::new();
  bless $This, ref($Class) || $Class;
  $This->_InitializeMACCSKeys();

  $This->_InitializeMACCSKeysProperties(%NamesAndValues);

  return $This;
}

# Initialize object data...
#
sub _InitializeMACCSKeys {
  my($This) = @_;

  # Type of fingerprint to generate:
  #
  # MACCSKeyBits - A bit vector indicating presence/absence of keys
  # MACCSKeyCount - A vector containing count of keys
  #
  $This->{Type} = '';
  $This->{KeyBits} = '';

  # Size of key set: 166 or 322...
  $This->{Size} = '';
}

# Initialize class ...
sub _InitializeClass {
  #Class name...
  $ClassName = __PACKAGE__;
}

# Initialize object properties....
sub _InitializeMACCSKeysProperties {
  my($This, %NamesAndValues) = @_;

  my($Name, $Value, $MethodName);
  while (($Name, $Value) = each  %NamesAndValues) {
    $MethodName = "Set${Name}";
    $This->$MethodName($Value);
  }

  # Make sure molecule object was specified...
  if (!exists $NamesAndValues{Molecule}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying molecule...";
  }

  # Make sure type and size were specified...
  if (!exists $NamesAndValues{Type}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying type...";
  }
  if (!exists $NamesAndValues{Size}) {
    croak "Error: ${ClassName}->New: Object can't be instantiated without specifying size...";
  }

  # Make sure approriate size is specified...
  if ($NamesAndValues{Size} !~ /^(166|322)$/) {
    croak "Error: ${ClassName}->New: The current release of MayaChemTools doesn't support MDL MACCS $NamesAndValues{Size} keys...";
  }

  if ($This->{Type} =~ /^MACCSKeyBits$/i) {
    $This->_InitializeMACCSKeyBits();
  }
  elsif ($This->{Type} =~ /^MACCSKeyCount$/i) {
    $This->_InitializeMACCSKeyCounts();
  }
  else {
    croak "Error: ${ClassName}->_InitializeMACCSKeysProperties: Unknown MACCS keys type: $This->{Type}; Supported type keys: MACCSKeyBits or MACCSKeyCount......";
  }

  return $This;
}

# Initialize MACCS key bits...
#
sub _InitializeMACCSKeyBits {
  my($This) = @_;

  $This->{KeyBits} = 1;

  # Vector type...
  $This->{VectorType} = 'FingerprintsBitVector';

  $This->_InitializeFingerprintsBitVector();

  return $This;
}

# Initialize MACCS key counts...
#
sub _InitializeMACCSKeyCounts {
  my($This) = @_;

  $This->{KeyBits} = 0;

  # Vector type and type of values...
  $This->{VectorType} = 'FingerprintsVector';
  $This->{FingerprintsVectorType} = 'OrderedNumericalValues';

  $This->_InitializeFingerprintsVector();

  # Initialize values to zero...
  my(@Values);
  @Values = (0) x $This->{Size};
  $This->{FingerprintsVector}->AddValues(\@Values);

  return $This;
}

# Set type...
#
sub SetType {
  my($This, $Type) = @_;

  if ($This->{Type}) {
    croak "Error: ${ClassName}->SetType: Can't change type:  It's already set...";
  }

  if ($Type =~ /^MACCSKeyBits$/i) {
    $This->{Type} = 'MACCSKeyBits';;
    $This->{KeyBits} = 1;
  }
  elsif ($Type =~ /^MACCSKeyCount$/i) {
    $This->{Type} = 'MACCSKeyCount';;
    $This->{KeyBits} = 0;
  }
  else {
    croak "Error: ${ClassName}->SetType: Unknown type MACCS keys: $Type; Supported type keys: MACCSKeyBits or MACCSKeyCount...";
  }
  return $This;
}

# Set size...
#
sub SetSize {
  my($This, $Value) = @_;

  if ($This->{Size}) {
    croak "Error: ${ClassName}->SetSize: Can't change size:  It's already set...";
  }
  if (!TextUtil::IsPositiveInteger($Value)) {
    croak "Error: ${ClassName}->SetSize: Size value, $Value, is not valid:  It must be a positive integer...";
  }
  if ($Value !~ /^(166|322)/i) {
    croak "Error: ${ClassName}->Size: The current release of MayaChemTools doesn't support MDL MACCS $Value keys...";
  }
  $This->{Size} = $Value;

  return $This;
}

# Generate description...
#
sub GetDescription {
  my($This) = @_;

  # Is description explicity set?
  if (exists $This->{Description}) {
    return $This->{Description};
  }

  return "$This->{Type}";
}

# Generate MDL MACCS keys..
#
sub GenerateMACCSKeys {
  my($This) = @_;

  # Cache appropriate molecule data...
  $This->_SetupMoleculeDataCache();

  if ($This->{Size} == 166) {
    $This->_GenerateMACCS166Keys();
  }
  elsif ($This->{Size} == 322) {
    $This->_GenerateMACCS322Keys();
  }
  else {
    croak "Error: ${ClassName}->GenerateMACCSKeys: The current release of MayaChemTools doesn't support MDL MACCS  $This->{Size} keys...";
  }

  $This->{FingerprintsGenerated} = 1;

  # Clear cached molecule data...
  $This->_ClearMoleculeDataCache();

  return $This;
}

# Setup GenerateFingerprints method in order to be consistent with all other
# fingerprints classes implemented in the current release of MayaChemTools...
#
sub GenerateFingerprints {
  my($This) = @_;

  return $This->GenerateMACCSKeys();
}

# Generate MDL MACCS 166 keys...
#
# Information on the 166 keys [ Ref. 45-47 ]:
#
# Atom symbols:
#
# A : Any valid perodic table element symbol
# Q  : Hetro atoms; any non-C or non-H atom
# X  : Halogens; F, Cl, Br, I
# Z  : Others; other than H, C, N, O, Si, P, S, F, Cl, Br, I
#
# Bond types:
#
# -  : Single
# =  : Double
# T  : Triple
# #  : Triple
# ~  : Single or double query bond
# %  : An aromatic query bond
#
# None : Any bond type; no explict bond specified
#
# $  : Ring bond; $ before a bond type specifies ring bond
# !  : Chain or non-ring bond; ! before a bond type specifies chain bond
#
# @  : A ring linkage and the number following it specifies the
#      atoms position in the line, thus @1 means linked back to the first atom in
#      the list.
#
# Aromatic: Kekule or Arom5
#
# Kekule: Bonds in 6-membered rings with alternalte single/double bonds or perimeter
#         bonds
#
# Arom5: Bonds in 5-membered rings with two double bonds and a hetro atom at
#        the apex of the ring.
#
# Index Key Description
# 1	ISOTOPE
# 2	103 < ATOMIC NO. < 256
# 3	GROUP IVA,VA,VIA PERIODS 4-6 (Ge...)
# 4	ACTINIDE
# 5	GROUP IIIB,IVB (Sc...)
# 6	LANTHANIDE
# 7	GROUP VB,VIB,VIIB (V...)
# 8	QAAA@1
# 9	GROUP VIII (Fe...)
# 10	GROUP IIA (ALKALINE EARTH)
# 11	4M RING
# 12	GROUP IB,IIB (Cu...)
# 13	ON(C)C
# 14	S-S
# 15	OC(O)O
# 16	QAA@1
# 17	CTC
# 18	GROUP IIIA (B...)
# 19	7M RING
# 20	SI
# 21	C=C(Q)Q
# 22	3M RING
# 23	NC(O)O
# 24	N-O
# 25	NC(N)N
# 26	C$=C($A)$A
# 27	I
# 28	QCH2Q
# 29	P
# 30	CQ(C)(C)A
# 31	QX
# 32	CSN
# 33	NS
# 34	CH2=A
# 35	GROUP IA (ALKALI METAL)
# 36	S HETEROCYCLE
# 37	NC(O)N
# 38	NC(C)N
# 39	OS(O)O
# 40	S-O
# 41	CTN
# 42	F
# 43	QHAQH
# 44	OTHER
# 45	C=CN
# 46	BR
# 47	SAN
# 48	OQ(O)O
# 49	CHARGE
# 50	C=C(C)C
# 51	CSO
# 52	NN
# 53	QHAAAQH
# 54	QHAAQH
# 55	OSO
# 56	ON(O)C
# 57	O HETEROCYCLE
# 58	QSQ
# 59	Snot%A%A
# 60	S=O
# 61	AS(A)A
# 62	A$A!A$A
# 63	N=O
# 64	A$A!S
# 65	C%N
# 66	CC(C)(C)A
# 67	QS
# 68	QHQH (&...)
# 69	QQH
# 70	QNQ
# 71	NO
# 72	OAAO
# 73	S=A
# 74	CH3ACH3
# 75	A!N$A
# 76	C=C(A)A
# 77	NAN
# 78	C=N
# 79	NAAN
# 80	NAAAN
# 81	SA(A)A
# 82	ACH2QH
# 83 	QAAAA@1
# 84	NH2
# 85	CN(C)C
# 86	CH2QCH2
# 87	X!A$A
# 88	S
# 89	OAAAO
# 90	QHAACH2A
# 91	QHAAACH2A
# 92	OC(N)C
# 93	QCH3
# 94	QN
# 95	NAAO
# 96	5M RING
# 97	NAAAO
# 98	QAAAAA@1
# 99	C=C
# 100	ACH2N
# 101	8M RING
# 102	QO
# 103	CL
# 104	QHACH2A
# 105	A$A($A)$A
# 106	QA(Q)Q
# 107	XA(A)A
# 108	CH3AAACH2A
# 109	ACH2O
# 110	NCO
# 111	NACH2A
# 112	AA(A)(A)A
# 113	Onot%A%A
# 114	CH3CH2A
# 115	CH3ACH2A
# 116	CH3AACH2A
# 117	NAO
# 118	ACH2CH2A > 1
# 119	N=A
# 120	HETEROCYCLIC ATOM > 1 (&...)
# 121	N HETEROCYCLE
# 122	AN(A)A
# 123	OCO
# 124	QQ
# 125	AROMATIC RING > 1
# 126	A!O!A
# 127	A$A!O > 1 (&...)
# 128	ACH2AAACH2A
# 129	ACH2AACH2A
# 130	QQ > 1 (&...)
# 131	QH > 1
# 132	OACH2A
# 133	A$A!N
# 134	X (HALOGEN)
# 135	Nnot%A%A
# 136	O=A > 1
# 137	HETEROCYCLE
# 138	QCH2A > 1 (&...)
# 139	OH
# 140	O > 3 (&...)
# 141	CH3 > 2 (&...)
# 142	N > 1
# 143	A$A!O
# 144	Anot%A%Anot%A
# 145	6M RING > 1
# 146	O > 2
# 147	ACH2CH2A
# 148	AQ(A)A
# 149	CH3 > 1
# 150	A!A$A!A
# 151	NH
# 152	OC(C)C
# 153	QCH2A
# 154	C=O
# 155	A!CH2!A
# 156	NA(A)A
# 157	C-O
# 158	C-N
# 159	O > 1
# 160	CH3
# 161	N
# 162	AROMATIC
# 163	6M RING
# 164	O
# 165	RING
# 166 	FRAGMENTS
#
sub _GenerateMACCS166Keys {
  my($This) = @_;
  my($KeyNum, $KeyIndex, $MethodName, $KeyValue, $SkipPosCheck);

  $SkipPosCheck = 1;

  # Generate and set key values...
  KEYNUM: for $KeyNum (1 .. 166) {
    $MethodName = "_Generate166KeySetKey${KeyNum}";
    $KeyValue = $This->$MethodName();

    if (!$KeyValue) {
      next KEYNUM;
    }
    $KeyIndex = $KeyNum - 1;
    if ($This->{KeyBits}) {
      $This->{FingerprintsBitVector}->SetBit($KeyIndex, $SkipPosCheck);
    }
    else {
      $This->{FingerprintsVector}->SetValue($KeyIndex, $KeyValue, $SkipPosCheck);
    }
  }

  # Add key labels for MACCSKeyCount...
  if (!$This->{KeyBits}) {
    $This->_SetMACCSKeyCountValueIDs();
  }

  return $This;
}

# Generate MDL MACCS 322 keys...
#
# MDL MACCS 322 key set is defined in tables 1, 2 and 3 by: Joseph L. Durant; Burton A. Leland;
# Douglas R. Henry; James G. Nourse. Reoptimization of MDL Keys for Use in Drug Discovery [ Ref. 46 ].
#
# Atom symbols:
#
# A : Any valid perodic table element symbol
# Q  : Hetro atoms; any non-C or non-H atom
# X  : Others; other than H, C, N, O, Si, P, S, F, Cl, Br, I
# Z is neither defined nor used
#
# Atom symbol, X, used for 322 keys [ Ref 46 ] doesn't refer to Halogens as it does for 166 keys. In
# order to keep the definition of 322 keys consistent with the published definitions, the symbol X is
# used to imply "others" atoms, but it's internally mapped to symbol X as defined for 166 keys
# during the generation of key values.
#
# The keys include:
#
# o 26 atom properties of type P, as listed in Table 1
# o 32 one-atom environments, as listed in Table 3
# o 264 atom-bond-atom combinations listed in Table 4
#
# Total number of keys in three tables: 322
#
# Removal of two rare properties in Table 1 number 21 and 22 results in a 320 keyset.
#
# Atom properties-based keys (26):
#
# Index Description
# 1     A(AAA) or AA(A)A - atom with at least three neighbors
# 2     Q - heteroatom
# 3     Anot%not-A - atom involved in one or more multiple bonds, not aromatic
# 4     A(AAAA) or AA(A)(A)A - atom with at least four neighbors
# 5     A(QQ) or QA(Q) - atom with at least two heteroatom neighbors
# 6     A(QQQ) or QA(Q)Q - atom with at least three heteroatom neighbors
# 7     QH - heteroatom with at least one hydrogen attached
# 8     CH2(AA) or ACH2A - carbon with at least two single bonds and at least two hydrogens attached
# 9     CH3(A) or ACH3 - carbon with at least one single bond and at least three hydrogens attached
# 10    Halogen
# 11    A(-A-A-A) or A-A(-A)-A - atom has at least three single bonds
# 12    AAAAAA@1 > 2 - atom is in at least two different six-membered rings
# 13    A($A$A$A) or A$A($A)$A - atom has more than two ring bonds
# 14    A$A!A$A - atom is at a ring/chain boundary. When a comparison is done
#       with another atom the path passes through the chain bond.
# 15    Anot%A%Anot%A - atom is at an aromatic/nonaromatic boundary. When a
#       comparison is done with another atom the path
#       passes through the aromatic bond.
# 16    A!A!A  - atom with more than one chain bond
# 17    A!A$A!A - atom is at a ring/chain boundary. When a comparison is done
#       with another atom the path passes through the ring bond.
# 18    A%Anot%A%A - atom is at an aromatic/nonaromatic boundary. When a
#       comparison is done with another atom the
#       path passes through the nonaromatic bond.
# 19    HETEROCYCLE - atom is a heteroatom in a ring.
# 20    rare properties: atom with five or more neighbors, atom in
#       four or more rings, or atom types other than
#       H, C, N, O, S, F, Cl, Br, or I
# 21    rare properties: atom has a charge, is an isotope, has two or
#       more multiple bonds, or has a triple bond.
# 22    N - nitrogen
# 23    S - sulfur
# 24    O - oxygen
# 25    A(AA)A(A)A(AA) - atom has two neighbors, each with three or more neighbors
#       (including the central atom).
# 26    CHACH2 - atom has two hydrocarbon (CH2) neighbors
#
#
# Atomic environments properties-based keys (32):
#
# Index Key Description
# 27    C(CC)
# 28    C(CCC)
# 29    C(CN)
# 30    C(CCN)
# 31    C(NN)
# 32    C(NNC)
# 33    C(NNN)
# 34    C(CO)
# 35    C(CCO)
# 36    C(NO)
# 37    C(NCO)
# 38    C(NNO)
# 39    C(OO)
# 40    C(COO)
# 41    C(NOO)
# 42    C(OOO)
# 43    Q(CC)
# 44    Q(CCC)
# 45    Q(CN)
# 46    Q(CCN)
# 47    Q(NN)
# 48    Q(CNN)
# 49    Q(NNN)
# 50    Q(CO)
# 51    Q(CCO)
# 52    Q(NO)
# 53    Q(CNO)
# 54    Q(NNO)
# 55    Q(OO)
# 56    Q(COO)
# 57    Q(NOO)
# 58    Q(OOO)
#
# Note: The first symbol is the central atom, with atoms bonded to the
# central atom listed in parentheses. Q is any non-C, non-H atom. If
# only two atoms are in parentheses, there is no implication concerning
# the other atoms bonded to the central atom.
#
# Atom-Bond-Atom properties-based keys: (264)
#
# Index Key Description
# 59    C-C
# 60    C-N
# 61    C-O
# 62    C-S
# 63    C-Cl
# 64    C-P
# 65    C-F
# 66    C-Br
# 67    C-Si
# 68    C-I
# 69    C-X
# 70    N-N
# 71    N-O
# 72    N-S
# 73    N-Cl
# 74    N-P
# 75    N-F
# 76    N-Br
# 77    N-Si
# 78    N-I
# 79    N-X
# 80    O-O
# 81    O-S
# 82    O-Cl
# 83    O-P
# 84    O-F
# 85    O-Br
# 86    O-Si
# 87    O-I
# 88    O-X
# 89    S-S
# 90    S-Cl
# 91    S-P
# 92    S-F
# 93    S-Br
# 94    S-Si
# 95    S-I
# 96    S-X
# 97    Cl-Cl
# 98    Cl-P
# 99    Cl-F
# 100   Cl-Br
# 101   Cl-Si
# 102   Cl-I
# 103   Cl-X
# 104   P-P
# 105   P-F
# 106   P-Br
# 107   P-Si
# 108   P-I
# 109   P-X
# 110   F-F
# 111   F-Br
# 112   F-Si
# 113   F-I
# 114   F-X
# 115   Br-Br
# 116   Br-Si
# 117   Br-I
# 118   Br-X
# 119   Si-Si
# 120   Si-I
# 121   Si-X
# 122   I-I
# 123   I-X
# 124   X-X
# 125   C=C
# 126   C=N
# 127   C=O
# 128   C=S
# 129   C=Cl
# 130   C=P
# 131   C=F
# 132   C=Br
# 133   C=Si
# 134   C=I
# 135   C=X
# 136   N=N
# 137   N=O
# 138   N=S
# 139   N=Cl
# 140   N=P
# 141   N=F
# 142   N=Br
# 143   N=Si
# 144   N=I
# 145   N=X
# 146   O=O
# 147   O=S
# 148   O=Cl
# 149   O=P
# 150   O=F
# 151   O=Br
# 152   O=Si
# 153   O=I
# 154   O=X
# 155   S=S
# 156   S=Cl
# 157   S=P
# 158   S=F
# 159   S=Br
# 160   S=Si
# 161   S=I
# 162   S=X
# 163   Cl=Cl
# 164   Cl=P
# 165   Cl=F
# 166   Cl=Br
# 167   Cl=Si
# 168   Cl=I
# 169   Cl=X
# 170   P=P
# 171   P=F
# 172   P=Br
# 173   P=Si
# 174   P=I
# 175   P=X
# 176   F=F
# 177   F=Br
# 178   F=Si
# 179   F=I
# 180   F=X
# 181   Br=Br
# 182   Br=Si
# 183   Br=I
# 184   Br=X
# 185   Si=Si
# 186   Si=I
# 187   Si=X
# 188   I=I
# 189   I=X
# 190   X=X
# 191   C#C
# 192   C#N
# 193   C#O
# 194   C#S
# 195   C#Cl
# 196   C#P
# 197   C#F
# 198   C#Br
# 199   C#Si
# 200   C#I
# 201   C#X
# 202   N#N
# 203   N#O
# 204   N#S
# 205   N#Cl
# 206   N#P
# 207   N#F
# 208   N#Br
# 209   N#Si
# 210   N#I
# 211   N#X
# 212   O#O
# 213   O#S
# 214   O#Cl
# 215   O#P
# 216   O#F
# 217   O#Br
# 218   O#Si
# 219   O#I
# 220   O#X
# 221   S#S
# 222   S#Cl
# 223   S#P
# 224   S#F
# 225   S#Br
# 226   S#Si
# 227   S#I
# 228   S#X
# 229   Cl#Cl
# 230   Cl#P
# 231   Cl#F
# 232   Cl#Br
# 233   Cl#Si
# 234   Cl#I
# 235   Cl#X
# 236   P#P
# 237   P#F
# 238   P#Br
# 239   P#Si
# 240   P#I
# 241   P#X
# 242   F#F
# 243   F#Br
# 244   F#Si
# 245   F#I
# 246   F#X
# 247   Br#Br
# 248   Br#Si
# 249   Br#I
# 250   Br#X
# 251   Si#Si
# 252   Si#I
# 253   Si#X
# 254   I#I
# 255   I#X
# 256   X#X
# 257   C$C
# 258   C$N
# 259   C$O
# 260   C$S
# 261   C$Cl
# 262   C$P
# 263   C$F
# 264   C$Br
# 265   C$Si
# 266   C$I
# 267   C$X
# 268   N$N
# 269   N$O
# 270   N$S
# 271   N$Cl
# 272   N$P
# 273   N$F
# 274   N$Br
# 275   N$Si
# 276   N$I
# 277   N$X
# 278   O$O
# 279   O$S
# 280   O$Cl
# 281   O$P
# 282   O$F
# 283   O$Br
# 284   O$Si
# 285   O$I
# 286   O$X
# 287   S$S
# 288   S$Cl
# 289   S$P
# 290   S$F
# 291   S$Br
# 292   S$Si
# 293   S$I
# 294   S$X
# 295   Cl$Cl
# 296   Cl$P
# 297   Cl$F
# 298   Cl$Br
# 299   Cl$Si
# 300   Cl$I
# 301   Cl$X
# 302   P$P
# 303   P$F
# 304   P$Br
# 305   P$Si
# 306   P$I
# 307   P$X
# 308   F$F
# 309   F$Br
# 310   F$Si
# 311   F$I
# 312   F$X
# 313   Br$Br
# 314   Br$Si
# 315   Br$I
# 316   Br$X
# 317   Si$Si
# 318   Si$I
# 319   Si$X
# 320   I$I
# 321   I$X
# 322   X$X
#
# Note: Instead of using '%' as rind bond as mentioned in the article [ Ref. 46 ], MayaChemTools
# used '$' as a symbol for ring bond to follow conventions used for MACCS 166 keys; the symbol '%'
# is used to indicate an aromatic query bond.
#
sub _GenerateMACCS322Keys {
  my($This) = @_;
  my($KeyNum, $KeyIndex, $MethodName, $KeyValue, $SkipPosCheck);

  $SkipPosCheck = 1;

  # Generate and set key values...
  KEYNUM: for $KeyNum (1 .. 322) {
    $MethodName = "_Generate322KeySetKey${KeyNum}";
    $KeyValue = $This->$MethodName();

    if (!$KeyValue) {
      next KEYNUM;
    }
    $KeyIndex = $KeyNum - 1;
    if ($This->{KeyBits}) {
      $This->{FingerprintsBitVector}->SetBit($KeyIndex, $SkipPosCheck);
    }
    else {
      $This->{FingerprintsVector}->SetValue($KeyIndex, $KeyValue, $SkipPosCheck);
    }
  }

  # Add key labels for MACCSKeyCount...
  if (!$This->{KeyBits}) {
    $This->_SetMACCSKeyCountValueIDs();
  }
  return $This;
}

# Set MACCS key count value IDs for fingerprint vector. The value IDs labels format
# is: Key<KeyNum>.
#
# By default, no value IDs are set for fingerprint vector values.
#
sub _SetMACCSKeyCountValueIDs {
  my($This) = @_;

  if (!$This->{FingerprintsVector}) {
    return;
  }
  my(@ValueIDs);

  @ValueIDs = map { "Key$_"; } (1 .. $This->{Size});
  $This->{FingerprintsVector}->AddValueIDs(\@ValueIDs);

  return $This;
}

##################################
#
#  Implementation of MDL MACCS 166 keys...
#
##################################

# Generate key 1 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 1 description: ISOTOPE
#
sub _Generate166KeySetKey1 {
  my($This) = @_;
  my($Atom, $KeyValue);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ($Atom->IsIsotope()) {
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
      $KeyValue++;
    }
  }
  return $KeyValue;
}

# Generate key 2 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 2 description: 103 < ATOMIC NO. < 256
#
sub _Generate166KeySetKey2 {
  my($This) = @_;
  my($Atom, $AtomicNumber, $KeyValue);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber > 103 && $AtomicNumber < 256) {
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
      $KeyValue++;
    }
  }
  return $KeyValue;
}

# Generate key 3 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 3 description: GROUP IVA,VA,VIA (GroupNumber: 14, 15, 16) PERIODS 4-6 (Ge...)
#
sub _Generate166KeySetKey3 {
  my($This) = @_;
  my($Atom, $KeyValue, $AtomicNumber, $GroupNumber, $PeriodNumber);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber) {
      $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);
      $PeriodNumber = PeriodicTable::GetElementPeriodNumber($AtomicNumber);
      if ($PeriodNumber =~ /^(4|5|6)$/ && $GroupNumber =~ /^(14|15|16)$/) {
	if ($This->{KeyBits}) {
	  $KeyValue = 1;
	  last ATOM;
	}
	$KeyValue++;
      }
    }
  }
  return $KeyValue;
}

# Generate key 4 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 4 description: ACTINIDE
#
sub _Generate166KeySetKey4 {
  my($This) = @_;
  my($Atom, $AtomicNumber, $KeyValue);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber >= 89 && $AtomicNumber <= 103) {
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
      $KeyValue++;
    }
  }
  return $KeyValue;
}

# Generate key 5 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 5 description: GROUP IIIB,IVB (Sc...)
#
sub _Generate166KeySetKey5 {
  my($This) = @_;
  my($Atom, $KeyValue, $AtomicNumber, $GroupNumber);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber) {
      $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);
      if ($GroupNumber =~ /^(3|4)$/) {
	if ($This->{KeyBits}) {
	  $KeyValue = 1;
	  last ATOM;
	}
	$KeyValue++;
      }
    }
  }
  return $KeyValue;
}

# Generate key 6 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 6 description: LANTHANIDE
#
sub _Generate166KeySetKey6 {
  my($This) = @_;
  my($Atom, $AtomicNumber, $KeyValue);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber >= 57 && $AtomicNumber <= 71) {
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
      $KeyValue++;
    }
  }
  return $KeyValue;
}

# Generate key 7 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 7 description: GROUP VB,VIB,VIIB (V...)
#
sub _Generate166KeySetKey7 {
  my($This) = @_;
  my($Atom, $KeyValue, $AtomicNumber, $GroupNumber);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber) {
      $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);
      if ($GroupNumber =~ /^(5|6|7)$/) {
	if ($This->{KeyBits}) {
	  $KeyValue = 1;
	  last ATOM;
	}
	$KeyValue++;
      }
    }
  }
  return $KeyValue;
}

# Generate key 8 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 8 description: QAAA@1
#
sub _Generate166KeySetKey8 {
  my($This) = @_;
  my($Atom, $KeyValue, $RingSize);

  $RingSize = 4;
  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ($This->_IsHeteroAtom($Atom) && $Atom->IsInRingOfSize($RingSize)) {
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
      $KeyValue++;
    }
  }
  return $KeyValue;
}

# Generate key 9 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 9 description: GROUP VIII (Fe...)
#
sub _Generate166KeySetKey9 {
  my($This) = @_;
  my($Atom, $KeyValue, $AtomicNumber, $GroupNumber);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber) {
      $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);
      if ($GroupNumber =~ /^(8|9|10)$/) {
	if ($This->{KeyBits}) {
	  $KeyValue = 1;
	  last ATOM;
	}
	$KeyValue++;
      }
    }
  }
  return $KeyValue;
}

# Generate key 10 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 10 description: GROUP IIA (ALKALINE EARTH)
#
sub _Generate166KeySetKey10 {
  my($This) = @_;
  my($Atom, $KeyValue, $AtomicNumber, $GroupNumber);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber) {
      $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);
      if ($GroupNumber =~ /^2$/) {
	if ($This->{KeyBits}) {
	  $KeyValue = 1;
	  last ATOM;
	}
	$KeyValue++;
      }
    }
  }
  return $KeyValue;
}

# Generate key 11 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 11 description: 4M RING
#
sub _Generate166KeySetKey11 {
  my($This) = @_;
  my($Molecule, $KeyValue, $RingSize, $NumOfRings);

  $RingSize = 4;
  $Molecule = $This->GetMolecule();
  $NumOfRings = $Molecule->GetNumOfRingsWithSize($RingSize);

  if ($This->{KeyBits}) {
    $KeyValue = $NumOfRings ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfRings;
  }
  return $KeyValue;
}

# Generate key 12 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 12 description: GROUP IB,IIB (Cu...)
#
sub _Generate166KeySetKey12 {
  my($This) = @_;
  my($Atom, $KeyValue, $AtomicNumber, $GroupNumber);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber) {
      $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);
      if ($GroupNumber =~ /^(11|12)$/) {
	if ($This->{KeyBits}) {
	  $KeyValue = 1;
	  last ATOM;
	}
	$KeyValue++;
      }
    }
  }
  return $KeyValue;
}

# Generate key 13 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 13 description: ON(C)C
#
sub _Generate166KeySetKey13 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'N';
  @NbrAtomSymbols = ('O', 'C', 'C');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 14 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 14 description: S-S
#
sub _Generate166KeySetKey14 {
  my($This) = @_;
  my($BondOrder) = 1;

  return $This->_DetectBondKeys('S', 'S', $BondOrder);
}

# Generate key 15 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 15 description: OC(O)O
#
sub _Generate166KeySetKey15 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('O', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 16 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 16 description: QAA@1
#
sub _Generate166KeySetKey16 {
  my($This) = @_;
  my($Atom, $KeyValue, $RingSize);

  $RingSize = 3;
  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ($This->_IsHeteroAtom($Atom) && $Atom->IsInRingOfSize($RingSize)) {
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
      $KeyValue++;
    }
  }
  return $KeyValue;
}

# Generate key 17 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 17 description: CTC
#
sub _Generate166KeySetKey17 {
  my($This) = @_;
  my($BondOrder) = 3;

  return $This->_DetectBondKeys('C', 'C', $BondOrder);
}

# Generate key 18 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 18 description: GROUP IIIA (B...)
#
sub _Generate166KeySetKey18 {
  my($This) = @_;
  my($Atom, $KeyValue, $AtomicNumber, $GroupNumber);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber) {
      $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);
      if ($GroupNumber =~ /^13$/) {
	if ($This->{KeyBits}) {
	  $KeyValue = 1;
	  last ATOM;
	}
	$KeyValue++;
      }
    }
  }
  return $KeyValue;
}

# Generate key 19 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 19 description: 7M RING
#
sub _Generate166KeySetKey19 {
  my($This) = @_;
  my($Molecule, $KeyValue, $RingSize, $NumOfRings);

  $RingSize = 7;
  $Molecule = $This->GetMolecule();
  $NumOfRings = $Molecule->GetNumOfRingsWithSize($RingSize);

  $KeyValue = 0;
  if ($NumOfRings) {
    $KeyValue = ($This->{KeyBits}) ? 1 : $NumOfRings;
  }
  return $KeyValue;
}

# Generate key 20 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 20 description: SI
#
sub _Generate166KeySetKey20 {
  my($This) = @_;

  return $This->_DetectAtomKeys('Si');
}

# Generate key 21 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 21 description: C=C(Q)Q
#
sub _Generate166KeySetKey21 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'Q', 'Q');
  @NbrBondSymbols = ('=', undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 22 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 22 description: 3M RING
#
sub _Generate166KeySetKey22 {
  my($This) = @_;
  my($Molecule, $KeyValue, $RingSize, $NumOfRings);

  $RingSize = 3;
  $Molecule = $This->GetMolecule();
  $NumOfRings = $Molecule->GetNumOfRingsWithSize($RingSize);

  if ($This->{KeyBits}) {
    $KeyValue = $NumOfRings ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfRings;
  }
  return $KeyValue;
}

# Generate key 23 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 23 description: NC(O)O
#
sub _Generate166KeySetKey23 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 24 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 24 description: N-O
#
sub _Generate166KeySetKey24 {
  my($This) = @_;
  my($BondOrder) = 1;

  return $This->_DetectBondKeys('N', 'O', $BondOrder);
}

# Generate key 25 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 25 description: NC(N)N
#
sub _Generate166KeySetKey25 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'N', 'N');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 26 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 26 description: C$=C($A)$A
#
sub _Generate166KeySetKey26 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'A', 'A');
  @NbrBondSymbols = ('$=', '$', '$');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 27 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 27 description: I
#
sub _Generate166KeySetKey27 {
  my($This) = @_;

  return $This->_DetectAtomKeys('I');
}

# Generate key 28 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 28 description: QCH2Q
#
sub _Generate166KeySetKey28 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('Q', 'Q');
  @NbrBondSymbols = (undef, undef);
  $MinKeyCount = undef;
  $CentralAtomMinHydrogenCount = 2;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate key 29 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 29 description: P
#
sub _Generate166KeySetKey29 {
  my($This) = @_;

  return $This->_DetectAtomKeys('P');
}

# Generate key 30 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 30 description: CQ(C)(C)A
#
sub _Generate166KeySetKey30 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'C', 'C', 'A');
  @NbrBondSymbols = (undef, undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 31 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 31 description: QX
#
sub _Generate166KeySetKey31 {
  my($This) = @_;

  return $This->_DetectBondKeys('Q', 'X');
}

# Generate key 32 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 32 description: CSN
#
sub _Generate166KeySetKey32 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'S';
  @NbrAtomSymbols = ('C', 'N');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 33 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 33 description: NS
#
sub _Generate166KeySetKey33 {
  my($This) = @_;

  return $This->_DetectBondKeys('N', 'S');
}

# Generate key 34 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 34 description: CH2=A
#
sub _Generate166KeySetKey34 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('A');
  @NbrBondSymbols = ('=');
  $MinKeyCount = undef;
  $CentralAtomMinHydrogenCount = 2;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate key 35 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 35 description: GROUP IA (ALKALI METAL)
#
sub _Generate166KeySetKey35 {
  my($This) = @_;
  my($Atom, $KeyValue, $AtomicNumber, $GroupNumber);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    $AtomicNumber = $Atom->GetAtomicNumber();
    if ($AtomicNumber) {
      $GroupNumber = PeriodicTable::GetElementGroupNumber($AtomicNumber);
      if ($GroupNumber =~ /^1$/) {
	if ($This->{KeyBits}) {
	  $KeyValue = 1;
	  last ATOM;
	}
	$KeyValue++;
      }
    }
  }
  return $KeyValue;
}

# Generate key 36 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 36 description: S HETEROCYCLE
#
sub _Generate166KeySetKey36 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing) = (1, 1);

  return $This->_DetectAtomKeys('S', $MinKeyCount, $IsInRing);
}

# Generate key 37 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 37 description: NC(O)N
#
sub _Generate166KeySetKey37 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'O', 'N');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 38 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 38 description: NC(C)N
#
sub _Generate166KeySetKey38 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'C', 'N');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 39 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 39 description: OS(O)O
#
sub _Generate166KeySetKey39 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'S';
  @NbrAtomSymbols = ('O', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 40 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 40 description: S-O
#
sub _Generate166KeySetKey40 {
  my($This) = @_;
  my($BondOrder) = 1;

  return $This->_DetectBondKeys('S', 'O', $BondOrder);
}

# Generate key 41 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 41 description: CTN
#
sub _Generate166KeySetKey41 {
  my($This) = @_;
  my($BondOrder) = 3;

  return $This->_DetectBondKeys('C', 'N', $BondOrder);
}

# Generate key 42 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 42 description: F
#
sub _Generate166KeySetKey42 {
  my($This) = @_;

  return $This->_DetectAtomKeys('F');
}

# Generate key 43 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 43 description: QHAQH
#
sub _Generate166KeySetKey43 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols, @NbrAtomMinHydrogenCount);

  $CentralAtomSymbol = 'A';
  $CentralAtomMinHydrogenCount = undef;

  @NbrAtomSymbols = ('Q', 'Q');
  @NbrBondSymbols = (undef, undef);
  @NbrAtomMinHydrogenCount = (1, 1);

  $MinKeyCount = undef;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount, \@NbrAtomMinHydrogenCount);
}

# Generate key 44 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 44 description: OTHER
#
sub _Generate166KeySetKey44 {
  my($This) = @_;

  return $This->_DetectAtomKeys('Z');
}

# Generate key 45 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 45 description: C=CN
#
sub _Generate166KeySetKey45 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'N');
  @NbrBondSymbols = ('=', undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 46 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 46 description: BR
#
sub _Generate166KeySetKey46 {
  my($This) = @_;

  return $This->_DetectAtomKeys('Br');
}

# Generate key 47 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 47 description: SAN
#
sub _Generate166KeySetKey47 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('S', 'N');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 48 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 48 description: OQ(O)O
#
sub _Generate166KeySetKey48 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('O', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 49 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 49 description: CHARGE
#
sub _Generate166KeySetKey49 {
  my($This) = @_;
  my($Molecule, $KeyValue);

  $Molecule = $This->GetMolecule();
  $KeyValue = $Molecule->GetFormalCharge() ? 1 : 0;

  return $KeyValue;
}

# Generate key 50 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 50 description: C=C(C)C
#
sub _Generate166KeySetKey50 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'C', 'C');
  @NbrBondSymbols = ('=', undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 51 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 51 description: CSO
#
sub _Generate166KeySetKey51 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'S';
  @NbrAtomSymbols = ('C', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 52 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 52 description: NN
#
sub _Generate166KeySetKey52 {
  my($This) = @_;

  return $This->_DetectBondKeys('N', 'N');
}

# Generate key 53 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 53 description: QHAAAQH
#
sub _Generate166KeySetKey53 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('Q', 'A', 'A', 'A', 'Q');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (1, undef, undef, undef, 1);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 54 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 54 description: QHAAQH
#
sub _Generate166KeySetKey54 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('Q', 'A', 'A', 'Q');
  @CentralAtomsBondSymbols = (undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (1, undef, undef, 1);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 55 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 55 description: OSO
#
sub _Generate166KeySetKey55 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'S';
  @NbrAtomSymbols = ('O', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 56 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 56 description: ON(O)C
#
sub _Generate166KeySetKey56 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'N';
  @NbrAtomSymbols = ('O', 'O', 'C');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 57 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 57 description: O HETEROCYCLE
#
sub _Generate166KeySetKey57 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing) = (undef, 1);

  return $This->_DetectAtomKeys('O', $MinKeyCount, $IsInRing);
}

# Generate key 58 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 58 description: QSQ
#
sub _Generate166KeySetKey58 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'S';
  @NbrAtomSymbols = ('Q', 'Q');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 59 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 59 description: Snot%A%A
#
sub _Generate166KeySetKey59 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('S', 'A');
  @NbrBondSymbols = ('not%', '%');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 60 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 60 description: S=O
#
sub _Generate166KeySetKey60 {
  my($This) = @_;
  my($BondOrder) = 2;

  return $This->_DetectBondKeys('S', 'O', $BondOrder);
}

# Generate key 61 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 61 description: AS(A)A
#
sub _Generate166KeySetKey61 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'S';
  @NbrAtomSymbols = ('A', 'A', 'A');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 62 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 62 description: A$A!A$A
#
sub _Generate166KeySetKey62 {
  my($This) = @_;
  my($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, @NbrAtomsSymbols, @NbrAtomsBondSymbols);

  ($BondAtomSymbol1, $BondAtomSymbol2) = ('A', 'A');
  $BondSymbol = '!';

  @NbrAtomsSymbols = (['A'], ['A']);
  @NbrAtomsBondSymbols = (['$'], ['$']);
  return $This->_DetectBondNeighborhoodKeys($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, \@NbrAtomsSymbols, \@NbrAtomsBondSymbols);
}

# Generate key 63 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 63 description: N=O
#
sub _Generate166KeySetKey63 {
  my($This) = @_;
  my($BondOrder) = 2;

  return $This->_DetectBondKeys('N', 'O', $BondOrder);
}

# Generate key 64 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 64 description: A$A!S
#
sub _Generate166KeySetKey64 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'S');
  @NbrBondSymbols = ('$', '!');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 65 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 65 description: C%N
#
sub _Generate166KeySetKey65 {
  my($This) = @_;
  my($BondSymbol) = '%';

  return $This->_DetectBondKeys('C', 'N', $BondSymbol);
}

# Generate key 66 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 66 description: CC(C)(C)A
#
sub _Generate166KeySetKey66 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'C', 'C', 'A');
  @NbrBondSymbols = (undef, undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 67 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 67 description: QS
#
sub _Generate166KeySetKey67 {
  my($This) = @_;

  return $This->_DetectBondKeys('Q', 'S');
}

# Generate key 68 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 68 description: QHQH (&...)
#
sub _Generate166KeySetKey68 {
  my($This) = @_;
  my($AtomSymbol1, $AtomSymbol2, $BondSymbol) = ('Q', 'Q', undef);
  my($MinKeyCount) = undef;
  my($Atom1MinHydrogenCount, $Atom2MinHydrogenCount) = (1, 1);

  return $This->_DetectBondKeys($AtomSymbol1, $AtomSymbol2, $BondSymbol, $MinKeyCount, $Atom1MinHydrogenCount, $Atom2MinHydrogenCount);
}

# Generate key 69 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 69 description: QQH
#
sub _Generate166KeySetKey69 {
  my($This) = @_;
  my($AtomSymbol1, $AtomSymbol2, $BondSymbol) = ('Q', 'Q', undef);
  my($MinKeyCount) = undef;
  my($Atom1MinHydrogenCount, $Atom2MinHydrogenCount) = (undef, 1);

  return $This->_DetectBondKeys($AtomSymbol1, $AtomSymbol2, $BondSymbol, $MinKeyCount, $Atom1MinHydrogenCount, $Atom2MinHydrogenCount);
}

# Generate key 70 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 70 description: QNQ
#
sub _Generate166KeySetKey70 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'N';
  @NbrAtomSymbols = ('Q', 'Q');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 71 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 71 description: NO
#
sub _Generate166KeySetKey71 {
  my($This) = @_;

  return $This->_DetectBondKeys('N', 'O');
}

# Generate key 72 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 72 description: OAAO
#
sub _Generate166KeySetKey72 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols);

  @CentralAtomsSymbols = ('O', 'A', 'A', 'O');
  @CentralAtomsBondSymbols = (undef, undef, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols);
}

# Generate key 73 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 73 description: S=A
#
sub _Generate166KeySetKey73 {
  my($This) = @_;
  my($BondOrder) = 2;

  return $This->_DetectBondKeys('S', 'A', $BondOrder);
}

# Generate key 74 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 74 description: CH3ACH3
#
sub _Generate166KeySetKey74 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols, @NbrAtomMinHydrogenCount);

  $CentralAtomSymbol = 'A';
  $CentralAtomMinHydrogenCount = undef;

  @NbrAtomSymbols = ('C', 'C');
  @NbrBondSymbols = (undef, undef);
  @NbrAtomMinHydrogenCount = (3, 3);

  $MinKeyCount = undef;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount, \@NbrAtomMinHydrogenCount);
}

# Generate key 75 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 75 description: A!N$A
#
sub _Generate166KeySetKey75 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'N';
  @NbrAtomSymbols = ('A', 'A');
  @NbrBondSymbols = ('!', '$');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 76 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 76 description: C=C(A)A
#
sub _Generate166KeySetKey76 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'A', 'A');
  @NbrBondSymbols = ('=', undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 77 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 77 description: NAN
#
sub _Generate166KeySetKey77 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('N', 'N');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 78 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 78 description: C=N
#
sub _Generate166KeySetKey78 {
  my($This) = @_;
  my($BondOrder) = 2;

  return $This->_DetectBondKeys('C', 'N', $BondOrder);
}

# Generate key 79 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 79 description: NAAN
#
sub _Generate166KeySetKey79 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols);

  @CentralAtomsSymbols = ('N', 'A', 'A', 'N');
  @CentralAtomsBondSymbols = (undef, undef, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols);
}

# Generate key 80 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 80 description: NAAAN
#
sub _Generate166KeySetKey80 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols);

  @CentralAtomsSymbols = ('N', 'A', 'A', 'A', 'N');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols);
}

# Generate key 81 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 81 description: SA(A)A
#
sub _Generate166KeySetKey81 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('S', 'A', 'A');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 82 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 82 description: ACH2QH
#
sub _Generate166KeySetKey82 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols, @NbrAtomMinHydrogenCount);

  $CentralAtomSymbol = 'C';
  $CentralAtomMinHydrogenCount = 2;

  @NbrAtomSymbols = ('A', 'Q');
  @NbrBondSymbols = (undef, undef);
  @NbrAtomMinHydrogenCount = (undef, 1);

  $MinKeyCount = undef;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount, \@NbrAtomMinHydrogenCount);
}

# Generate key 83 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 83 description: QAAAA@1
#
sub _Generate166KeySetKey83 {
  my($This) = @_;
  my($Atom, $KeyValue, $RingSize);

  $RingSize = 5;
  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ($This->_IsHeteroAtom($Atom) && $Atom->IsInRingOfSize($RingSize)) {
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
      $KeyValue++;
    }
  }
  return $KeyValue;
}

# Generate key 84 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 84 description: NH2
#
sub _Generate166KeySetKey84 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing, $MinHydrogenCount) = (undef, undef, 2);

  return $This->_DetectAtomKeys('N', $MinKeyCount, $IsInRing, $MinHydrogenCount);
}

# Generate key 85 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 85 description: CN(C)C
#
sub _Generate166KeySetKey85 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'N';
  @NbrAtomSymbols = ('C', 'C', 'C');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 86 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 86 description: CH2QCH2
#
sub _Generate166KeySetKey86 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols, @NbrAtomMinHydrogenCount);

  $CentralAtomSymbol = 'Q';
  $CentralAtomMinHydrogenCount = undef;

  @NbrAtomSymbols = ('C', 'C');
  @NbrBondSymbols = (undef, undef);
  @NbrAtomMinHydrogenCount = (2, 2);

  $MinKeyCount = undef;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount, \@NbrAtomMinHydrogenCount);
}

# Generate key 87 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 87 description: X!A$A
#
sub _Generate166KeySetKey87 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('X', 'A');
  @NbrBondSymbols = ('!', '$');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 88 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 88 description: S
#
sub _Generate166KeySetKey88 {
  my($This) = @_;

  return $This->_DetectAtomKeys('S');
}

# Generate key 89 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 89 description: OAAAO
#
sub _Generate166KeySetKey89 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols);

  @CentralAtomsSymbols = ('O', 'A', 'A', 'A', 'O');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols);
}

# Generate key 90 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 90 description: QHAACH2A
#
sub _Generate166KeySetKey90 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('Q', 'A', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (1, undef, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 91 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 91 description: QHAAACH2A
#
sub _Generate166KeySetKey91 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('Q', 'A', 'A', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (1, undef, undef, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 92 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 92 description: OC(N)C
#
sub _Generate166KeySetKey92 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('O', 'N', 'C');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 93 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 93 description: QCH3
#
sub _Generate166KeySetKey93 {
  my($This) = @_;
  my($AtomSymbol1, $AtomSymbol2, $BondSymbol) = ('Q', 'C', undef);
  my($MinKeyCount) = undef;
  my($Atom1MinHydrogenCount, $Atom2MinHydrogenCount) = (undef, 3);

  return $This->_DetectBondKeys($AtomSymbol1, $AtomSymbol2, $BondSymbol, $MinKeyCount, $Atom1MinHydrogenCount, $Atom2MinHydrogenCount);
}

# Generate key 94 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 94 description: QN
#
sub _Generate166KeySetKey94 {
  my($This) = @_;

  return $This->_DetectBondKeys('Q', 'N');
}

# Generate key 95 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 95 description: NAAO
#
sub _Generate166KeySetKey95 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols);

  @CentralAtomsSymbols = ('N', 'A', 'A', 'O');
  @CentralAtomsBondSymbols = (undef, undef, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols);
}

# Generate key 96 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 96 description: 5M RING
#
sub _Generate166KeySetKey96 {
  my($This) = @_;
  my($Molecule, $KeyValue, $RingSize, $NumOfRings);

  $RingSize = 5;
  $Molecule = $This->GetMolecule();
  $NumOfRings = $Molecule->GetNumOfRingsWithSize($RingSize);

  if ($This->{KeyBits}) {
    $KeyValue = $NumOfRings ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfRings;
  }
  return $KeyValue;
}

# Generate key 97 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 97 description: NAAAO
#
sub _Generate166KeySetKey97 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols);

  @CentralAtomsSymbols = ('N', 'A', 'A', 'A', 'O');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols);
}

# Generate key 98 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 98 description: QAAAAA@1
#
sub _Generate166KeySetKey98 {
  my($This) = @_;
  my($Atom, $KeyValue, $RingSize);

  $RingSize = 6;
  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ($This->_IsHeteroAtom($Atom) && $Atom->IsInRingOfSize($RingSize)) {
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
      $KeyValue++;
    }
  }
  return $KeyValue;
}

# Generate key 99 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 99 description: C=C
#
sub _Generate166KeySetKey99 {
  my($This) = @_;
  my($BondOrder) = 2;

  return $This->_DetectBondKeys('C', 'C', $BondOrder);
}

# Generate key 100 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 100 description: ACH2N
#
sub _Generate166KeySetKey100 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  $CentralAtomMinHydrogenCount = 2;

  @NbrAtomSymbols = ('A', 'N');
  @NbrBondSymbols = (undef, undef);

  $MinKeyCount = undef;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate key 101 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 101 description: 8M RING
#
sub _Generate166KeySetKey101 {
  my($This) = @_;
  my($Molecule, $KeyValue, $RingSize, $NumOfRings);

  $RingSize = 8;
  $Molecule = $This->GetMolecule();
  $NumOfRings = $Molecule->GetNumOfRingsWithSize($RingSize);

  if ($This->{KeyBits}) {
    $KeyValue = $NumOfRings ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfRings;
  }
  return $KeyValue;
}

# Generate key 102 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 102 description: QO
#
sub _Generate166KeySetKey102 {
  my($This) = @_;

  return $This->_DetectBondKeys('Q', 'O');
}

# Generate key 103 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 103 description: CL
#
sub _Generate166KeySetKey103 {
  my($This) = @_;

  return $This->_DetectAtomKeys('Cl');
}

# Generate key 104 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 104 description: QHACH2A
#
sub _Generate166KeySetKey104 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('Q', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (1, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 105 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 105 description: A$A($A)$A
#
sub _Generate166KeySetKey105 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'A', 'A');
  @NbrBondSymbols = ('$', '$', '$');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 106 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 106 description: QA(Q)Q
#
sub _Generate166KeySetKey106 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('Q', 'Q', 'Q');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 107 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 107 description: XA(A)A
#
sub _Generate166KeySetKey107 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('X', 'A', 'A');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 108 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 108 description: CH3AAACH2A
#
sub _Generate166KeySetKey108 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('C', 'A', 'A', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (3, undef, undef, undef, 1, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 109 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 109 description: ACH2O
#
sub _Generate166KeySetKey109 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  $CentralAtomMinHydrogenCount = 2;

  @NbrAtomSymbols = ('A', 'O');
  @NbrBondSymbols = (undef, undef);

  $MinKeyCount = undef;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate key 110 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 110 description: NCO
#
sub _Generate166KeySetKey110 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 111 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 111 description: NACH2A
#
sub _Generate166KeySetKey111 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('N', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (undef, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 112 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 112 description: AA(A)(A)A
#
sub _Generate166KeySetKey112 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'A', 'A', 'A');
  @NbrBondSymbols = (undef, undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 113 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 113 description: Onot%A%A
#
sub _Generate166KeySetKey113 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('O', 'A');
  @NbrBondSymbols = ('not%', '%');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 114 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 114 description: CH3CH2A
#
sub _Generate166KeySetKey114 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols, @NbrAtomMinHydrogenCount);

  $CentralAtomSymbol = 'C';
  $CentralAtomMinHydrogenCount = 2;

  @NbrAtomSymbols = ('C', 'A');
  @NbrBondSymbols = (undef, undef);
  @NbrAtomMinHydrogenCount = (3, undef);

  $MinKeyCount = undef;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount, \@NbrAtomMinHydrogenCount);
}

# Generate key 115 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 115 description: CH3ACH2A
#
sub _Generate166KeySetKey115 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('C', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (3, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 116 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 116 description: CH3AACH2A
#
sub _Generate166KeySetKey116 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('C', 'A', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (3, undef, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 117 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 117 description: NAO
#
sub _Generate166KeySetKey117 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('N', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 118 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 118 description: ACH2CH2A > 1
#
sub _Generate166KeySetKey118 {
  my($This) = @_;
  my($MinKeyCount, @CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  $MinKeyCount = 2;
  @CentralAtomsSymbols = ('A', 'C', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (undef, 2, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount, $MinKeyCount);
}

# Generate key 119 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 119 description: N=A
#
sub _Generate166KeySetKey119 {
  my($This) = @_;
  my($BondOrder) = 2;

  return $This->_DetectBondKeys('N', 'A', $BondOrder);
}

# Generate key 120 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 120 description: HETEROCYCLIC ATOM > 1 (&...)
#
sub _Generate166KeySetKey120 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing) = (2, 1);

  return $This->_DetectAtomKeys('Q', $MinKeyCount, $IsInRing);
}

# Generate key 121 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 121 description: N HETEROCYCLE
#
sub _Generate166KeySetKey121 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing) = (undef, 1);

  return $This->_DetectAtomKeys('N', $MinKeyCount, $IsInRing);
}

# Generate key 122 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 122 description: AN(A)A
#
sub _Generate166KeySetKey122 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'N';
  @NbrAtomSymbols = ('A', 'A', 'A');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 123 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 123 description: OCO
#
sub _Generate166KeySetKey123 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('O', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 124 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 124 description: QQ
#
sub _Generate166KeySetKey124 {
  my($This) = @_;

  return $This->_DetectBondKeys('Q', 'Q');
}

# Generate key 125 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 125 description: AROMATIC RING > 1
#
sub _Generate166KeySetKey125 {
  my($This) = @_;
  my($Molecule, $NumOfAromaticRings, $KeyValue);

  $Molecule = $This->GetMolecule();
  $NumOfAromaticRings = $Molecule->GetNumOfAromaticRings();

  if ($This->{KeyBits}) {
    $KeyValue = ($NumOfAromaticRings > 1) ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfAromaticRings;
  }
  return $KeyValue;
}

# Generate key 126 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 126 description: A!O!A
#
sub _Generate166KeySetKey126 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'O';
  @NbrAtomSymbols = ('A', 'A');
  @NbrBondSymbols = ('!', '!');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 127 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 127 description: A$A!O > 1 (&...)
#
sub _Generate166KeySetKey127 {
  my($This) = @_;
  my($CentralAtomSymbol, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'O');
  @NbrBondSymbols = ('$', '!');
  $MinKeyCount = 2;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount);
}

# Generate key 128 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 128 description: ACH2AAACH2A
#
sub _Generate166KeySetKey128 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('A', 'C', 'A', 'A', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (undef, 2, undef, undef, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 129 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 129 description: ACH2AACH2A
#
sub _Generate166KeySetKey129 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('A', 'C', 'A', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (undef, 2, undef, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 130 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 130 description: QQ > 1 (&...)
#
sub _Generate166KeySetKey130 {
  my($This) = @_;
  my($BondOrder, $MinKeyCount) = (undef, 2);

  return $This->_DetectBondKeys('Q', 'Q', $BondOrder, $MinKeyCount);
}

# Generate key 131 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 131 description: QH > 1
#
sub _Generate166KeySetKey131 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing, $MinHydrogenCount) = (2, undef, 1);

  return $This->_DetectAtomKeys('Q', $MinKeyCount, $IsInRing, $MinHydrogenCount);
}

# Generate key 132 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 132 description: OACH2A
#
sub _Generate166KeySetKey132 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('O', 'A', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (undef, undef, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 133 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 133 description: A$A!N
#
sub _Generate166KeySetKey133 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'N');
  @NbrBondSymbols = ('$', '!');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 134 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 134 description: X (HALOGEN)
#
sub _Generate166KeySetKey134 {
  my($This) = @_;

  return $This->_DetectAtomKeys('X');
}

# Generate key 135 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 135 description: Nnot%A%A
#
sub _Generate166KeySetKey135 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('N', 'A');
  @NbrBondSymbols = ('not%', '%');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 136 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 136 description: O=A > 1
#
sub _Generate166KeySetKey136 {
  my($This) = @_;
  my($BondOrder, $MinKeyCount) = (2, 2);

  return $This->_DetectBondKeys('O', 'A', $BondOrder, $MinKeyCount);
}

# Generate key 137 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 137 description: HETEROCYCLE
#
sub _Generate166KeySetKey137 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing) = (1, 1);

  return $This->_DetectAtomKeys('Q', $MinKeyCount, $IsInRing);
}

# Generate key 138 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 138 description: QCH2A > 1 (&...)
#
sub _Generate166KeySetKey138 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('Q', 'A');
  @NbrBondSymbols = (undef, undef);
  $MinKeyCount = 2;
  $CentralAtomMinHydrogenCount = 2;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate key 139 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 139 description: OH
#
sub _Generate166KeySetKey139 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing, $MinHydrogenCount) = (undef, undef, 1);

  return $This->_DetectAtomKeys('O', $MinKeyCount, $IsInRing, $MinHydrogenCount);
}

# Generate key 140 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 140 description: O > 3 (&...)
#
sub _Generate166KeySetKey140 {
  my($This) = @_;
  my($MinKeyCount) = 4;

  return $This->_DetectAtomKeys('O', $MinKeyCount);
}

# Generate key 141 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 141 description: CH3 > 2 (&...)
#
sub _Generate166KeySetKey141 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing, $MinHydrogenCount) = (3, undef, 3);

  return $This->_DetectAtomKeys('C', $MinKeyCount, $IsInRing, $MinHydrogenCount);
}

# Generate key 142 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 142 description: N > 1
#
sub _Generate166KeySetKey142 {
  my($This) = @_;
  my($MinKeyCount) = 2;

  return $This->_DetectAtomKeys('N', $MinKeyCount);
}

# Generate key 143 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 143 description: A$A!O
#
sub _Generate166KeySetKey143 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'O');
  @NbrBondSymbols = ('$', '!');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 144 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 144 description: Anot%A%Anot%A
#
sub _Generate166KeySetKey144 {
  my($This) = @_;
  my($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, @NbrAtomsSymbols, @NbrAtomsBondSymbols);

  ($BondAtomSymbol1, $BondAtomSymbol2) = ('A', 'A');
  $BondSymbol = '%';

  @NbrAtomsSymbols = (['A'], ['A']);
  @NbrAtomsBondSymbols = (['not%'], ['not%']);
  return $This->_DetectBondNeighborhoodKeys($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, \@NbrAtomsSymbols, \@NbrAtomsBondSymbols);
}

# Generate key 145 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 145 description: 6M RING > 1
#
sub _Generate166KeySetKey145 {
  my($This) = @_;
  my($Molecule, $KeyValue, $RingSize, $NumOfRings);

  $RingSize = 6;
  $Molecule = $This->GetMolecule();
  $NumOfRings = $Molecule->GetNumOfRingsWithSize($RingSize);

  if ($This->{KeyBits}) {
    $KeyValue = ($NumOfRings > 1) ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfRings;
  }
  return $KeyValue;
}

# Generate key 146 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 146 description: O > 2
#
sub _Generate166KeySetKey146 {
  my($This) = @_;
  my($MinKeyCount) = 3;

  return $This->_DetectAtomKeys('O', $MinKeyCount);
}

# Generate key 147 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 147 description: ACH2CH2A
#
sub _Generate166KeySetKey147 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount);

  @CentralAtomsSymbols = ('A', 'C', 'C', 'A');
  @CentralAtomsBondSymbols = (undef, undef, undef);
  @CentralAtomsMinHydrogenCount = (undef, 2, 2, undef);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount);
}

# Generate key 148 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 148 description: AQ(A)A
#
sub _Generate166KeySetKey148 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('A', 'A', 'A');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 149 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 149 description: CH3 > 1
#
sub _Generate166KeySetKey149 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing, $MinHydrogenCount) = (2, undef, 3);

  return $This->_DetectAtomKeys('C', $MinKeyCount, $IsInRing, $MinHydrogenCount);
}

# Generate key 150 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 150 description: A!A$A!A
#
sub _Generate166KeySetKey150 {
  my($This) = @_;
  my(@CentralAtomsSymbols, @CentralAtomsBondSymbols);

  @CentralAtomsSymbols = ('A', 'A', 'A', 'A');
  @CentralAtomsBondSymbols = ('!', '$', '!');

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols);
}

# Generate key 151 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 151 description: NH
#
sub _Generate166KeySetKey151 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing, $MinHydrogenCount) = (undef, undef, 1);

  return $This->_DetectAtomKeys('N', $MinKeyCount, $IsInRing, $MinHydrogenCount);
}

# Generate key 152 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 152 description: OC(C)C
#
sub _Generate166KeySetKey152 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('O', 'C', 'C');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate key 153 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 153 description: QCH2A
#
sub _Generate166KeySetKey153 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('Q', 'A');
  @NbrBondSymbols = (undef, undef);
  $MinKeyCount = undef;
  $CentralAtomMinHydrogenCount = 2;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate key 154 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 154 description: C=O
#
sub _Generate166KeySetKey154 {
  my($This) = @_;
  my($BondOrder) = 2;

  return $This->_DetectBondKeys('C', 'O', $BondOrder);
}

# Generate key 155 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 155 description: A!CH2!A
#
sub _Generate166KeySetKey155 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('A', 'A');
  @NbrBondSymbols = ('!', '!');
  $MinKeyCount = undef;
  $CentralAtomMinHydrogenCount = 2;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate key 156 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 156 description: NA(A)A
#
sub _Generate166KeySetKey156 {
  my($This) = @_;
  my($MinKeyCount, @CentralAtomsSymbols, @CentralAtomsBondSymbols, @CentralAtomsMinHydrogenCount, @CentralAtomNbrsAtomSymbols, @CentralAtomNbrsBondSymbols);

  @CentralAtomsSymbols = ('N', 'A', 'A');
  @CentralAtomsBondSymbols = (undef, undef);
  @CentralAtomsMinHydrogenCount = (undef, undef, undef);

  @CentralAtomNbrsAtomSymbols = (undef, ['A'], undef);
  @CentralAtomNbrsBondSymbols = (undef, undef, undef);
  $MinKeyCount = undef;

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, \@CentralAtomsBondSymbols, \@CentralAtomsMinHydrogenCount, $MinKeyCount, \@CentralAtomNbrsAtomSymbols, \@CentralAtomNbrsBondSymbols);
}

# Generate key 157 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 157 description: C-O
#
sub _Generate166KeySetKey157 {
  my($This) = @_;
  my($BondOrder) = 1;

  return $This->_DetectBondKeys('C', 'O', $BondOrder);
}

# Generate key 158 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 158 description: C-N
#
sub _Generate166KeySetKey158 {
  my($This) = @_;
  my($BondOrder) = 1;

  return $This->_DetectBondKeys('C', 'N', $BondOrder);
}

# Generate key 159 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 159 description: O > 1
#
sub _Generate166KeySetKey159 {
  my($This) = @_;
  my($MinKeyCount) = 2;

  return $This->_DetectAtomKeys('O', $MinKeyCount);
}

# Generate key 160 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 160 description: CH3
#
sub _Generate166KeySetKey160 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing, $MinHydrogenCount) = (undef, undef, 3);

  return $This->_DetectAtomKeys('C', $MinKeyCount, $IsInRing, $MinHydrogenCount);
}

# Generate key 161 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 161 description: N
#
sub _Generate166KeySetKey161 {
  my($This) = @_;

  return $This->_DetectAtomKeys('N');
}

# Generate key 162 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 162 description: AROMATIC
#
sub _Generate166KeySetKey162 {
  my($This) = @_;
  my($Atom, $Molecule, $KeyValue);

  # Check molecule aromatic property...
  $Molecule = $This->GetMolecule();
  if ($Molecule->IsAromatic()) {
    return 1;
  }

  # Check aromatic property of each atom...
  $KeyValue = 1;
  ATOM: for $Atom (@{$This->{Atoms}}) {
      if (!$Atom->IsAromatic()) {
	$KeyValue = 0;
	last ATOM;
      }
  }
  return $KeyValue;
}

# Generate key 163 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 163 description: 6M RING
#
sub _Generate166KeySetKey163 {
  my($This) = @_;
  my($Molecule, $KeyValue, $RingSize, $NumOfRings);

  $RingSize = 6;
  $Molecule = $This->GetMolecule();
  $NumOfRings = $Molecule->GetNumOfRingsWithSize($RingSize);

  if ($This->{KeyBits}) {
    $KeyValue = $NumOfRings ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfRings;
  }
  return $KeyValue;
}

# Generate key 164 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 164 description: O
#
sub _Generate166KeySetKey164 {
  my($This) = @_;

  return $This->_DetectAtomKeys('O');
}

# Generate key 165 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 165 description: RING
#
sub _Generate166KeySetKey165 {
  my($This) = @_;
  my($Molecule, $KeyValue, $NumOfRings);

  $Molecule = $This->GetMolecule();
  $NumOfRings = $Molecule->GetNumOfRings();

  if ($This->{KeyBits}) {
    $KeyValue = $NumOfRings ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfRings;
  }
  return $KeyValue;
}

# Generate key 166 value as 1/0 indicating its presence/absence or count of its
# presence in a molecule.
#
# Key 166 description: FRAGMENTS
#
sub _Generate166KeySetKey166 {
  my($This) = @_;
  my($Molecule, $KeyValue, $NumOfComponents);

  $Molecule = $This->GetMolecule();
  $NumOfComponents = $Molecule->GetNumOfConnectedComponents();

  if ($This->{KeyBits}) {
    $KeyValue = ($NumOfComponents > 1) ? 1 : 0;
  }
  else {
    $KeyValue = $NumOfComponents;
  }
  return $KeyValue;
}

##################################
#
#  Implementation of MDL MACCS 322 keys...
#
##################################

# Generate 322 keyset key 1 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 1 description: A(AAA) or AA(A)A - atom with at least three neighbors
#
sub _Generate322KeySetKey1 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'A', 'A');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 2 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 2 description: Q - heteroatom
#
sub _Generate322KeySetKey2 {
  my($This) = @_;

  return $This->_DetectAtomKeys('Q');
}

# Generate 322 keyset key 3 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 3 description: Anot%not-A - atom involved in one or more multiple bonds, not aromatic
#
sub _Generate322KeySetKey3 {
  my($This) = @_;
  my($BondSymbol) = 'not%not-';

  return $This->_DetectBondKeys('A', 'A', $BondSymbol);
}

# Generate 322 keyset key 4 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 4 description:      A(AAAA) or AA(A)(A)A - atom with at least four neighbors
#
sub _Generate322KeySetKey4 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'A', 'A', 'A');
  @NbrBondSymbols = (undef, undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 5 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 5 description: A(QQ) or QA(Q) - atom with at least two heteroatom neighbors
#
sub _Generate322KeySetKey5 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('Q', 'Q');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 6 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 6 description: A(QQQ) or QA(Q)Q - atom with at least three heteroatom neighbors
#
sub _Generate322KeySetKey6 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('Q', 'Q', 'Q');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 7 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 7 description:      QH - heteroatom with at least one hydrogen attached
#
sub _Generate322KeySetKey7 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing, $MinHydrogenCount) = (undef, undef, 1);

  return $This->_DetectAtomKeys('Q', $MinKeyCount, $IsInRing, $MinHydrogenCount);
}

# Generate 322 keyset key 8 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 8 description: CH2(AA) or ACH2A - carbon with at least two single bonds and at least two hydrogens attached
#
sub _Generate322KeySetKey8 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('A', 'A');
  @NbrBondSymbols = (undef, undef);
  $MinKeyCount = undef;
  $CentralAtomMinHydrogenCount = 2;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate 322 keyset key 9 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 9 description: CH3(A) or ACH3 - carbon with at least one single bond and at least three hydrogens attached
#
sub _Generate322KeySetKey9 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('A');
  @NbrBondSymbols = (undef);
  $MinKeyCount = undef;
  $CentralAtomMinHydrogenCount = 3;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount);
}

# Generate 322 keyset key 10 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 10 description: Halogen
#
sub _Generate322KeySetKey10 {
  my($This) = @_;

  return $This->_DetectAtomKeys('X');
}

# Generate 322 keyset key 11 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 11 description: A(-A-A-A) or A-A(-A)-A - atom has at least three single bonds
#
sub _Generate322KeySetKey11 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'A', 'A');
  @NbrBondSymbols = ('-', '-', '-');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 12 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 12 description: AAAAAA@1 >= 2 - atom is in at least two different six-membered rings
#
sub _Generate322KeySetKey12 {
  my($This) = @_;
  my($Atom, $KeyValue, $RingSize, $NumOfRings);

  $RingSize = 6;
  $KeyValue = 0;

  ATOM: for $Atom (@{$This->{Atoms}}) {
    if (!$This->_IsAtom($Atom)) {
      next ATOM;
    }
    $NumOfRings = $Atom->GetNumOfRingsWithSize($RingSize);
    if ($NumOfRings >= 2) {
      $KeyValue++;
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
    }
  }
  return $KeyValue;
}

# Generate 322 keyset key 13 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 13 description: A($A$A$A) or A$A($A)$A - atom has more than two ring bonds (at least three ring bonds)
#
sub _Generate322KeySetKey13 {
  my($This) = @_;
  my($CentralAtomSymbol, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'A', 'A');
  @NbrBondSymbols = ('$', '$', '$');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 14 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 14 description: A$A!A$A - atom is at a ring/chain boundary. When a comparison is
#                     done with another atom the path passes through the chain bond.
#
sub _Generate322KeySetKey14 {
  my($This) = @_;
  my($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, @NbrAtomsSymbols, @NbrAtomsBondSymbols);

  ($BondAtomSymbol1, $BondAtomSymbol2) = ('A', 'A');
  $BondSymbol = '!';

  @NbrAtomsSymbols = (['A'], ['A']);
  @NbrAtomsBondSymbols = (['$'], ['$']);
  return $This->_DetectBondNeighborhoodKeys($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, \@NbrAtomsSymbols, \@NbrAtomsBondSymbols);
}

# Generate 322 keyset key 15 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 15 description:  Anot%A%Anot%A - atom is at an aromatic/nonaromatic boundary.
#                      When a comparison is done with another atom the path passes through the aromatic bond.
#
sub _Generate322KeySetKey15 {
  my($This) = @_;
  my($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, @NbrAtomsSymbols, @NbrAtomsBondSymbols);

  ($BondAtomSymbol1, $BondAtomSymbol2) = ('A', 'A');
  $BondSymbol = '%';

  @NbrAtomsSymbols = (['A'], ['A']);
  @NbrAtomsBondSymbols = (['not%'], ['not%']);
  return $This->_DetectBondNeighborhoodKeys($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, \@NbrAtomsSymbols, \@NbrAtomsBondSymbols);
}

# Generate 322 keyset key 16 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 16 description:     A!A!A  - atom with more than one chain bond
#
sub _Generate322KeySetKey16 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'A';
  @NbrAtomSymbols = ('A', 'A');
  @NbrBondSymbols = ('!', '!');

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 17 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 17 description: A!A$A!A - atom is at a ring/chain boundary. When a comparison
#                     is done  with another atom the path passes through the ring bond.
#
sub _Generate322KeySetKey17 {
  my($This) = @_;
  my($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, @NbrAtomsSymbols, @NbrAtomsBondSymbols);

  ($BondAtomSymbol1, $BondAtomSymbol2) = ('A', 'A');
  $BondSymbol = '$';

  @NbrAtomsSymbols = (['A'], ['A']);
  @NbrAtomsBondSymbols = (['!'], ['!']);
  return $This->_DetectBondNeighborhoodKeys($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, \@NbrAtomsSymbols, \@NbrAtomsBondSymbols);
}

# Generate 322 keyset key 18 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 18 description: A%Anot%A%A - atom is at an aromatic/nonaromatic boundary.
#                     When a comparison is done with another atom the path passes through
#                     the nonaromatic bond
#
sub _Generate322KeySetKey18 {
  my($This) = @_;
  my($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, @NbrAtomsSymbols, @NbrAtomsBondSymbols);

  ($BondAtomSymbol1, $BondAtomSymbol2) = ('A', 'A');
  $BondSymbol = 'not%';

  @NbrAtomsSymbols = (['A'], ['A']);
  @NbrAtomsBondSymbols = (['%'], ['%']);
  return $This->_DetectBondNeighborhoodKeys($BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, \@NbrAtomsSymbols, \@NbrAtomsBondSymbols);
}

# Generate 322 keyset key 19 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 19 description: HETEROCYCLE - atom is a heteroatom in a ring.
#
sub _Generate322KeySetKey19 {
  my($This) = @_;
  my($MinKeyCount, $IsInRing) = (undef, 1);

  return $This->_DetectAtomKeys('Q', $MinKeyCount, $IsInRing);
}

# Generate 322 keyset key 20 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 20 description: rare properties: atom with five or more neighbors, atom in four
#                     or more rings, or atom types other than  H, C, N, O, S, F, Cl, Br, or I
#
sub _Generate322KeySetKey20 {
  my($This) = @_;
  my($Atom, $KeyValue);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if (!($Atom->GetAtomicNumber() !~ /^(1|6|7|8|9|16|17|35|53)$/) || ($Atom->GetNumOfRings() >= 4) || ($Atom->GetNumOfNeighbors() >= 5) ) {
      next ATOM;
    }
    $KeyValue++;
    if ($This->{KeyBits}) {
      $KeyValue = 1;
      last ATOM;
    }
  }
  return $KeyValue;
}

# Generate 322 keyset key 21 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 21 description: rare properties: atom has a charge, is an isotope, has
#                     two or more multiple bonds, or has a triple bond.
#
sub _Generate322KeySetKey21 {
  my($This) = @_;
  my($Atom, $KeyValue);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if ( !($Atom->IsIsotope() || $Atom->GetFormalCharge()) ) {
      # Look for multiple and triple bonds...
      my($Bond, $NumOfTripleBonds, $NumOfMultipleBonds);

      ($NumOfTripleBonds, $NumOfMultipleBonds) = (0, 0);
      BOND: for $Bond ($Atom->GetBonds()) {
	if ($Bond->IsSingle()) { next BOND; }
	if ($Bond->IsDouble()) { $NumOfMultipleBonds++; next BOND; }
	if ($Bond->IsTriple()) { $NumOfTripleBonds++; next BOND; }
      }
      if ( !($NumOfTripleBonds || ($NumOfMultipleBonds >= 2)) ) {
	next ATOM;
      }
    }
    $KeyValue++;
    if ($This->{KeyBits}) {
      $KeyValue = 1;
      last ATOM;
    }
  }
  return $KeyValue;
}

# Generate 322 keyset key 22 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 22 description:  N - nitrogen
#
sub _Generate322KeySetKey22 {
  my($This) = @_;

  return $This->_DetectAtomKeys('N');
}

# Generate 322 keyset key 23 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 23 description: S - sulfur
#
sub _Generate322KeySetKey23 {
  my($This) = @_;

  return $This->_DetectAtomKeys('S');
}

# Generate 322 keyset key 24 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 24 description: O - oxygen
#
sub _Generate322KeySetKey24 {
  my($This) = @_;

  return $This->_DetectAtomKeys('O');
}

# Generate 322 keyset key 25 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 25 description: A(AA)A(A)A(AA) - atom has two neighbors, each with
#                     three or more neighbors (including the central atom).
#
sub _Generate322KeySetKey25 {
  my($This) = @_;
  my($MinKeyCount, @CentralAtomsSymbols, @NbrAtomsSymbols);

  @CentralAtomsSymbols = ('A', 'A', 'A');
  @NbrAtomsSymbols = (['A', 'A'], ['A'], ['A', 'A']);

  return $This->_DetectExtendedAtomNeighborhoodKeys(\@CentralAtomsSymbols, undef, undef, undef, \@NbrAtomsSymbols);
}

# Generate 322 keyset key 26 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 26 description:     CH2ACH2 - atom has two hydrocarbon (CH2) neighbors
#
sub _Generate322KeySetKey26 {
  my($This) = @_;
  my($CentralAtomSymbol, $CentralAtomMinHydrogenCount, $MinKeyCount, @NbrAtomSymbols, @NbrBondSymbols, @NbrAtomMinHydrogenCount);

  $CentralAtomSymbol = 'A';
  $CentralAtomMinHydrogenCount = undef;

  @NbrAtomSymbols = ('C', 'C');
  @NbrBondSymbols = (undef, undef);
  @NbrAtomMinHydrogenCount = (2, 2);

  $MinKeyCount = undef;

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinKeyCount, $CentralAtomMinHydrogenCount, \@NbrAtomMinHydrogenCount);
}

# Generate 322 keyset key 27 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 27 description: C(CC)
#
sub _Generate322KeySetKey27 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'C');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 28 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 28 description: C(CCC)
#
sub _Generate322KeySetKey28 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'C', 'C');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 29 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 29 description: C(CN)
#
sub _Generate322KeySetKey29 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'N');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 30 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 30 description: C(CCN)
#
sub _Generate322KeySetKey30 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'C', 'N');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 31 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 31 description: C(NN)
#
sub _Generate322KeySetKey31 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'N');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 32 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 32 description: C(NNC)
#
sub _Generate322KeySetKey32 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'N', 'C');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 33 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 33 description: C(NNN)
#
sub _Generate322KeySetKey33 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'N', 'N');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 34 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 34 description: C(CO)
#
sub _Generate322KeySetKey34 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 35 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 35 description: C(CCO)
#
sub _Generate322KeySetKey35 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'C', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 36 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 36 description: C(NO)
#
sub _Generate322KeySetKey36 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 37 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 37 description: C(NCO)
#
sub _Generate322KeySetKey37 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'C', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 38 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 38 description: C(NNO)
#
sub _Generate322KeySetKey38 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'N', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 39 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 39 description: C(OO)
#
sub _Generate322KeySetKey39 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('O', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 40 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 40 description: C(COO)
#
sub _Generate322KeySetKey40 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('C', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 41 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 41 description: C(NOO)
#
sub _Generate322KeySetKey41 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('N', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 42 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 42 description: C(OOO)
#
sub _Generate322KeySetKey42 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'C';
  @NbrAtomSymbols = ('O', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 43 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 43 description: Q(CC)
#
sub _Generate322KeySetKey43 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'C');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 44 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 44 description: Q(CCC)
#
sub _Generate322KeySetKey44 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'C', 'C');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 45 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 45 description: Q(CN)
#
sub _Generate322KeySetKey45 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'N');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 46 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 46 description: Q(CCN)
#
sub _Generate322KeySetKey46 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'C', 'N');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 47 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 47 description: Q(NN)
#
sub _Generate322KeySetKey47 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('N', 'N');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 48 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 48 description: Q(CNN)
#
sub _Generate322KeySetKey48 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'N', 'N');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 49 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 49 description: Q(NNN)
#
sub _Generate322KeySetKey49 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('N', 'N', 'N');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 50 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 50 description: Q(CO)
#
sub _Generate322KeySetKey50 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 51 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 51 description: Q(CCO)
#
sub _Generate322KeySetKey51 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'C', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 52 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 52 description: Q(NO)
#
sub _Generate322KeySetKey52 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('N', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 53 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 53 description: Q(CNO)
#
sub _Generate322KeySetKey53 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'N', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 54 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 54 description: Q(NNO)
#
sub _Generate322KeySetKey54 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('N', 'N', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 55 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 55 description: Q(OO)
#
sub _Generate322KeySetKey55 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('O', 'O');
  @NbrBondSymbols = (undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 56 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 56 description: Q(COO)
#
sub _Generate322KeySetKey56 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('C', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 57 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 57 description: Q(NOO)
#
sub _Generate322KeySetKey57 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('N', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 58 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 58 description: Q(OOO)
#
sub _Generate322KeySetKey58 {
  my($This) = @_;
  my($CentralAtomSymbol, @NbrAtomSymbols, @NbrBondSymbols);

  $CentralAtomSymbol = 'Q';
  @NbrAtomSymbols = ('O', 'O', 'O');
  @NbrBondSymbols = (undef, undef, undef);

  return $This->_DetectAtomNeighborhoodKeys($CentralAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols);
}

# Generate 322 keyset key 59 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 59 description: C-C
#
sub _Generate322KeySetKey59 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'C', $BondSymbol);
}

# Generate 322 keyset key 60 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 60 description: C-N
#
sub _Generate322KeySetKey60 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'N', $BondSymbol);
}

# Generate 322 keyset key 61 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 61 description: C-O
#
sub _Generate322KeySetKey61 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'O', $BondSymbol);
}

# Generate 322 keyset key 62 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 62 description: C-S
#
sub _Generate322KeySetKey62 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'S', $BondSymbol);
}

# Generate 322 keyset key 63 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 63 description: C-Cl
#
sub _Generate322KeySetKey63 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 64 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 64 description: C-P
#
sub _Generate322KeySetKey64 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'P', $BondSymbol);
}

# Generate 322 keyset key 65 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 65 description: C-F
#
sub _Generate322KeySetKey65 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'F', $BondSymbol);
}

# Generate 322 keyset key 66 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 66 description: C-Br
#
sub _Generate322KeySetKey66 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'Br', $BondSymbol);
}

# Generate 322 keyset key 67 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 67 description: C-Si
#
sub _Generate322KeySetKey67 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'Si', $BondSymbol);
}

# Generate 322 keyset key 68 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 68 description: C-I
#
sub _Generate322KeySetKey68 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'I', $BondSymbol);
}

# Generate 322 keyset key 69 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 69 description: C-X
#
sub _Generate322KeySetKey69 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('C', 'Z', $BondSymbol);
}

# Generate 322 keyset key 70 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 70 description: N-N
#
sub _Generate322KeySetKey70 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'N', $BondSymbol);
}

# Generate 322 keyset key 71 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 71 description: N-O
#
sub _Generate322KeySetKey71 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'O', $BondSymbol);
}

# Generate 322 keyset key 72 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 72 description: N-S
#
sub _Generate322KeySetKey72 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'S', $BondSymbol);
}

# Generate 322 keyset key 73 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 73 description: N-Cl
#
sub _Generate322KeySetKey73 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 74 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 74 description: N-P
#
sub _Generate322KeySetKey74 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'P', $BondSymbol);
}

# Generate 322 keyset key 75 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 75 description: N-F
#
sub _Generate322KeySetKey75 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'F', $BondSymbol);
}

# Generate 322 keyset key 76 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 76 description: N-Br
#
sub _Generate322KeySetKey76 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'Br', $BondSymbol);
}

# Generate 322 keyset key 77 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 77 description: N-Si
#
sub _Generate322KeySetKey77 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'Si', $BondSymbol);
}

# Generate 322 keyset key 78 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 78 description: N-I
#
sub _Generate322KeySetKey78 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'I', $BondSymbol);
}

# Generate 322 keyset key 79 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 79 description: N-X
#
sub _Generate322KeySetKey79 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('N', 'Z', $BondSymbol);
}

# Generate 322 keyset key 80 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 80 description: O-O
#
sub _Generate322KeySetKey80 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'O', $BondSymbol);
}

# Generate 322 keyset key 81 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 81 description: O-S
#
sub _Generate322KeySetKey81 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'S', $BondSymbol);
}

# Generate 322 keyset key 82 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 82 description: O-Cl
#
sub _Generate322KeySetKey82 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 83 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 83 description: O-P
#
sub _Generate322KeySetKey83 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'P', $BondSymbol);
}

# Generate 322 keyset key 84 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 84 description: O-F
#
sub _Generate322KeySetKey84 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'F', $BondSymbol);
}

# Generate 322 keyset key 85 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 85 description: O-Br
#
sub _Generate322KeySetKey85 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'Br', $BondSymbol);
}

# Generate 322 keyset key 86 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 86 description: O-Si
#
sub _Generate322KeySetKey86 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'Si', $BondSymbol);
}

# Generate 322 keyset key 87 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 87 description: O-I
#
sub _Generate322KeySetKey87 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'I', $BondSymbol);
}

# Generate 322 keyset key 88 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 88 description: O-X
#
sub _Generate322KeySetKey88 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('O', 'Z', $BondSymbol);
}

# Generate 322 keyset key 89 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 89 description: S-S
#
sub _Generate322KeySetKey89 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('S', 'S', $BondSymbol);
}

# Generate 322 keyset key 90 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 90 description: S-Cl
#
sub _Generate322KeySetKey90 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('S', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 91 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 91 description: S-P
#
sub _Generate322KeySetKey91 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('S', 'P', $BondSymbol);
}

# Generate 322 keyset key 92 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 92 description: S-F
#
sub _Generate322KeySetKey92 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('S', 'F', $BondSymbol);
}

# Generate 322 keyset key 93 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 93 description: S-Br
#
sub _Generate322KeySetKey93 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('S', 'Br', $BondSymbol);
}

# Generate 322 keyset key 94 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 94 description: S-Si
#
sub _Generate322KeySetKey94 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('S', 'Si', $BondSymbol);
}

# Generate 322 keyset key 95 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 95 description: S-I
#
sub _Generate322KeySetKey95 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('S', 'I', $BondSymbol);
}

# Generate 322 keyset key 96 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 96 description: S-X
#
sub _Generate322KeySetKey96 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('S', 'Z', $BondSymbol);
}

# Generate 322 keyset key 97 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 97 description: Cl-Cl
#
sub _Generate322KeySetKey97 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Cl', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 98 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 98 description: Cl-P
#
sub _Generate322KeySetKey98 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Cl', 'P', $BondSymbol);
}

# Generate 322 keyset key 99 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 99 description: Cl-F
#
sub _Generate322KeySetKey99 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Cl', 'F', $BondSymbol);
}

# Generate 322 keyset key 100 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 100 description: Cl-Br
#
sub _Generate322KeySetKey100 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Cl', 'Br', $BondSymbol);
}

# Generate 322 keyset key 101 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 101 description: Cl-Si
#
sub _Generate322KeySetKey101 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Cl', 'Si', $BondSymbol);
}

# Generate 322 keyset key 102 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 102 description: Cl-I
#
sub _Generate322KeySetKey102 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Cl', 'I', $BondSymbol);
}

# Generate 322 keyset key 103 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 103 description: Cl-X
#
sub _Generate322KeySetKey103 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Cl', 'Z', $BondSymbol);
}

# Generate 322 keyset key 104 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 104 description: P-P
#
sub _Generate322KeySetKey104 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('P', 'P', $BondSymbol);
}

# Generate 322 keyset key 105 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 105 description: P-F
#
sub _Generate322KeySetKey105 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('P', 'F', $BondSymbol);
}

# Generate 322 keyset key 106 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 106 description: P-Br
#
sub _Generate322KeySetKey106 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('P', 'Br', $BondSymbol);
}

# Generate 322 keyset key 107 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 107 description: P-Si
#
sub _Generate322KeySetKey107 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('P', 'Si', $BondSymbol);
}

# Generate 322 keyset key 108 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 108 description: P-I
#
sub _Generate322KeySetKey108 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('P', 'I', $BondSymbol);
}

# Generate 322 keyset key 109 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 109 description: P-X
#
sub _Generate322KeySetKey109 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('P', 'Z', $BondSymbol);
}

# Generate 322 keyset key 110 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 110 description: F-F
#
sub _Generate322KeySetKey110 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('F', 'F', $BondSymbol);
}

# Generate 322 keyset key 111 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 111 description: F-Br
#
sub _Generate322KeySetKey111 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('F', 'Br', $BondSymbol);
}

# Generate 322 keyset key 112 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 112 description: F-Si
#
sub _Generate322KeySetKey112 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('F', 'Si', $BondSymbol);
}

# Generate 322 keyset key 113 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 113 description: F-I
#
sub _Generate322KeySetKey113 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('F', 'I', $BondSymbol);
}

# Generate 322 keyset key 114 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 114 description: F-X
#
sub _Generate322KeySetKey114 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('F', 'Z', $BondSymbol);
}

# Generate 322 keyset key 115 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 115 description: Br-Br
#
sub _Generate322KeySetKey115 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Br', 'Br', $BondSymbol);
}

# Generate 322 keyset key 116 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 116 description: Br-Si
#
sub _Generate322KeySetKey116 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Br', 'Si', $BondSymbol);
}

# Generate 322 keyset key 117 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 117 description: Br-I
#
sub _Generate322KeySetKey117 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Br', 'I', $BondSymbol);
}

# Generate 322 keyset key 118 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 118 description: Br-X
#
sub _Generate322KeySetKey118 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Br', 'Z', $BondSymbol);
}

# Generate 322 keyset key 119 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 119 description: Si-Si
#
sub _Generate322KeySetKey119 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Si', 'Si', $BondSymbol);
}

# Generate 322 keyset key 120 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 120 description: Si-I
#
sub _Generate322KeySetKey120 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Si', 'I', $BondSymbol);
}

# Generate 322 keyset key 121 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 121 description: Si-X
#
sub _Generate322KeySetKey121 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Si', 'Z', $BondSymbol);
}

# Generate 322 keyset key 122 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 122 description: I-I
#
sub _Generate322KeySetKey122 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('I', 'I', $BondSymbol);
}

# Generate 322 keyset key 123 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 123 description: I-X
#
sub _Generate322KeySetKey123 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('I', 'Z', $BondSymbol);
}

# Generate 322 keyset key 124 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 124 description: X-X
#
sub _Generate322KeySetKey124 {
  my($This) = @_;
  my($BondSymbol) = '-';

  return $This->_DetectBondKeys('Z', 'Z', $BondSymbol);
}

# Generate 322 keyset key 125 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 125 description: C=C
#
sub _Generate322KeySetKey125 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'C', $BondSymbol);
}

# Generate 322 keyset key 126 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 126 description: C=N
#
sub _Generate322KeySetKey126 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'N', $BondSymbol);
}

# Generate 322 keyset key 127 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 127 description: C=O
#
sub _Generate322KeySetKey127 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'O', $BondSymbol);
}

# Generate 322 keyset key 128 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 128 description: C=S
#
sub _Generate322KeySetKey128 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'S', $BondSymbol);
}

# Generate 322 keyset key 129 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 129 description: C=Cl
#
sub _Generate322KeySetKey129 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 130 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 130 description: C=P
#
sub _Generate322KeySetKey130 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'P', $BondSymbol);
}

# Generate 322 keyset key 131 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 131 description: C=F
#
sub _Generate322KeySetKey131 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'F', $BondSymbol);
}

# Generate 322 keyset key 132 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 132 description: C=Br
#
sub _Generate322KeySetKey132 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'Br', $BondSymbol);
}

# Generate 322 keyset key 133 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 133 description: C=Si
#
sub _Generate322KeySetKey133 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'Si', $BondSymbol);
}

# Generate 322 keyset key 134 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 134 description: C=I
#
sub _Generate322KeySetKey134 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'I', $BondSymbol);
}

# Generate 322 keyset key 135 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 135 description: C=X
#
sub _Generate322KeySetKey135 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('C', 'Z', $BondSymbol);
}

# Generate 322 keyset key 136 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 136 description: N=N
#
sub _Generate322KeySetKey136 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'N', $BondSymbol);
}

# Generate 322 keyset key 137 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 137 description: N=O
#
sub _Generate322KeySetKey137 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'O', $BondSymbol);
}

# Generate 322 keyset key 138 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 138 description: N=S
#
sub _Generate322KeySetKey138 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'S', $BondSymbol);
}

# Generate 322 keyset key 139 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 139 description: N=Cl
#
sub _Generate322KeySetKey139 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 140 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 140 description: N=P
#
sub _Generate322KeySetKey140 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'P', $BondSymbol);
}

# Generate 322 keyset key 141 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 141 description: N=F
#
sub _Generate322KeySetKey141 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'F', $BondSymbol);
}

# Generate 322 keyset key 142 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 142 description: N=Br
#
sub _Generate322KeySetKey142 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'Br', $BondSymbol);
}

# Generate 322 keyset key 143 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 143 description: N=Si
#
sub _Generate322KeySetKey143 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'Si', $BondSymbol);
}

# Generate 322 keyset key 144 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 144 description: N=I
#
sub _Generate322KeySetKey144 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'I', $BondSymbol);
}

# Generate 322 keyset key 145 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 145 description: N=X
#
sub _Generate322KeySetKey145 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('N', 'Z', $BondSymbol);
}

# Generate 322 keyset key 146 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 146 description: O=O
#
sub _Generate322KeySetKey146 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'O', $BondSymbol);
}

# Generate 322 keyset key 147 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 147 description: O=S
#
sub _Generate322KeySetKey147 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'S', $BondSymbol);
}

# Generate 322 keyset key 148 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 148 description: O=Cl
#
sub _Generate322KeySetKey148 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 149 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 149 description: O=P
#
sub _Generate322KeySetKey149 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'P', $BondSymbol);
}

# Generate 322 keyset key 150 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 150 description: O=F
#
sub _Generate322KeySetKey150 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'F', $BondSymbol);
}

# Generate 322 keyset key 151 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 151 description: O=Br
#
sub _Generate322KeySetKey151 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'Br', $BondSymbol);
}

# Generate 322 keyset key 152 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 152 description: O=Si
#
sub _Generate322KeySetKey152 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'Si', $BondSymbol);
}

# Generate 322 keyset key 153 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 153 description: O=I
#
sub _Generate322KeySetKey153 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'I', $BondSymbol);
}

# Generate 322 keyset key 154 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 154 description: O=X
#
sub _Generate322KeySetKey154 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('O', 'Z', $BondSymbol);
}

# Generate 322 keyset key 155 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 155 description: S=S
#
sub _Generate322KeySetKey155 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('S', 'S', $BondSymbol);
}

# Generate 322 keyset key 156 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 156 description: S=Cl
#
sub _Generate322KeySetKey156 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('S', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 157 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 157 description: S=P
#
sub _Generate322KeySetKey157 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('S', 'P', $BondSymbol);
}

# Generate 322 keyset key 158 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 158 description: S=F
#
sub _Generate322KeySetKey158 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('S', 'F', $BondSymbol);
}

# Generate 322 keyset key 159 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 159 description: S=Br
#
sub _Generate322KeySetKey159 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('S', 'Br', $BondSymbol);
}

# Generate 322 keyset key 160 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 160 description: S=Si
#
sub _Generate322KeySetKey160 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('S', 'Si', $BondSymbol);
}

# Generate 322 keyset key 161 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 161 description: S=I
#
sub _Generate322KeySetKey161 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('S', 'I', $BondSymbol);
}

# Generate 322 keyset key 162 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 162 description: S=X
#
sub _Generate322KeySetKey162 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('S', 'Z', $BondSymbol);
}

# Generate 322 keyset key 163 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 163 description: Cl=Cl
#
sub _Generate322KeySetKey163 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Cl', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 164 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 164 description: Cl=P
#
sub _Generate322KeySetKey164 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Cl', 'P', $BondSymbol);
}

# Generate 322 keyset key 165 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 165 description: Cl=F
#
sub _Generate322KeySetKey165 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Cl', 'F', $BondSymbol);
}

# Generate 322 keyset key 166 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 166 description: Cl=Br
#
sub _Generate322KeySetKey166 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Cl', 'Br', $BondSymbol);
}

# Generate 322 keyset key 167 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 167 description: Cl=Si
#
sub _Generate322KeySetKey167 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Cl', 'Si', $BondSymbol);
}

# Generate 322 keyset key 168 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 168 description: Cl=I
#
sub _Generate322KeySetKey168 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Cl', 'I', $BondSymbol);
}

# Generate 322 keyset key 169 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 169 description: Cl=X
#
sub _Generate322KeySetKey169 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Cl', 'Z', $BondSymbol);
}

# Generate 322 keyset key 170 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 170 description: P=P
#
sub _Generate322KeySetKey170 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('P', 'P', $BondSymbol);
}

# Generate 322 keyset key 171 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 171 description: P=F
#
sub _Generate322KeySetKey171 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('P', 'F', $BondSymbol);
}

# Generate 322 keyset key 172 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 172 description: P=Br
#
sub _Generate322KeySetKey172 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('P', 'Br', $BondSymbol);
}

# Generate 322 keyset key 173 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 173 description: P=Si
#
sub _Generate322KeySetKey173 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('P', 'Si', $BondSymbol);
}

# Generate 322 keyset key 174 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 174 description: P=I
#
sub _Generate322KeySetKey174 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('P', 'I', $BondSymbol);
}

# Generate 322 keyset key 175 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 175 description: P=X
#
sub _Generate322KeySetKey175 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('P', 'Z', $BondSymbol);
}

# Generate 322 keyset key 176 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 176 description: F=F
#
sub _Generate322KeySetKey176 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('F', 'F', $BondSymbol);
}

# Generate 322 keyset key 177 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 177 description: F=Br
#
sub _Generate322KeySetKey177 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('F', 'Br', $BondSymbol);
}

# Generate 322 keyset key 178 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 178 description: F=Si
#
sub _Generate322KeySetKey178 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('F', 'Si', $BondSymbol);
}

# Generate 322 keyset key 179 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 179 description: F=I
#
sub _Generate322KeySetKey179 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('F', 'I', $BondSymbol);
}

# Generate 322 keyset key 180 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 180 description: F=X
#
sub _Generate322KeySetKey180 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('F', 'Z', $BondSymbol);
}

# Generate 322 keyset key 181 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 181 description: Br=Br
#
sub _Generate322KeySetKey181 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Br', 'Br', $BondSymbol);
}

# Generate 322 keyset key 182 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 182 description: Br=Si
#
sub _Generate322KeySetKey182 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Br', 'Si', $BondSymbol);
}

# Generate 322 keyset key 183 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 183 description: Br=I
#
sub _Generate322KeySetKey183 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Br', 'I', $BondSymbol);
}

# Generate 322 keyset key 184 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 184 description: Br=X
#
sub _Generate322KeySetKey184 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Br', 'Z', $BondSymbol);
}

# Generate 322 keyset key 185 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 185 description: Si=Si
#
sub _Generate322KeySetKey185 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Si', 'Si', $BondSymbol);
}

# Generate 322 keyset key 186 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 186 description: Si=I
#
sub _Generate322KeySetKey186 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Si', 'I', $BondSymbol);
}

# Generate 322 keyset key 187 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 187 description: Si=X
#
sub _Generate322KeySetKey187 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Si', 'Z', $BondSymbol);
}

# Generate 322 keyset key 188 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 188 description: I=I
#
sub _Generate322KeySetKey188 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('I', 'I', $BondSymbol);
}

# Generate 322 keyset key 189 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 189 description: I=X
#
sub _Generate322KeySetKey189 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('I', 'Z', $BondSymbol);
}

# Generate 322 keyset key 190 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 190 description: X=X
#
sub _Generate322KeySetKey190 {
  my($This) = @_;
  my($BondSymbol) = '=';

  return $This->_DetectBondKeys('Z', 'Z', $BondSymbol);
}

# Generate 322 keyset key 191 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 191 description: C#C
#
sub _Generate322KeySetKey191 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'C', $BondSymbol);
}

# Generate 322 keyset key 192 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 192 description: C#N
#
sub _Generate322KeySetKey192 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'N', $BondSymbol);
}

# Generate 322 keyset key 193 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 193 description: C#O
#
sub _Generate322KeySetKey193 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'O', $BondSymbol);
}

# Generate 322 keyset key 194 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 194 description: C#S
#
sub _Generate322KeySetKey194 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'S', $BondSymbol);
}

# Generate 322 keyset key 195 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 195 description: C#Cl
#
sub _Generate322KeySetKey195 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 196 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 196 description: C#P
#
sub _Generate322KeySetKey196 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'P', $BondSymbol);
}

# Generate 322 keyset key 197 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 197 description: C#F
#
sub _Generate322KeySetKey197 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'F', $BondSymbol);
}

# Generate 322 keyset key 198 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 198 description: C#Br
#
sub _Generate322KeySetKey198 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'Br', $BondSymbol);
}

# Generate 322 keyset key 199 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 199 description: C#Si
#
sub _Generate322KeySetKey199 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'Si', $BondSymbol);
}

# Generate 322 keyset key 200 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 200 description: C#I
#
sub _Generate322KeySetKey200 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'I', $BondSymbol);
}

# Generate 322 keyset key 201 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 201 description: C#X
#
sub _Generate322KeySetKey201 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('C', 'Z', $BondSymbol);
}

# Generate 322 keyset key 202 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 202 description: N#N
#
sub _Generate322KeySetKey202 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'N', $BondSymbol);
}

# Generate 322 keyset key 203 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 203 description: N#O
#
sub _Generate322KeySetKey203 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'O', $BondSymbol);
}

# Generate 322 keyset key 204 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 204 description: N#S
#
sub _Generate322KeySetKey204 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'S', $BondSymbol);
}

# Generate 322 keyset key 205 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 205 description: N#Cl
#
sub _Generate322KeySetKey205 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 206 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 206 description: N#P
#
sub _Generate322KeySetKey206 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'P', $BondSymbol);
}

# Generate 322 keyset key 207 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 207 description: N#F
#
sub _Generate322KeySetKey207 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'F', $BondSymbol);
}

# Generate 322 keyset key 208 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 208 description: N#Br
#
sub _Generate322KeySetKey208 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'Br', $BondSymbol);
}

# Generate 322 keyset key 209 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 209 description: N#Si
#
sub _Generate322KeySetKey209 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'Si', $BondSymbol);
}

# Generate 322 keyset key 210 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 210 description: N#I
#
sub _Generate322KeySetKey210 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'I', $BondSymbol);
}

# Generate 322 keyset key 211 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 211 description: N#X
#
sub _Generate322KeySetKey211 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('N', 'Z', $BondSymbol);
}

# Generate 322 keyset key 212 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 212 description: O#O
#
sub _Generate322KeySetKey212 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'O', $BondSymbol);
}

# Generate 322 keyset key 213 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 213 description: O#S
#
sub _Generate322KeySetKey213 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'S', $BondSymbol);
}

# Generate 322 keyset key 214 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 214 description: O#Cl
#
sub _Generate322KeySetKey214 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 215 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 215 description: O#P
#
sub _Generate322KeySetKey215 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'P', $BondSymbol);
}

# Generate 322 keyset key 216 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 216 description: O#F
#
sub _Generate322KeySetKey216 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'F', $BondSymbol);
}

# Generate 322 keyset key 217 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 217 description: O#Br
#
sub _Generate322KeySetKey217 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'Br', $BondSymbol);
}

# Generate 322 keyset key 218 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 218 description: O#Si
#
sub _Generate322KeySetKey218 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'Si', $BondSymbol);
}

# Generate 322 keyset key 219 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 219 description: O#I
#
sub _Generate322KeySetKey219 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'I', $BondSymbol);
}

# Generate 322 keyset key 220 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 220 description: O#X
#
sub _Generate322KeySetKey220 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('O', 'Z', $BondSymbol);
}

# Generate 322 keyset key 221 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 221 description: S#S
#
sub _Generate322KeySetKey221 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('S', 'S', $BondSymbol);
}

# Generate 322 keyset key 222 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 222 description: S#Cl
#
sub _Generate322KeySetKey222 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('S', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 223 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 223 description: S#P
#
sub _Generate322KeySetKey223 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('S', 'P', $BondSymbol);
}

# Generate 322 keyset key 224 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 224 description: S#F
#
sub _Generate322KeySetKey224 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('S', 'F', $BondSymbol);
}

# Generate 322 keyset key 225 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 225 description: S#Br
#
sub _Generate322KeySetKey225 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('S', 'Br', $BondSymbol);
}

# Generate 322 keyset key 226 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 226 description: S#Si
#
sub _Generate322KeySetKey226 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('S', 'Si', $BondSymbol);
}

# Generate 322 keyset key 227 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 227 description: S#I
#
sub _Generate322KeySetKey227 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('S', 'I', $BondSymbol);
}

# Generate 322 keyset key 228 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 228 description: S#X
#
sub _Generate322KeySetKey228 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('S', 'Z', $BondSymbol);
}

# Generate 322 keyset key 229 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 229 description: Cl#Cl
#
sub _Generate322KeySetKey229 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Cl', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 230 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 230 description: Cl#P
#
sub _Generate322KeySetKey230 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Cl', 'P', $BondSymbol);
}

# Generate 322 keyset key 231 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 231 description: Cl#F
#
sub _Generate322KeySetKey231 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Cl', 'F', $BondSymbol);
}

# Generate 322 keyset key 232 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 232 description: Cl#Br
#
sub _Generate322KeySetKey232 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Cl', 'Br', $BondSymbol);
}

# Generate 322 keyset key 233 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 233 description: Cl#Si
#
sub _Generate322KeySetKey233 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Cl', 'Si', $BondSymbol);
}

# Generate 322 keyset key 234 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 234 description: Cl#I
#
sub _Generate322KeySetKey234 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Cl', 'I', $BondSymbol);
}

# Generate 322 keyset key 235 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 235 description: Cl#X
#
sub _Generate322KeySetKey235 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Cl', 'Z', $BondSymbol);
}

# Generate 322 keyset key 236 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 236 description: P#P
#
sub _Generate322KeySetKey236 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('P', 'P', $BondSymbol);
}

# Generate 322 keyset key 237 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 237 description: P#F
#
sub _Generate322KeySetKey237 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('P', 'F', $BondSymbol);
}

# Generate 322 keyset key 238 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 238 description: P#Br
#
sub _Generate322KeySetKey238 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('P', 'Br', $BondSymbol);
}

# Generate 322 keyset key 239 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 239 description: P#Si
#
sub _Generate322KeySetKey239 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('P', 'Si', $BondSymbol);
}

# Generate 322 keyset key 240 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 240 description: P#I
#
sub _Generate322KeySetKey240 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('P', 'I', $BondSymbol);
}

# Generate 322 keyset key 241 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 241 description: P#X
#
sub _Generate322KeySetKey241 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('P', 'Z', $BondSymbol);
}

# Generate 322 keyset key 242 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 242 description: F#F
#
sub _Generate322KeySetKey242 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('F', 'F', $BondSymbol);
}

# Generate 322 keyset key 243 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 243 description: F#Br
#
sub _Generate322KeySetKey243 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('F', 'Br', $BondSymbol);
}

# Generate 322 keyset key 244 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 244 description: F#Si
#
sub _Generate322KeySetKey244 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('F', 'Si', $BondSymbol);
}

# Generate 322 keyset key 245 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 245 description: F#I
#
sub _Generate322KeySetKey245 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('F', 'I', $BondSymbol);
}

# Generate 322 keyset key 246 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 246 description: F#X
#
sub _Generate322KeySetKey246 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('F', 'Z', $BondSymbol);
}

# Generate 322 keyset key 247 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 247 description: Br#Br
#
sub _Generate322KeySetKey247 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Br', 'Br', $BondSymbol);
}

# Generate 322 keyset key 248 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 248 description: Br#Si
#
sub _Generate322KeySetKey248 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Br', 'Si', $BondSymbol);
}

# Generate 322 keyset key 249 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 249 description: Br#I
#
sub _Generate322KeySetKey249 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Br', 'I', $BondSymbol);
}

# Generate 322 keyset key 250 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 250 description: Br#X
#
sub _Generate322KeySetKey250 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Br', 'Z', $BondSymbol);
}

# Generate 322 keyset key 251 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 251 description: Si#Si
#
sub _Generate322KeySetKey251 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Si', 'Si', $BondSymbol);
}

# Generate 322 keyset key 252 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 252 description: Si#I
#
sub _Generate322KeySetKey252 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Si', 'I', $BondSymbol);
}

# Generate 322 keyset key 253 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 253 description: Si#X
#
sub _Generate322KeySetKey253 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Si', 'Z', $BondSymbol);
}

# Generate 322 keyset key 254 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 254 description: I#I
#
sub _Generate322KeySetKey254 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('I', 'I', $BondSymbol);
}

# Generate 322 keyset key 255 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 255 description: I#X
#
sub _Generate322KeySetKey255 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('I', 'Z', $BondSymbol);
}

# Generate 322 keyset key 256 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 256 description: X#X
#
sub _Generate322KeySetKey256 {
  my($This) = @_;
  my($BondSymbol) = '#';

  return $This->_DetectBondKeys('Z', 'Z', $BondSymbol);
}

# Generate 322 keyset key 257 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 257 description: C$C
#
sub _Generate322KeySetKey257 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'C', $BondSymbol);
}

# Generate 322 keyset key 258 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 258 description: C$N
#
sub _Generate322KeySetKey258 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'N', $BondSymbol);
}

# Generate 322 keyset key 259 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 259 description: C$O
#
sub _Generate322KeySetKey259 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'O', $BondSymbol);
}

# Generate 322 keyset key 260 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 260 description: C$S
#
sub _Generate322KeySetKey260 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'S', $BondSymbol);
}

# Generate 322 keyset key 261 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 261 description: C$Cl
#
sub _Generate322KeySetKey261 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 262 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 262 description: C$P
#
sub _Generate322KeySetKey262 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'P', $BondSymbol);
}

# Generate 322 keyset key 263 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 263 description: C$F
#
sub _Generate322KeySetKey263 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'F', $BondSymbol);
}

# Generate 322 keyset key 264 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 264 description: C$Br
#
sub _Generate322KeySetKey264 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'Br', $BondSymbol);
}

# Generate 322 keyset key 265 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 265 description: C$Si
#
sub _Generate322KeySetKey265 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'Si', $BondSymbol);
}

# Generate 322 keyset key 266 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 266 description: C$I
#
sub _Generate322KeySetKey266 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'I', $BondSymbol);
}

# Generate 322 keyset key 267 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 267 description: C$X
#
sub _Generate322KeySetKey267 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('C', 'Z', $BondSymbol);
}

# Generate 322 keyset key 268 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 268 description: N$N
#
sub _Generate322KeySetKey268 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'N', $BondSymbol);
}

# Generate 322 keyset key 269 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 269 description: N$O
#
sub _Generate322KeySetKey269 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'O', $BondSymbol);
}

# Generate 322 keyset key 270 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 270 description: N$S
#
sub _Generate322KeySetKey270 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'S', $BondSymbol);
}

# Generate 322 keyset key 271 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 271 description: N$Cl
#
sub _Generate322KeySetKey271 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 272 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 272 description: N$P
#
sub _Generate322KeySetKey272 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'P', $BondSymbol);
}

# Generate 322 keyset key 273 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 273 description: N$F
#
sub _Generate322KeySetKey273 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'F', $BondSymbol);
}

# Generate 322 keyset key 274 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 274 description: N$Br
#
sub _Generate322KeySetKey274 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'Br', $BondSymbol);
}

# Generate 322 keyset key 275 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 275 description: N$Si
#
sub _Generate322KeySetKey275 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'Si', $BondSymbol);
}

# Generate 322 keyset key 276 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 276 description: N$I
#
sub _Generate322KeySetKey276 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'I', $BondSymbol);
}

# Generate 322 keyset key 277 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 277 description: N$X
#
sub _Generate322KeySetKey277 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('N', 'Z', $BondSymbol);
}

# Generate 322 keyset key 278 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 278 description: O$O
#
sub _Generate322KeySetKey278 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'O', $BondSymbol);
}

# Generate 322 keyset key 279 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 279 description: O$S
#
sub _Generate322KeySetKey279 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'S', $BondSymbol);
}

# Generate 322 keyset key 280 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 280 description: O$Cl
#
sub _Generate322KeySetKey280 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 281 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 281 description: O$P
#
sub _Generate322KeySetKey281 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'P', $BondSymbol);
}

# Generate 322 keyset key 282 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 282 description: O$F
#
sub _Generate322KeySetKey282 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'F', $BondSymbol);
}

# Generate 322 keyset key 283 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 283 description: O$Br
#
sub _Generate322KeySetKey283 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'Br', $BondSymbol);
}

# Generate 322 keyset key 284 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 284 description: O$Si
#
sub _Generate322KeySetKey284 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'Si', $BondSymbol);
}

# Generate 322 keyset key 285 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 285 description: O$I
#
sub _Generate322KeySetKey285 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'I', $BondSymbol);
}

# Generate 322 keyset key 286 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 286 description: O$X
#
sub _Generate322KeySetKey286 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('O', 'Z', $BondSymbol);
}

# Generate 322 keyset key 287 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 287 description: S$S
#
sub _Generate322KeySetKey287 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('S', 'S', $BondSymbol);
}

# Generate 322 keyset key 288 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 288 description: S$Cl
#
sub _Generate322KeySetKey288 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('S', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 289 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 289 description: S$P
#
sub _Generate322KeySetKey289 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('S', 'P', $BondSymbol);
}

# Generate 322 keyset key 290 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 290 description: S$F
#
sub _Generate322KeySetKey290 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('S', 'F', $BondSymbol);
}

# Generate 322 keyset key 291 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 291 description: S$Br
#
sub _Generate322KeySetKey291 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('S', 'Br', $BondSymbol);
}

# Generate 322 keyset key 292 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 292 description: S$Si
#
sub _Generate322KeySetKey292 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('S', 'Si', $BondSymbol);
}

# Generate 322 keyset key 293 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 293 description: S$I
#
sub _Generate322KeySetKey293 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('S', 'I', $BondSymbol);
}

# Generate 322 keyset key 294 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 294 description: S$X
#
sub _Generate322KeySetKey294 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('S', 'Z', $BondSymbol);
}

# Generate 322 keyset key 295 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 295 description: Cl$Cl
#
sub _Generate322KeySetKey295 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Cl', 'Cl', $BondSymbol);
}

# Generate 322 keyset key 296 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 296 description: Cl$P
#
sub _Generate322KeySetKey296 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Cl', 'P', $BondSymbol);
}

# Generate 322 keyset key 297 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 297 description: Cl$F
#
sub _Generate322KeySetKey297 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Cl', 'F', $BondSymbol);
}

# Generate 322 keyset key 298 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 298 description: Cl$Br
#
sub _Generate322KeySetKey298 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Cl', 'Br', $BondSymbol);
}

# Generate 322 keyset key 299 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 299 description: Cl$Si
#
sub _Generate322KeySetKey299 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Cl', 'Si', $BondSymbol);
}

# Generate 322 keyset key 300 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 300 description: Cl$I
#
sub _Generate322KeySetKey300 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Cl', 'I', $BondSymbol);
}

# Generate 322 keyset key 301 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 301 description: Cl$X
#
sub _Generate322KeySetKey301 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Cl', 'Z', $BondSymbol);
}

# Generate 322 keyset key 302 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 302 description: P$P
#
sub _Generate322KeySetKey302 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('P', 'P', $BondSymbol);
}

# Generate 322 keyset key 303 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 303 description: P$F
#
sub _Generate322KeySetKey303 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('P', 'F', $BondSymbol);
}

# Generate 322 keyset key 304 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 304 description: P$Br
#
sub _Generate322KeySetKey304 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('P', 'Br', $BondSymbol);
}

# Generate 322 keyset key 305 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 305 description: P$Si
#
sub _Generate322KeySetKey305 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('P', 'Si', $BondSymbol);
}

# Generate 322 keyset key 306 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 306 description: P$I
#
sub _Generate322KeySetKey306 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('P', 'I', $BondSymbol);
}

# Generate 322 keyset key 307 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 307 description: P$X
#
sub _Generate322KeySetKey307 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('P', 'Z', $BondSymbol);
}

# Generate 322 keyset key 308 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 308 description: F$F
#
sub _Generate322KeySetKey308 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('F', 'F', $BondSymbol);
}

# Generate 322 keyset key 309 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 309 description: F$Br
#
sub _Generate322KeySetKey309 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('F', 'Br', $BondSymbol);
}

# Generate 322 keyset key 310 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 310 description: F$Si
#
sub _Generate322KeySetKey310 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('F', 'Si', $BondSymbol);
}

# Generate 322 keyset key 311 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 311 description: F$I
#
sub _Generate322KeySetKey311 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('F', 'I', $BondSymbol);
}

# Generate 322 keyset key 312 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 312 description: F$X
#
sub _Generate322KeySetKey312 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('F', 'Z', $BondSymbol);
}

# Generate 322 keyset key 313 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 313 description: Br$Br
#
sub _Generate322KeySetKey313 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Br', 'Br', $BondSymbol);
}

# Generate 322 keyset key 314 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 314 description: Br$Si
#
sub _Generate322KeySetKey314 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Br', 'Si', $BondSymbol);
}

# Generate 322 keyset key 315 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 315 description: Br$I
#
sub _Generate322KeySetKey315 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Br', 'I', $BondSymbol);
}

# Generate 322 keyset key 316 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 316 description: Br$X
#
sub _Generate322KeySetKey316 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Br', 'Z', $BondSymbol);
}

# Generate 322 keyset key 317 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 317 description: Si$Si
#
sub _Generate322KeySetKey317 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Si', 'Si', $BondSymbol);
}

# Generate 322 keyset key 318 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 318 description: Si$I
#
sub _Generate322KeySetKey318 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Si', 'I', $BondSymbol);
}

# Generate 322 keyset key 319 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 319 description: Si$X
#
sub _Generate322KeySetKey319 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Si', 'Z', $BondSymbol);
}

# Generate 322 keyset key 320 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 320 description: I$I
#
sub _Generate322KeySetKey320 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('I', 'I', $BondSymbol);
}

# Generate 322 keyset key 321 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 321 description: I$X
#
sub _Generate322KeySetKey321 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('I', 'Z', $BondSymbol);
}

# Generate 322 keyset key 322 value as 1/0 indicating its presence/absence or
# count of its presence in a molecule.
#
# Key 322 description: X$X
#
sub _Generate322KeySetKey322 {
  my($This) = @_;
  my($BondSymbol) = '$';

  return $This->_DetectBondKeys('Z', 'Z', $BondSymbol);
}

# A : Any valid perodic table elemnet symbol
sub _IsAtom {
  my($This, $Atom) = @_;

  return $Atom->GetAtomicNumber() ? 1 : 0;
}

# Q  : Hetro atoms; any non-C or non-H atom
sub _IsHeteroAtom {
  my($This, $Atom) = @_;

  return ($Atom->GetAtomicNumber() =~ /^(1|6)$/) ? 0 : 1;
}

# X  : Halogens; F, Cl, Br, I
sub _IsHalogenAtom {
  my($This, $Atom) = @_;

  return ($Atom->GetAtomicNumber() =~ /^(9|17|35|53)$/) ? 1 : 0;
}

# Z  : Others; other than H, C, N, O, Si, P, S, F, Cl, Br, I
sub _IsOtherAtom {
  my($This, $Atom) = @_;

  return ($Atom->GetAtomicNumber() =~ /^(1|6|7|8|9|14|15|16|17|35|53)$/) ? 0 : 1;
}

# Detect atom keys like Cl, Br and so on...
#
sub _DetectAtomKeys {
  my($This, $AtomSymbol, $MinKeyCount, $IsInRing, $MinHydrogenCount) = @_;
  my($Atom, $KeyValue);

  $KeyValue = 0;
  ATOM: for $Atom (@{$This->{Atoms}}) {
    if (!$This->_DoesAtomMatchesSymbol($Atom, $AtomSymbol)) {
      next ATOM;
    }
    if (defined($IsInRing) && $IsInRing && !$Atom->IsInRing()) {
      next ATOM;
    }
    if (defined $MinHydrogenCount) {
      if (!$This->_DoesAtomMinHydrogenCountMatch($Atom, $MinHydrogenCount)) {
	next ATOM;
      }
    }
    $KeyValue++;
    if (defined($MinKeyCount) && $KeyValue < $MinKeyCount) {
      next ATOM;
    }
    if ($This->{KeyBits}) {
      $KeyValue = 1;
      last ATOM;
    }
  }
  return $KeyValue;
}

# Detect bond keys like S-S, N-O and so on...
#
sub _DetectBondKeys {
  my($This, $AtomSymbol1, $AtomSymbol2, $BondSymbol, $MinKeyCount, $Atom1MinHydrogenCount, $Atom2MinHydrogenCount) = @_;
  my($Atom1, $Atom2, $Bond, $KeyValue, $MatchSpecifiedAtomOrder);

  $MatchSpecifiedAtomOrder = 0;

  $KeyValue = 0;
  BOND: for $Bond (@{$This->{Bonds}}) {
    ($Atom1, $Atom2) = $Bond->GetAtoms();
    if (!$This->_DoBondAtomsMatchBondSymbols($Atom1, $Atom2, $AtomSymbol1, $AtomSymbol2, $BondSymbol, $MatchSpecifiedAtomOrder, $Atom1MinHydrogenCount, $Atom2MinHydrogenCount)) {
      next BOND;
    }
    $KeyValue++;
    if (defined($MinKeyCount) && $KeyValue < $MinKeyCount) {
      next BOND;
    }
    if ($This->{KeyBits}) {
      $KeyValue = 1;
      last BOND;
    }
  }
  return $KeyValue;
}

# Detect atom neighborhood keys like ON(C)C, OC(O)O and so on.
#
sub _DetectAtomNeighborhoodKeys {
  my($This, $CentralAtomSymbol, $NbrAtomSymbolsRef, $NbrBondSymbolsRef, $MinKeyCount, $CentralAtomMinHydrogenCount, $NbrAtomMinHydrogenCountRef) = @_;
  my($KeyValue, $CentralAtom);

  $KeyValue = 0;

  CENTRALATOM: for $CentralAtom (@{$This->{Atoms}}) {
    if (!$This->_DoesAtomNeighborhoodMatch($CentralAtom, $CentralAtomSymbol, $NbrAtomSymbolsRef, $NbrBondSymbolsRef, $CentralAtomMinHydrogenCount, $NbrAtomMinHydrogenCountRef)) {
      next CENTRALATOM;
    }
    $KeyValue++;
    if (defined($MinKeyCount) && $KeyValue < $MinKeyCount) {
      next CENTRALATOM;
    }
    if ($This->{KeyBits}) {
      $KeyValue = 1;
      last CENTRALATOM;
    }
  }
  return $KeyValue;
}

# Detect bond neighborhood keys like A%Anot%A%A and so on.
#
sub _DetectBondNeighborhoodKeys {
  my($This, $BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, $NbrAtomSymbolsRef, $NbrBondSymbolsRef, $MinKeyCount, $BondAtomMinHydrogenCountRef, $NbrsMinHydrogenCountRef) = @_;
  my($BondAtomIndex, $BondAtom1, $BondAtom2, $MatchedBondAtom1, $MatchedBondAtom2, $BondAtom, $Bond, $KeyValue, $BondAtomSymbol, $MatchSpecifiedAtomOrder, $BondAtom1MinHydrogenCount, $BondAtom2MinHydrogenCount, $MinHydrogenCount, @NbrsToExcludeFromMatch, @NbrAtomSymbols, @NbrBondSymbols, @NbrMinHydrogenCount, );

  $MatchSpecifiedAtomOrder = 1;
  ($BondAtom1MinHydrogenCount, $BondAtom2MinHydrogenCount) = defined($BondAtomMinHydrogenCountRef) ? ( @{$BondAtomMinHydrogenCountRef} ) : (undef, undef);

  $KeyValue = 0;
  BOND: for $Bond (@{$This->{Bonds}}) {
    ($BondAtom1, $BondAtom2) = $Bond->GetAtoms();

    # Match bond first...
    if ($This->_DoBondAtomsMatchBondSymbols($BondAtom1, $BondAtom2, $BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, $MatchSpecifiedAtomOrder, $BondAtom1MinHydrogenCount, $BondAtom2MinHydrogenCount)) {
      ($MatchedBondAtom1, $MatchedBondAtom2) = ($BondAtom1, $BondAtom2);
    }
    elsif ($This->_DoBondAtomsMatchBondSymbols($BondAtom2, $BondAtom1, $BondAtomSymbol1, $BondAtomSymbol2, $BondSymbol, $MatchSpecifiedAtomOrder, $BondAtom1MinHydrogenCount, $BondAtom2MinHydrogenCount)) {
      ($MatchedBondAtom1, $MatchedBondAtom2) = ($BondAtom2, $BondAtom1);
    }
    else {
      next BOND;
    }
    # Match neighbors of bonded atoms...
    for $BondAtomIndex (0 .. 1) {
      $MinHydrogenCount = undef;
      @NbrsToExcludeFromMatch = ();

      if ($BondAtomIndex == 0) {
	$BondAtom = $MatchedBondAtom1;
	$BondAtomSymbol = $BondAtomSymbol1;
	push @NbrsToExcludeFromMatch, $MatchedBondAtom2;
      }
      elsif ($BondAtomIndex == 1) {
	$BondAtom = $MatchedBondAtom2;
	$BondAtomSymbol = $BondAtomSymbol2;
	push @NbrsToExcludeFromMatch, $MatchedBondAtom1;
      }

      @NbrAtomSymbols = (defined($NbrAtomSymbolsRef) && defined($NbrAtomSymbolsRef->[$BondAtomIndex])) ? @{$NbrAtomSymbolsRef->[$BondAtomIndex]} : ();
      @NbrBondSymbols = (defined($NbrBondSymbolsRef) && defined($NbrBondSymbolsRef->[$BondAtomIndex]) ) ? @{$NbrBondSymbolsRef->[$BondAtomIndex]} : ();
      @NbrMinHydrogenCount = (defined($NbrsMinHydrogenCountRef) && defined($NbrsMinHydrogenCountRef->[$BondAtomIndex]) ) ? @{$NbrsMinHydrogenCountRef->[$BondAtomIndex]} : ();
      if (!$This->_DoesAtomNeighborhoodMatch($BondAtom, $BondAtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinHydrogenCount, \@NbrMinHydrogenCount, \@NbrsToExcludeFromMatch)) {
	next BOND;
      }
    }
    $KeyValue++;
    if (defined($MinKeyCount) && $KeyValue < $MinKeyCount) {
      next BOND;
    }
    if ($This->{KeyBits}) {
      $KeyValue = 1;
      last BOND;
    }
  }
  return $KeyValue;
}

# Detect extended atom neighborhood keys like QHAQH, QHAAQH, and so on...
#
sub _DetectExtendedAtomNeighborhoodKeys {
  my($This, $CentralAtomsSymbolsRef, $CentralAtomsBondSymbolsRef, $CentralAtomsMinHydrogenCountRef, $MinKeyCount, $NbrAtomSymbolsRef, $NbrBondSymbolsRef, $NbrsMinHydrogenCountRef) = @_;
  my($KeyValue, $Molecule, $FirstCentralAtomIndex, $LastCentralAtomIndex, $NumOfCentralAtoms);

  $KeyValue = 0;

  $Molecule = $This->GetMolecule();
  $NumOfCentralAtoms = @{$CentralAtomsSymbolsRef};
  $FirstCentralAtomIndex = 0;
  $LastCentralAtomIndex = $NumOfCentralAtoms - 1;

  # Retrieve first central atom information...
  my($FirstCentralAtomSymbol, $FirstCentralAtomMinHydrogenCount);
  $FirstCentralAtomSymbol = $CentralAtomsSymbolsRef->[$FirstCentralAtomIndex];
  $FirstCentralAtomMinHydrogenCount = defined($CentralAtomsMinHydrogenCountRef) ? $CentralAtomsMinHydrogenCountRef->[$FirstCentralAtomIndex] : undef;

  # Retrieve last central atom information...
  my($LastCentralAtomSymbol, $LastCentralAtomMinHydrogenCount);
  $LastCentralAtomSymbol = $CentralAtomsSymbolsRef->[$LastCentralAtomIndex];
  $LastCentralAtomMinHydrogenCount = defined($CentralAtomsMinHydrogenCountRef) ? $CentralAtomsMinHydrogenCountRef->[$LastCentralAtomIndex] : undef;

  my($Atom, $AtomPathRef, $AtomPathsRef, $FirstAtomIndex, $LastAtomIndex, $AtomIndex, $FirstPathAtom, $LastPathAtom, $FirstPathAtomID, $LastPathAtomID, $DetectedPathID, $PathAtom, $NextPathAtom, $PreviousPathAtom, $AtomSymbol, $NextAtomSymbol, $BondSymbol, $MatchSpecifiedAtomOrder, $MinHydrogenCount, @NbrsToExcludeFromMatch, @NbrAtomSymbols, @NbrBondSymbols, @NbrMinHydrogenCount, %AlreadyDetectedPaths);

  # Go over all the atoms...
  #
  ATOM: for $Atom (@{$This->{Atoms}}) {
    # Match first central atom...
    if (!$This->_DoesAtomNeighborhoodMatch($Atom, $FirstCentralAtomSymbol, undef, undef, $FirstCentralAtomMinHydrogenCount, undef)) {
      next ATOM;
    }
    # Get atom paths starting from matched central atom with length equal to NumOfCentralAtoms...
    #
    $AtomPathsRef = $Molecule->GetAllAtomPathsStartingAtWithLength($Atom, $NumOfCentralAtoms);
    if (!(defined($AtomPathsRef) && @{$AtomPathsRef})) {
      next ATOM;
    }
    ATOMPATH: for $AtomPathRef (@{$AtomPathsRef}) {
      $FirstAtomIndex = 0;
      $FirstPathAtom = $AtomPathRef->[$FirstAtomIndex];
      $LastAtomIndex = @{$AtomPathRef} - 1;
      $LastPathAtom = $AtomPathRef->[$LastAtomIndex];

      # Match last central atom to the last atom in path...
      if (!$This->_DoesAtomNeighborhoodMatch($LastPathAtom, $LastCentralAtomSymbol, undef, undef, $LastCentralAtomMinHydrogenCount, undef)) {
	next ATOMPATH;
      }

      # Match other path atoms with central atoms..
      for $AtomIndex ($FirstAtomIndex .. $LastAtomIndex) {
	$PathAtom = $AtomPathRef->[$AtomIndex];
	$AtomSymbol = $CentralAtomsSymbolsRef->[$AtomIndex];
	$MinHydrogenCount = defined($CentralAtomsMinHydrogenCountRef) ? $CentralAtomsMinHydrogenCountRef->[$AtomIndex] : undef;

	@NbrsToExcludeFromMatch = ();
	if ($AtomIndex == $FirstAtomIndex) {
	  $NextPathAtom = $AtomPathRef->[$AtomIndex + 1]; $PreviousPathAtom = undef;
	  push @NbrsToExcludeFromMatch, $NextPathAtom;
	}
	elsif ($AtomIndex == $LastAtomIndex) {
	  $NextPathAtom = undef; $PreviousPathAtom = $AtomPathRef->[$AtomIndex - 1];
	  push @NbrsToExcludeFromMatch, $PreviousPathAtom;
	}
	else {
	  $NextPathAtom = $AtomPathRef->[$AtomIndex + 1]; $PreviousPathAtom = $AtomPathRef->[$AtomIndex - 1];
	  push @NbrsToExcludeFromMatch, ($PreviousPathAtom, $NextPathAtom);
	}

	@NbrAtomSymbols = (defined($NbrAtomSymbolsRef) && defined($NbrAtomSymbolsRef->[$AtomIndex])) ? @{$NbrAtomSymbolsRef->[$AtomIndex]} : ();
	@NbrBondSymbols = (defined($NbrBondSymbolsRef) && defined($NbrBondSymbolsRef->[$AtomIndex]) ) ? @{$NbrBondSymbolsRef->[$AtomIndex]} : ();
	@NbrMinHydrogenCount = (defined($NbrsMinHydrogenCountRef) && defined($NbrsMinHydrogenCountRef->[$AtomIndex]) ) ? @{$NbrsMinHydrogenCountRef->[$AtomIndex]} : ();

	if (!$This->_DoesAtomNeighborhoodMatch($PathAtom, $AtomSymbol, \@NbrAtomSymbols, \@NbrBondSymbols, $MinHydrogenCount, \@NbrMinHydrogenCount, \@NbrsToExcludeFromMatch)) {
	  next ATOMPATH;
	}
	# Match path bond symbols...
	if (defined($CentralAtomsBondSymbolsRef) && ($AtomIndex < $LastAtomIndex)) {
	  $NextAtomSymbol = $CentralAtomsSymbolsRef->[$AtomIndex + 1];
	  $BondSymbol = $CentralAtomsBondSymbolsRef->[$AtomIndex];
	  $MatchSpecifiedAtomOrder = 1;
	  if (!$This->_DoBondAtomsMatchBondSymbols($PathAtom, $NextPathAtom, $AtomSymbol, $NextAtomSymbol, $BondSymbol, $MatchSpecifiedAtomOrder)) {
	    next ATOMPATH;
	  }
	}
      }
      # Keep track of the first and last atom ID in the matched path to avoid double counting of paths...
      if (defined($MinKeyCount) || !$This->{KeyBits}) {
	$FirstPathAtomID = $FirstPathAtom->GetID(); $LastPathAtomID = $LastPathAtom->GetID();
	$DetectedPathID = ($FirstPathAtomID < $LastPathAtomID) ? "${FirstPathAtomID}-${LastPathAtomID}" : "${LastPathAtomID}-${FirstPathAtomID}";
	if (exists $AlreadyDetectedPaths{$DetectedPathID}) {
	  $AlreadyDetectedPaths{$DetectedPathID} += 1;
	  next ATOMPATH;
	}
	$AlreadyDetectedPaths{$DetectedPathID} = 1;
      }

      $KeyValue++;
      if (defined($MinKeyCount) && $KeyValue < $MinKeyCount) {
	next ATOMPATH;
      }
      if ($This->{KeyBits}) {
	$KeyValue = 1;
	last ATOM;
      }
    }
  }
  return $KeyValue;
}

# Go over the atoms attached to central atom and match 'em against specified
# neighborhood atom symbol and bond symbols...
#
sub _DoesAtomNeighborhoodMatch {
  my($This, $CentralAtom, $CentralAtomSymbol, $NbrAtomSymbolsRef, $NbrBondSymbolsRef, $CentralAtomMinHydrogenCount, $NbrAtomMinHydrogenCountRef, $NbrsToExcludeRef) = @_;

  # Match central atom first...
  if (!$This->_DoesAtomMatchesSymbol($CentralAtom, $CentralAtomSymbol)) {
    return 0;
  }
  if (defined $CentralAtomMinHydrogenCount) {
    if (!$This->_DoesAtomMinHydrogenCountMatch($CentralAtom, $CentralAtomMinHydrogenCount)) {
      return 0;
    }
  }
  if (!defined $NbrAtomSymbolsRef) {
    # No neighbors to match...
    return 1;
  }

  # Match neighbors...
  my($NbrAtom, $Index, $NbrAtomSymbol, $NbrBondSymbol, $NbrAtomMinHydrogenCount, $NbrAtomMatchCount, $MinNbrAtomMatchCount, $MatchSpecifiedAtomOrder, @CentralAtomNeighbors, %NbrAtomAlreadyMatchedMap);

  $MinNbrAtomMatchCount = @$NbrAtomSymbolsRef;
  if (!$MinNbrAtomMatchCount) {
    # No neighbors to match...
    return 1;
  }

  $NbrAtomMatchCount = 0;

  %NbrAtomAlreadyMatchedMap = ();
  $MatchSpecifiedAtomOrder = 1;

  @CentralAtomNeighbors = ();
  if (defined($NbrsToExcludeRef) && @{$NbrsToExcludeRef}) {
    push @CentralAtomNeighbors, $CentralAtom->GetNeighbors(@{$NbrsToExcludeRef});
  }
  else {
    push @CentralAtomNeighbors, $CentralAtom->GetNeighbors();
  }

  NBRATOM: for $NbrAtom (@CentralAtomNeighbors) {
    NBRATOMSYMBOL: for $Index (0 .. ($MinNbrAtomMatchCount - 1)) {
      if (exists $NbrAtomAlreadyMatchedMap{$Index}) {
	next NBRATOMSYMBOL;
      }
      $NbrAtomSymbol = $NbrAtomSymbolsRef->[$Index];
      $NbrBondSymbol = $NbrBondSymbolsRef->[$Index];
      if (!$This->_DoBondAtomsMatchBondSymbols($CentralAtom, $NbrAtom, $CentralAtomSymbol, $NbrAtomSymbol, $NbrBondSymbol, $MatchSpecifiedAtomOrder)) {
	next NBRATOMSYMBOL;
      }

      if (defined($NbrAtomMinHydrogenCountRef) && $NbrAtomMinHydrogenCountRef->[$Index]) {
	$NbrAtomMinHydrogenCount = $NbrAtomMinHydrogenCountRef->[$Index];
	if (!$This->_DoesAtomMinHydrogenCountMatch($NbrAtom, $NbrAtomMinHydrogenCount)) {
	  next NBRATOMSYMBOL;
	}
      }
      $NbrAtomAlreadyMatchedMap{$Index} = $Index;
      $NbrAtomMatchCount++;

      if ($NbrAtomMatchCount == $MinNbrAtomMatchCount) {
	last NBRATOM;
      }
      next NBRATOM;
    }
  }

  return ($NbrAtomMatchCount == $MinNbrAtomMatchCount) ? 1 : 0;
}

# Checks whether bond atoms match bond symbols...
#
sub _DoBondAtomsMatchBondSymbols {
  my($This, $Atom1, $Atom2, $AtomSymbol1, $AtomSymbol2, $BondSymbol, $MatchSpecifiedAtomOrder, $Atom1MinHydrogenCount, $Atom2MinHydrogenCount) = @_;
  my($Status, $ReverseMinHydrogenCountMatch);

  $ReverseMinHydrogenCountMatch = 0;

  if (defined($MatchSpecifiedAtomOrder) && $MatchSpecifiedAtomOrder) {
    if (!($This->_DoesAtomMatchesSymbol($Atom1, $AtomSymbol1) && $This->_DoesAtomMatchesSymbol($Atom2, $AtomSymbol2))) {
      return 0;
    }
  }
  else {
    if ($This->_DoesAtomMatchesSymbol($Atom1, $AtomSymbol1) && $This->_DoesAtomMatchesSymbol($Atom2, $AtomSymbol2)) {
      $ReverseMinHydrogenCountMatch = 0;
    }
    elsif ($This->_DoesAtomMatchesSymbol($Atom1, $AtomSymbol2) && $This->_DoesAtomMatchesSymbol($Atom2, $AtomSymbol1)) {
      $ReverseMinHydrogenCountMatch = 1;
    }
    else {
      return 0;
    }
  }

  # Match any hydrogen count...
  if (defined($Atom1MinHydrogenCount) || defined($Atom2MinHydrogenCount)) {
    my($MinHydrogenCountMatchAtom1, $MinHydrogenCountMatchAtom2);

    ($MinHydrogenCountMatchAtom1, $MinHydrogenCountMatchAtom2) = $ReverseMinHydrogenCountMatch ? ($Atom2, $Atom1) : ($Atom1, $Atom2);
    if (defined $Atom1MinHydrogenCount ) {
      if (!$This->_DoesAtomMinHydrogenCountMatch($MinHydrogenCountMatchAtom1, $Atom1MinHydrogenCount)) {
	return 0;
      }
    }
    if (defined $Atom2MinHydrogenCount ) {
      if (!$This->_DoesAtomMinHydrogenCountMatch($MinHydrogenCountMatchAtom2, $Atom2MinHydrogenCount)) {
	return 0;
      }
    }
  }

  if (defined($BondSymbol) && $BondSymbol) {
    my($Bond);
    $Bond = $Atom1->GetBondToAtom($Atom2);
    if (!$This->_DoesBondMatchesSymbol($Bond, $BondSymbol)) {
      return 0;
    }
  }
  return 1;
}

# Match both implicit and explicit hydrogens on central atom...
sub _DoesAtomMinHydrogenCountMatch {
  my($This, $Atom, $MinHydrogenCount) = @_;

  if (!(defined($MinHydrogenCount) && $MinHydrogenCount)) {
    return 0;
  }
  return ($Atom->GetNumOfHydrogens() <  $MinHydrogenCount) ? 0 : 1;
}

# Checks whether atom matches supported symbol...
#
sub _DoesAtomMatchesSymbol {
  my($This, $Atom, $Symbol) = @_;
  my($Status);

  $Status = 0;
  SYMBOL: {
    if ($Symbol =~ /^Q$/i) { $Status = $This->_IsHeteroAtom($Atom) ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^X$/i) { $Status = $This->_IsHalogenAtom($Atom) ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^Z$/i) { $Status = $This->_IsOtherAtom($Atom) ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^A$/i) { $Status = $This->_IsAtom($Atom) ? 1 : 0; last SYMBOL; }
    $Status = ($Atom->GetAtomSymbol() =~ /^$Symbol$/i) ? 1 : 0;
  }
  return $Status;
}

# Checks whether bond matches supported symbol...
#
sub _DoesBondMatchesSymbol {
  my($This, $Bond, $Symbol) = @_;
  my($Status, $BondOrder);

  $Status = 0;
  SYMBOL: {
    if ($Symbol =~ /^(1|-)$/i) { $Status = $Bond->IsSingle() ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^(2|=)$/i) { $Status = $Bond->IsDouble() ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^(3|#|T)$/i) { $Status = $Bond->IsTriple() ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^(1.5|%)$/i) { $Status = $Bond->IsAromatic() ? 1 : 0; last SYMBOL; }

    if ($Symbol =~ /^\~$/i) { $Status = ($Bond->IsSingle() || $Bond->IsDouble()) ? 1 : 0; last SYMBOL; }

    if ($Symbol =~ /^\$$/i) { $Status = $Bond->IsInRing() ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^\!$/i) { $Status = $Bond->IsInRing() ? 0 : 1; last SYMBOL; }

    if ($Symbol =~ /^(\$-)$/i) { $Status = ($Bond->IsInRing() && $Bond->IsSingle()) ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^(\$=)$/i) { $Status = ($Bond->IsInRing() && $Bond->IsDouble()) ? 1 : 0; last SYMBOL; }
    if ($Symbol =~ /^(\$#|\$T)$/i) { $Status = ($Bond->IsInRing() && $Bond->IsTriple()) ? 1 : 0; last SYMBOL; }

    if ($Symbol =~ /^(not%)$/i) { $Status = $Bond->IsAromatic() ? 0 : 1; last SYMBOL; }
    if ($Symbol =~ /^(not%not-)$/i) { $Status = $Bond->IsAromatic() ? 0 : ($Bond->IsSingle() ? 0 : 1); last SYMBOL; }
    if ($Symbol =~ /^(not%not=)$/i) { $Status = $Bond->IsAromatic() ? 0 : ($Bond->IsDouble() ? 0 : 1); last SYMBOL; }

    $Status = 0;
  }
  return $Status;
}

# Cache  appropriate molecule data...
#
sub _SetupMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = $This->GetMolecule()->GetAtoms();
  @{$This->{Bonds}} = $This->GetMolecule()->GetBonds();

  return $This;
}

# Clear cached molecule data...
#
sub _ClearMoleculeDataCache {
  my($This) = @_;

  @{$This->{Atoms}} = ();
  @{$This->{Bonds}} = ();

  return $This;
}

# Return a string containg data for MACCSKeys object...
sub StringifyMACCSKeys {
  my($This) = @_;
  my($MACCSKeysString);

  # Type of Keys...
  $MACCSKeysString = "Type: $This->{Type}; Size: $This->{Size}";

  if ($This->{Type} =~ /^MACCSKeyBits$/i) {
    $MACCSKeysString .= "; FingerprintsBitVector: < $This->{FingerprintsBitVector} >";
  }
  elsif ($This->{Type} =~ /^MACCSKeyCount$/i) {
    $MACCSKeysString .= "; FingerprintsVector: < $This->{FingerprintsVector} >";
  }

  return $MACCSKeysString;
}

1;

__END__

=head1 NAME

MACCSKeys

=head1 SYNOPSIS

use Fingerprints::MACCSKeys;

use Fingerprints::MACCSKeys qw(:all);

=head1 DESCRIPTION

B<MACCSKeys> [ Ref 45-47 ] class provides the following methods:

new, GenerateFingerprints, GenerateMACCSKeys, GetDescription, SetSize, SetType,
StringifyMACCSKeys

B<MACCSKeys> is derived from B<Fingerprints> class which in turn is  derived from
B<ObjectProperty> base class that provides methods not explicitly defined in B<MACCSKeys>,
B<Fingerprints> or B<ObjectProperty> classes using Perl's AUTOLOAD functionality. These
methods are generated on-the-fly for a specified object property:

    Set<PropertyName>(<PropertyValue>);
    $PropertyValue = Get<PropertyName>();
    Delete<PropertyName>();

For each MACCS (Molecular ACCess System) keys definition, atoms are processed to
determine their membership to the key and the appropriate molecular fingerprints strings
are generated. An atom can belong to multiple MACCS keys.

For I<MACCSKeyBits> value of B<Type> option, a fingerprint bit-vector string containing
zeros and ones is generated and for I<MACCSKeyCount> value, a fingerprint vector string
corresponding to number of MACCS keys [ Ref 45-47 ] is generated.

I<MACCSKeyBits or MACCSKeyCount> values for B<Type> along with two possible
I<166 | 322>  values of B<Size> supports generation of four different types of MACCS
keys fingerprint: I<MACCS166KeyBits, MACCS166KeyCount, MACCS322KeyBits, MACCS322KeyCount>.

The current release of MayaChemTools generates the following types of MACCS keys
fingerprints bit-vector and vector strings:

    FingerprintsBitVector;MACCSKeyBits;166;BinaryString;Ascending;00000000
    0000000000000000000000000000000001001000010010000000010010000000011100
    0100101010111100011011000100110110000011011110100110111111111111011111
    11111111111110111000

    FingerprintsBitVector;MACCSKeyBits;166;HexadecimalString;Ascending;000
    000000021210210e845f8d8c60b79dffbffffd1

    FingerprintsBitVector;MACCSKeyBits;322;BinaryString;Ascending;11101011
    1110011111100101111111000111101100110000000000000011100010000000000000
    0000000000000000000000000000000000000000000000101000000000000000000000
    0000000000000000000000000000000000000000000000000000000000000000000000
    0000000000000000000000000000000000000011000000000000000000000000000000
    0000000000000000000000000000000000000000

    FingerprintsBitVector;MACCSKeyBits;322;HexadecimalString;Ascending;7d7
    e7af3edc000c1100000000000000500000000000000000000000000000000300000000
    000000000

    FingerprintsVector;MACCSKeyCount;166;OrderedNumericalValues;ValuesStri
    ng;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 3 0 0 0 0 4 0 0 2 0 0 0 0 0 0 0 0 2 0 0 2 0 0 0 0
    0 0 0 0 1 1 8 0 0 0 1 0 0 1 0 1 0 1 0 3 1 3 1 0 0 0 1 2 0 11 1 0 0 0
    5 0 0 1 2 0 1 1 0 0 0 0 0 1 1 0 1 1 1 1 0 4 0 0 1 1 0 4 6 1 1 1 2 1 1
    3 5 2 2 0 5 3 5 1 1 2 5 1 2 1 2 4 8 3 5 5 2 2 0 3 5 4 1

    FingerprintsVector;MACCSKeyCount;322;OrderedNumericalValues;ValuesStri
    ng;14 8 2 0 2 0 4 4 2 1 4 0 0 2 5 10 5 2 1 0 0 2 0 5 13 3 28 5 5 3 0 0
    0 4 2 1 1 0 1 1 0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 22 5 3 0 0 0 1 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11 0 2 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...

=head2 METHODS

=over 4

=item B<new>

    $NewMACCSKeys = new MACCSKeys(%NamesAndValues);

Using specified I<MACCSKeys> property names and values hash, B<new> method creates a new object
and returns a reference to newly created B<PathLengthFingerprints> object. By default, the
following properties are initialized:

    Molecule = '';
    Type = ''
    Size = ''

Examples:

    $MACCSKeys = new MACCSKeys('Molecule' => $Molecule,
                               'Type' => 'MACCSKeyBits',
                               'Size' => 166);

    $MACCSKeys = new MACCSKeys('Molecule' => $Molecule,
                               'Type' => 'MACCSKeyCount',
                               'Size' => 166);

    $MACCSKeys = new MACCSKeys('Molecule' => $Molecule,
                               'Type' => 'MACCSKeyBit',
                               'Size' => 322);

    $MACCSKeys = new MACCSKeys('Molecule' => $Molecule,
                               'Type' => 'MACCSKeyCount',
                               'Size' => 322);

    $MACCSKeys->GenerateMACCSKeys();
    print "$MACCSKeys\n";

=item B<GetDescription>

    $Description = $MACCSKeys->GetDescription();

Returns a string containing description of MACCS keys fingerprints.

=item B<GenerateMACCSKeys or GenerateFingerprints>

    $MACCSKeys = $MACCSKeys->GenerateMACCSKeys();

Generates MACCS keys fingerprints and returns I<MACCSKeys>.

For I<MACCSKeyBits> value of B<Type>, a fingerprint bit-vector string containing
zeros and ones is generated and for I<MACCSKeyCount> value, a fingerprint vector string
corresponding to number of MACCS keys is generated.

I<MACCSKeyBits or MACCSKeyCount> values for B<Type> option along with two possible
I<166 | 322>  values of B<Size> supports generation of four different types of MACCS
keys fingerprint: I<MACCS166KeyBits, MACCS166KeyCount, MACCS322KeyBits, MACCS322KeyCount>.

Definition of MACCS keys uses the following atom and bond symbols to define atom and
bond environments:

    Atom symbols for 166 keys [ Ref 47 ]:

    A : Any valid periodic table element symbol
    Q  : Hetro atoms; any non-C or non-H atom
    X  : Halogens; F, Cl, Br, I
    Z  : Others; other than H, C, N, O, Si, P, S, F, Cl, Br, I

    Atom symbols for 322 keys [ Ref 46 ]:

    A : Any valid periodic table element symbol
    Q  : Hetro atoms; any non-C or non-H atom
    X  : Others; other than H, C, N, O, Si, P, S, F, Cl, Br, I
    Z is neither defined nor used

    Bond types:

    -  : Single
    =  : Double
    T  : Triple
    #  : Triple
    ~  : Single or double query bond
    %  : An aromatic query bond

    None : Any bond type; no explicit bond specified

    $  : Ring bond; $ before a bond type specifies ring bond
    !  : Chain or non-ring bond; ! before a bond type specifies chain bond

    @  : A ring linkage and the number following it specifies the
         atoms position in the line, thus @1 means linked back to the first
         atom in the list.

    Aromatic: Kekule or Arom5

    Kekule: Bonds in 6-membered rings with alternate single/double bonds
            or perimeter bonds
    Arom5:  Bonds in 5-membered rings with two double bonds and a hetro
            atom at the apex of the ring.

MACCS 166 keys [ Ref 45-47 ] are defined as follows:

    Key Description

    1	ISOTOPE
    2	103 < ATOMIC NO. < 256
    3	GROUP IVA,VA,VIA PERIODS 4-6 (Ge...)
    4	ACTINIDE
    5	GROUP IIIB,IVB (Sc...)
    6	LANTHANIDE
    7	GROUP VB,VIB,VIIB (V...)
    8	QAAA@1
    9	GROUP VIII (Fe...)
    10	GROUP IIA (ALKALINE EARTH)
    11	4M RING
    12	GROUP IB,IIB (Cu...)
    13	ON(C)C
    14	S-S
    15	OC(O)O
    16	QAA@1
    17	CTC
    18	GROUP IIIA (B...)
    19	7M RING
    20	SI
    21	C=C(Q)Q
    22	3M RING
    23	NC(O)O
    24	N-O
    25	NC(N)N
    26	C$=C($A)$A
    27	I
    28	QCH2Q
    29	P
    30	CQ(C)(C)A
    31	QX
    32	CSN
    33	NS
    34	CH2=A
    35	GROUP IA (ALKALI METAL)
    36	S HETEROCYCLE
    37	NC(O)N
    38	NC(C)N
    39	OS(O)O
    40	S-O
    41	CTN
    42	F
    43	QHAQH
    44	OTHER
    45	C=CN
    46	BR
    47	SAN
    48	OQ(O)O
    49	CHARGE
    50	C=C(C)C
    51	CSO
    52	NN
    53	QHAAAQH
    54	QHAAQH
    55	OSO
    56	ON(O)C
    57	O HETEROCYCLE
    58	QSQ
    59	Snot%A%A
    60	S=O
    61	AS(A)A
    62	A$A!A$A
    63	N=O
    64	A$A!S
    65	C%N
    66	CC(C)(C)A
    67	QS
    68	QHQH (&...)
    69	QQH
    70	QNQ
    71	NO
    72	OAAO
    73	S=A
    74	CH3ACH3
    75	A!N$A
    76	C=C(A)A
    77	NAN
    78	C=N
    79	NAAN
    80	NAAAN
    81	SA(A)A
    82	ACH2QH
    83 	QAAAA@1
    84	NH2
    85	CN(C)C
    86	CH2QCH2
    87	X!A$A
    88	S
    89	OAAAO
    90	QHAACH2A
    91	QHAAACH2A
    92	OC(N)C
    93	QCH3
    94	QN
    95	NAAO
    96	5M RING
    97	NAAAO
    98	QAAAAA@1
    99	C=C
    100	ACH2N
    101	8M RING
    102	QO
    103	CL
    104	QHACH2A
    105	A$A($A)$A
    106	QA(Q)Q
    107	XA(A)A
    108	CH3AAACH2A
    109	ACH2O
    110	NCO
    111	NACH2A
    112	AA(A)(A)A
    113	Onot%A%A
    114	CH3CH2A
    115	CH3ACH2A
    116	CH3AACH2A
    117	NAO
    118	ACH2CH2A > 1
    119	N=A
    120	HETEROCYCLIC ATOM > 1 (&...)
    121	N HETEROCYCLE
    122	AN(A)A
    123	OCO
    124	QQ
    125	AROMATIC RING > 1
    126	A!O!A
    127	A$A!O > 1 (&...)
    128	ACH2AAACH2A
    129	ACH2AACH2A
    130	QQ > 1 (&...)
    131	QH > 1
    132	OACH2A
    133	A$A!N
    134	X (HALOGEN)
    135	Nnot%A%A
    136	O=A > 1
    137	HETEROCYCLE
    138	QCH2A > 1 (&...)
    139	OH
    140	O > 3 (&...)
    141	CH3 > 2 (&...)
    142	N > 1
    143	A$A!O
    144	Anot%A%Anot%A
    145	6M RING > 1
    146	O > 2
    147	ACH2CH2A
    148	AQ(A)A
    149	CH3 > 1
    150	A!A$A!A
    151	NH
    152	OC(C)C
    153	QCH2A
    154	C=O
    155	A!CH2!A
    156	NA(A)A
    157	C-O
    158	C-N
    159	O > 1
    160	CH3
    161	N
    162	AROMATIC
    163	6M RING
    164	O
    165	RING
    166 	FRAGMENTS

MACCS 322 keys set as defined in tables 1, 2 and 3 [ Ref 46 ] include:

    o 26 atom properties of type P, as listed in Table 1
    o 32 one-atom environments, as listed in Table 3
    o 264 atom-bond-atom combinations listed in Table 4

Total number of keys in three tables is : 322

Atom symbol, X, used for 322 keys [ Ref 46 ] doesn't refer to Halogens as it does for 166 keys. In
order to keep the definition of 322 keys consistent with the published definitions, the symbol X is
used to imply "others" atoms, but it's internally mapped to symbol X as defined for 166 keys
during the generation of key values.

Atom properties-based keys (26):

    Key   Description
    1     A(AAA) or AA(A)A - atom with at least three neighbors
    2     Q - heteroatom
    3     Anot%not-A - atom involved in one or more multiple bonds, not aromatic
    4     A(AAAA) or AA(A)(A)A - atom with at least four neighbors
    5     A(QQ) or QA(Q) - atom with at least two heteroatom neighbors
    6     A(QQQ) or QA(Q)Q - atom with at least three heteroatom neighbors
    7     QH - heteroatom with at least one hydrogen attached
    8     CH2(AA) or ACH2A - carbon with at least two single bonds and at least
          two hydrogens attached
    9     CH3(A) or ACH3 - carbon with at least one single bond and at least three
          hydrogens attached
    10    Halogen
    11    A(-A-A-A) or A-A(-A)-A - atom has at least three single bonds
    12    AAAAAA@1 > 2 - atom is in at least two different six-membered rings
    13    A($A$A$A) or A$A($A)$A - atom has more than two ring bonds
    14    A$A!A$A - atom is at a ring/chain boundary. When a comparison is done
          with another atom the path passes through the chain bond.
    15    Anot%A%Anot%A - atom is at an aromatic/nonaromatic boundary. When a
          comparison is done with another atom the path
          passes through the aromatic bond.
    16    A!A!A  - atom with more than one chain bond
    17    A!A$A!A - atom is at a ring/chain boundary. When a comparison is done
          with another atom the path passes through the ring bond.
    18    A%Anot%A%A - atom is at an aromatic/nonaromatic boundary. When a
          comparison is done with another atom the
          path passes through the nonaromatic bond.
    19    HETEROCYCLE - atom is a heteroatom in a ring.
    20    rare properties: atom with five or more neighbors, atom in
          four or more rings, or atom types other than
          H, C, N, O, S, F, Cl, Br, or I
    21    rare properties: atom has a charge, is an isotope, has two or
          more multiple bonds, or has a triple bond.
    22    N - nitrogen
    23    S - sulfur
    24    O - oxygen
    25    A(AA)A(A)A(AA) - atom has two neighbors, each with three or
          more neighbors (including the central atom).
    26    CHACH2 - atom has two hydrocarbon (CH2) neighbors

Atomic environments properties-based keys (32):

    Key   Description
    27    C(CC)
    28    C(CCC)
    29    C(CN)
    30    C(CCN)
    31    C(NN)
    32    C(NNC)
    33    C(NNN)
    34    C(CO)
    35    C(CCO)
    36    C(NO)
    37    C(NCO)
    38    C(NNO)
    39    C(OO)
    40    C(COO)
    41    C(NOO)
    42    C(OOO)
    43    Q(CC)
    44    Q(CCC)
    45    Q(CN)
    46    Q(CCN)
    47    Q(NN)
    48    Q(CNN)
    49    Q(NNN)
    50    Q(CO)
    51    Q(CCO)
    52    Q(NO)
    53    Q(CNO)
    54    Q(NNO)
    55    Q(OO)
    56    Q(COO)
    57    Q(NOO)
    58    Q(OOO)

Note: The first symbol is the central atom, with atoms bonded to the central atom listed in
parentheses. Q is any non-C, non-H atom. If only two atoms are in parentheses, there is
no implication concerning the other atoms bonded to the central atom.

Atom-Bond-Atom properties-based keys: (264)

    Key   Description
    59    C-C
    60    C-N
    61    C-O
    62    C-S
    63    C-Cl
    64    C-P
    65    C-F
    66    C-Br
    67    C-Si
    68    C-I
    69    C-X
    70    N-N
    71    N-O
    72    N-S
    73    N-Cl
    74    N-P
    75    N-F
    76    N-Br
    77    N-Si
    78    N-I
    79    N-X
    80    O-O
    81    O-S
    82    O-Cl
    83    O-P
    84    O-F
    85    O-Br
    86    O-Si
    87    O-I
    88    O-X
    89    S-S
    90    S-Cl
    91    S-P
    92    S-F
    93    S-Br
    94    S-Si
    95    S-I
    96    S-X
    97    Cl-Cl
    98    Cl-P
    99    Cl-F
    100   Cl-Br
    101   Cl-Si
    102   Cl-I
    103   Cl-X
    104   P-P
    105   P-F
    106   P-Br
    107   P-Si
    108   P-I
    109   P-X
    110   F-F
    111   F-Br
    112   F-Si
    113   F-I
    114   F-X
    115   Br-Br
    116   Br-Si
    117   Br-I
    118   Br-X
    119   Si-Si
    120   Si-I
    121   Si-X
    122   I-I
    123   I-X
    124   X-X
    125   C=C
    126   C=N
    127   C=O
    128   C=S
    129   C=Cl
    130   C=P
    131   C=F
    132   C=Br
    133   C=Si
    134   C=I
    135   C=X
    136   N=N
    137   N=O
    138   N=S
    139   N=Cl
    140   N=P
    141   N=F
    142   N=Br
    143   N=Si
    144   N=I
    145   N=X
    146   O=O
    147   O=S
    148   O=Cl
    149   O=P
    150   O=F
    151   O=Br
    152   O=Si
    153   O=I
    154   O=X
    155   S=S
    156   S=Cl
    157   S=P
    158   S=F
    159   S=Br
    160   S=Si
    161   S=I
    162   S=X
    163   Cl=Cl
    164   Cl=P
    165   Cl=F
    166   Cl=Br
    167   Cl=Si
    168   Cl=I
    169   Cl=X
    170   P=P
    171   P=F
    172   P=Br
    173   P=Si
    174   P=I
    175   P=X
    176   F=F
    177   F=Br
    178   F=Si
    179   F=I
    180   F=X
    181   Br=Br
    182   Br=Si
    183   Br=I
    184   Br=X
    185   Si=Si
    186   Si=I
    187   Si=X
    188   I=I
    189   I=X
    190   X=X
    191   C#C
    192   C#N
    193   C#O
    194   C#S
    195   C#Cl
    196   C#P
    197   C#F
    198   C#Br
    199   C#Si
    200   C#I
    201   C#X
    202   N#N
    203   N#O
    204   N#S
    205   N#Cl
    206   N#P
    207   N#F
    208   N#Br
    209   N#Si
    210   N#I
    211   N#X
    212   O#O
    213   O#S
    214   O#Cl
    215   O#P
    216   O#F
    217   O#Br
    218   O#Si
    219   O#I
    220   O#X
    221   S#S
    222   S#Cl
    223   S#P
    224   S#F
    225   S#Br
    226   S#Si
    227   S#I
    228   S#X
    229   Cl#Cl
    230   Cl#P
    231   Cl#F
    232   Cl#Br
    233   Cl#Si
    234   Cl#I
    235   Cl#X
    236   P#P
    237   P#F
    238   P#Br
    239   P#Si
    240   P#I
    241   P#X
    242   F#F
    243   F#Br
    244   F#Si
    245   F#I
    246   F#X
    247   Br#Br
    248   Br#Si
    249   Br#I
    250   Br#X
    251   Si#Si
    252   Si#I
    253   Si#X
    254   I#I
    255   I#X
    256   X#X
    257   C$C
    258   C$N
    259   C$O
    260   C$S
    261   C$Cl
    262   C$P
    263   C$F
    264   C$Br
    265   C$Si
    266   C$I
    267   C$X
    268   N$N
    269   N$O
    270   N$S
    271   N$Cl
    272   N$P
    273   N$F
    274   N$Br
    275   N$Si
    276   N$I
    277   N$X
    278   O$O
    279   O$S
    280   O$Cl
    281   O$P
    282   O$F
    283   O$Br
    284   O$Si
    285   O$I
    286   O$X
    287   S$S
    288   S$Cl
    289   S$P
    290   S$F
    291   S$Br
    292   S$Si
    293   S$I
    294   S$X
    295   Cl$Cl
    296   Cl$P
    297   Cl$F
    298   Cl$Br
    299   Cl$Si
    300   Cl$I
    301   Cl$X
    302   P$P
    303   P$F
    304   P$Br
    305   P$Si
    306   P$I
    307   P$X
    308   F$F
    309   F$Br
    310   F$Si
    311   F$I
    312   F$X
    313   Br$Br
    314   Br$Si
    315   Br$I
    316   Br$X
    317   Si$Si
    318   Si$I
    319   Si$X
    320   I$I
    321   I$X
    322   X$X

=item B<SetSize>

    $MACCSKeys->SetSize($Size);

Sets size of MACCS keys and returns I<MACCSKeys>. Possible values: I<166 or 322>.

=item B<SetType>

    $MACCSKeys->SetType($Type);

Sets type of MACCS keys and returns I<MACCSKeys>. Possible values: I<MACCSKeysBits or
MACCSKeysCount>.

=item B<StringifyMACCSKeys>

    $String = $MACCSKeys->StringifyMACCSKeys();

Returns a string containing information about I<MACCSKeys> object.

=back

=head1 AUTHOR

Manish Sud <msud@san.rr.com>

=head1 SEE ALSO

Fingerprints.pm, FingerprintsStringUtil.pm, AtomNeighborhoodsFingerprints.pm,
AtomTypesFingerprints.pm, EStateIndiciesFingerprints.pm, ExtendedConnectivityFingerprints.pm,
PathLengthFingerprints.pm, TopologicalAtomPairsFingerprints.pm, TopologicalAtomTripletsFingerprints.pm,
TopologicalAtomTorsionsFingerprints.pm, TopologicalPharmacophoreAtomPairsFingerprints.pm,
TopologicalPharmacophoreAtomTripletsFingerprints.pm

=head1 COPYRIGHT

Copyright (C) 2018 Manish Sud. All rights reserved.

This file is part of MayaChemTools.

MayaChemTools is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.

=cut
