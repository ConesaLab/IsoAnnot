// Signal : HSL3|U0001|Histone 3'UTR stem-loop structure (HSL3)|3|Pattern Similarity|*2759*
r1={au,ua,gc,cg,gu,ug}
0...1 mmmm p1=ggyyy u hhuh a r1~p1 mm 0...3
##
// Signal : IRE|U0002|Iron Responsive Element (IRE)|all|Pattern Similarity|*33208*
r1={au,ua,gc,cg,gu,ug}
(p1=2...8 c p2=5...5 CAGWGH r1~p2 r1~p1 | 
p3=2...8 nnc p4=5...5 CAGWGH r1~p4 n r1~p3 |
p5=6...8 c p6=2...2 p7=3...3 CAGWGH r1~p7 u r1~p6 r1~p5 )
##
// Signal : SECIS1|U0003|Selenocysteine Insertion Sequence - type 1 (SECIS1)|3|Pattern Similarity|*33208*
(R1={AU,UA,GC,CG,GU,UG} 
P1=4...19 2...9  R UGAN   
P2=8...12 0...3 AAV 7...10 R1~P2[1,0,0] 
P2:(((^RR | ^MC) | ^SU) 6...10 ) 
NGAN 2...9 R1~P1
|
P3=4...19 P10=2...9 A UGAN
P4=8...12  P4:(((^RR | ^MC) | ^SU) 6...10 )
0...3 AAV P5=11...14 P5:(0...1 P6=2...7 3...6 R1~P6 0...3$)
(((R1~P4[1,0,0] | R1~P4[1,1,0]) | R1~P4[1,0,1]) | R1~P4[0,1,1]) 
NGAN 2...9 R1~P3
)
##
// Signal : SECIS2|U0004|Selenocysteine Insertion Sequence - type 2 (SECIS2)|3|Pattern Similarity|*33208*
r1={au,ua,gc,cg,gu,ug}
p1=4...19
p2=2...9 a
ugan
p3=8...12  p3:(((^rr | ^mc) | ^su) 6...10 )
p4=0...3 p5=aav p6=11...14 p6:(0...1 p7=2...7 3...6 r1~p7 0...3$)
(((r1~p3[1,0,0] | r1~p3[1,1,0]) | r1~p3[1,0,1]) | r1~p3[0,1,1])
ngan
p10=2...9
r1~p1
##
// Signal : APP_SCE|U0005|Amyloid precursor protein mRNA stability control element (APP_SCE)|all|Pattern Similarity|*7742*
UCUCUUUACAUUUUGGUCUCUACACUACA[3,1,0]
##
// Signal : CPE|U0006|Cytoplasmic polyadenylation element (CPE)|3|Pattern Similarity|*33208*
uuuuuwau 1...200 aauaaa 1...100 $
##
// Signal : TGE|U0007|TGE translational regulation element (TGE)|3|Pattern Similarity|*33208*
suca 5...14 mm uuuc yy 4...9 uwucu
##
// Signal : NANOS_TCE|U0008|Nanos translation control element (NANOS_TCE)|3|Pattern Similarity|*33208*
r1={au,ua,cg,gc,gu,ug}
p1=4...10 p2=6...9 cuggc r1~p2 0...3
p3=10...12 7...8 r1~p3[1,0,0] 0...2
r1~p1 p3:(p4=uauaua)
##
// Signal : 15-LOX-DICE|U0009|15-Lipoxygenase Differentiation Control Element (15-LOX-DICE)|3|Pattern Similarity|*7742*
20 > repeat (cccyryccucy[1,1,0] p1=1...7 ang) 0...0 > 1
##
// Signal : ARE2|U0010|AU-rich class-2 Element (ARE2)|3|Pattern Similarity|*33208*
wwwwwwwwwwwwwwwwwwww[8,0,0]
p1=auuu p2=0...4 p1 p3=0...4 p1 p4=0...4 p1 
p1/p2/p1/p3/p1/p4/p1:(p5=auuuauuua)
##
// Signal : TOP|U0011|Terminal Oligopyrimidine Tract (TOP)|5|Pattern Similarity|*2759*
^C yyyyyyyyyyyyyy[0,11,0]G
##
// Signal : GLUT1|U0012|Glusose transporter type-1 3'UTR cis-acting element (GLUT1)|3|Pattern Similarity|*7742*
ccaaccactc
##
// Signal : TNF|U0013|Tumor necrosis factor alpha 3'UTR cis-acting element (TNF)|3|Pattern Similarity|*7742*
uuauguuu 0...8 ar auauuuau surayynrru ugucw
##
// Signal : VIM3|U0014|Vimentin 3'UTR cis-acting element (VIM3)|3|Pattern Similarity|*7742*
p1=yttrrrrraa[2,0,0] 0...4 
p2=cagctttcaagtgccttt[2,0,0] 0...2
p3=tscagtt[2,0,0] 6...7 
p4=gagcg[2,0,0] 0...1 p5=aagatw[2,0,0]
p1/p2/p3/p4/p5:(p6=yttrrrrraacagctttcaagtgccttttscagttgagcgaagatw[2,0,0])
##
// Signal : IRES|U0015|Internal Ribosome Entry Site (IRES)|5|Pattern Similarity|*2759*
r1={au,ua,cg,gc,gu,ug}
p1=5...6 0...6 p2=5...6 p3=0...2
p4=5...8 p5=3...8
r1~p4[1,0,0] p6=5...8 p7=3...5
r1~p6[1,0,0] p8=0...5
r1~p2[1,0,0] 0...6
r1~p1[1,0,0] 2...5 p9=5...6 p10=3...8
r1~p9[1,0,0] 3...10 $
##
// Signal : SXL_BS|U0016|SXL binding site (SXL_BS)|all|Pattern Similarity|*6656*
((((
UUUUUUUUUUUUUUUU[0,5,0]|
UUUUUUUUU 14...14 UUUUUUU 14...14 UUUUUUU 9...9 UUUUUUU)|
UCUUUUUGUUGUUUUUUUUCUAG)|
UAUUUUUUUUCACAG)|
UUUUUUGUUKUKUUUKUU)
##
// Signal : UNR_BS|U0017|UNR binding site (UNR_BS)|all|Pattern Similarity|*33208*
(RRRRRAAGUAR|RRRRRRRRAACRRR)
##
// Signal : RPMS12_TCE|U0018|Ribosomal S12 mitochondrial protein 5'UTR translation control element (RPMS12_TCE)|5|Pattern Similarity|*7742*
CCGCGACCTCACCTTTAGGTCCTGTG[3,0,0]
##
// Signal : BRE|U0019|Bruno 3'UTR responsive element (BRE)|3|Pattern Similarity|*6656*
p1=wuurururuu[1,0,0]  p1:(p2=wnurururu) p1:(nngngng[2,0,0]) 0...3
p3=wuurururuu[1,0,0]  p3:(p4=wnurururu) p3:(nngngng[2,0,0])
##
// Signal : ADH_DRE|U0020|Alcohol dehydrogenase 3'UTR downregulation control element (ADH_DRE)|3|Pattern Similarity|*6656*
AAGGCTGA
##
// Signal : BYDV_TE|U0021|Barley yellow dwarf virus translation control element (BYDV_TE)|3|Pattern Similarity|*33090*
r1={au,ua,gc,cg,gu,ug}
p1=7...7 3...3 p2=ggauccug[1,0,0] p2:(4...4 p3=4...4) gnrna ~p3
2...2 p4=4...4 4...4 ~p4 4...4 p5=nn nnn (p6=6...6 c (u | c) a[0,1,0]
guc a[0,1,0] (a | c) r1~p6 | p8=5...5 n p9=3...3 nn
p10=6...6  c (u | c) a[0,1,0] guc a[0,1,0] (a | c) r1~p10 n ~p9 nn ~p8)
nnn r1~p5 6...6 r1~p1[1,0,1]
##
// Signal : PRONEURAL-BOX|U0022|Proneural Box (PB)|3|Pattern Similarity|*33208*
AAUGGAAGACAAU
##
// Signal : K-BOX|U0023|K-Box (KB)|3|Pattern Similarity|*33208*
p1=cTGTGATa[1,0,0] p1:(p2=TGTGAT)
##
// Signal : BRD-BOX|U0024|Brd-Box (Brd)|3|Pattern Similarity|*33208*
AGCTTTA
##
// Signal : GY-BOX|U0025|GY-Box (GY)|3|Pattern Similarity|*33208*
GTCTTCC
##
// Signal : AR_CURE|U0026|Androgen receptor CU-rich element (AR_CURE)|3|Pattern Similarity|*33208*
p1=cuggg p2=10...40 p3=cccuccc
p2:(p4=yyyyyyyyyy)
p1/p2/p3:(c p5=uuuuuuuuuuuu[0,7,0] c)
##
// Signal : G3A|U0027|Elastin G3A 3'UTR stability motif (G3A)|3|Pattern Similarity|*7742*
tgnn p1=rrrrrrrrrrrrrrrrrrr[1,3,0]
p1:(ggga ggga 0...3 ggga 0...1 ggg)
##
// Signal : INS_SCE|U0028|Insulin 3'UTR stability element (INS_SCE)|3|Pattern Similarity|*33208*
p1=3...7 p1:(p3=RCC) p2=YCCYRCCCNCC[3,0,0] p2:(p4=NNCNRNNCNNN) 0...3 ((CC|AC)|CT)
0...3 TGSM 1...13 AATAAA
##
// Signal : ACTIN_ZIP3|U0029|Beta-actin 3'UTR zipcode (ACTIN_ZIP3)|3|Pattern Similarity|*33208*
gcggactgtta[2,0,0] 8...14 acaccctt[1,0,0] 3...5 tga 0...2 caaaaccta[1,0,0] 3...4 tgc
##
// Signal : GAP-43|U0030|Gap-43 Stabilization Element (GAP-43)|3|Pattern Similarity|*33208*
guccacuuuccu yyyyyy[0,6,0] muyucucucug
##
// Signal : CNDLE|U0031|CaMKII/Ng dendritic localization element (CNDLE)|3|Pattern Similarity|*33208*
yrggrnc 0...4 cs 6...10 cscagaga 0...1 tcc 0...1 ytcw 2...8 wgg
##
// Signal : AG-CRSD|U0032|alpha-globin 3'UTR C-rich stability determinant (AG-CRSD)|3|Pattern Similarity|*33208*
yyugcsns 0...1  ytkgsc 0...1 tcysnmmrnsccctyytycyctycyt 0...1 gcac
##
// Signal : uORF|U0033|Upstream Open Reading Frame (uORF)|5|Pattern Similarity|*2759*
aug p1=60...9999 ((uaa | uag) | uga) length(p1) mod 3
##
// Signal : GAIT|U0034|Gamma interferon activated inhibitor of translation (GAIT element)|3|Pattern Similarity|*33208*
AANGNUUACUUUGGAAUGRCUAUAMNCRUU[1,1,1]
##
// Signal : MBE|U0035|Musashi binding element (MBE)|3|Pattern Similarity|*33208*
R UUU[0,2,0] AGU
##
// Signal : HLE|U0036|Drosophila hairy mRNA localization element (HLE)|3|Pattern Similarity|*50557*
r1={au,ua,gc,cg,gu,ug}
p1=cauuuc acnn p2=cauu 0...1 cg ccrggaunn cgcn r1~p2 uugnuu 0...2 r1~p1
##
// Signal : MBP-A2RE11|U0037|Myelin Basic Protein Localization Element (MBP-A2RE11)|3|Pattern Similarity|*7742*
gswcacy GCCAwGGAGsC AsAsAGC
##
// Signal : Protamine-YRS|U0038|Protamine P1 3'UTR Y-Box recognition site (Protamine-YRS)|3|Pattern Similarity|*33208*
ca 0...4 wccaycnna 0...1 ctcctgcntga
##
// Signal : G-CSF_SLDE|U0039|Granulocyte colony-stimulating factor stem-loop destabilizing element (G-CSF_SLDE)|3|Pattern Similarity|*33208*
r1={au,ua,gc,cg}
ansym p1=tgtttaa tat r1~p1 smag
##
// Signal : Ren_SRE|U0040|Renin stability regulatory element (Ren_SRE)|3|Pattern Similarity|*7742*
uk 0...1  p1=cuuu  rrggrm 0...3 aaa 0...4 auaaagayuncauguu
##
// Signal : PTH1|U0041|PTH 3'UTR proximal cis-acting instability element (PTH1)|3|Pattern Similarity|*7742*
TGYTVTASACAGBRTARGGCAA
##
// Signal : PTH2|U0042|PTH 3'UTR distal cis-acting instability element (PTH2)|3|Pattern Similarity|*7742*
RTWKNYTMTTCTTTTTAAAGTA
##
// Signal : PAS|U0043|Polyadenylation Signal (PAS)|3|Pattern Similarity|*2759*
((((((((((((AAUAAA|(AUUAAA|AGUAAA))|UAUAAA)|CAUAAA)|GAUAAA)|AAUAUA)|AAUACA)|AAUACA)|AAUAGA)|AAUGAA)|ACUAAA)|AACAAA)|UUUAAA) p2=0...34 $
##
// Signal : PABP_ARS|U0044|PABP mRNA autoregulatoy repression sequence (PABP_ARS)|5|Pattern Similarity|*7742*
CCCC[0,1,0] AAAAAAAAAA[0,3,0] UUUAC AAAAAA[0,1,0]
##
// Signal : ApoB_CRE|U0045|ApoB 5'UTR cis-acting regulatory element (ApoB_CRE)|5|Pattern Similarity|*7742*
yragtnccynyctsrgyk 0...1 chgyckcnganganccngymsagcmagcmrs 0...1 gssgckaggccsrkg 0...6 ccaggccrcagccmrgragccrcccyacc 0...2 kcmryysgc
##
// Signal : TPP_riboswitch|U0046|Thiamin pyrophosphate riboswitch (TPP_riboswitch)|3|Pattern Similarity|*33090*
r1={au,ua,gc,cg,gu,ug}
p1=gg g p2=u p3=gyyug 5...50 r1~p3 ugagaa[1,0,0] r1~p2 gu r1~p1 cuuugaa[2,0,0]  p4=ccu gan p5=cagg[1,0,0] ruaaur[2,0,0] r1~p5 cgn r1~p4 ga gugugc[1,0,0]
##
// Signal : GRE|U0047|GU-rich destabilization element (GRE)|3|Pattern Similarity|*2759*
UGUUUGUUUGU[1,0,0]
##
// Signal : PG4|U0048|G-quadruplex structure (PG4)|5|Pattern Similarity|*2759*
3<repeat(ggg) 1...7
##
// Signal : POLB|U0049|DNA polymerase beta stem loopII regulatory element (POLB)|3|Pattern Similarity|*40674*
R1={AU,UA,CG,GC,UG,GU}
UUAUUK[0,1,0] P1=CY UA[0,0,1] P2=ASCUUUGC UVY R1~P2[1,0,0] CUUU[1,2,0] R1~P1[1,0,0] UGKUYU
##
// Signal : CAESAR|U0050|Cis-Acting Element of Structure-Anchored Repression (CAESAR)|3|Pattern Similarity|*7742*
r1={AU,UA,CG,GC,UG,GU}
W p1=YGGAGA[1,0,0] YA p2=UGGC VURRA[0,1,0] 0...2 r1~p2[1,0,0] (C|0...0) GAGA[2,4,0]
p3=GUGA[3,2,2] p4=GAG[2,1,0] 5...30 r1~p4[2,1,0] p5=AYYAG[2,2,1] 0...6 AACUUGAA[0,1,0] 0...6
r1~p5[3,1,2] r1~p3[3,2,2] AUCUC AU[1,0,3] r1~p1[1,3,0]
##
// Signal : U1A_PIE|U0051|U1A polyadenylation inhibition element (U1A_PIE)|3|Pattern Similarity|*7742*
R1={AU,UA,GC,CG,GU,UG}
P1=VAGC AUUGYAC P2=CCAG D P3=GUCU 3...8 Y R1~P3 AYUGHAC R1~P2 H R1~P1[1,0,0]
##
// Signal : bGlob_pre|U0052|Beta-globin pyrimidine rich element (bGlob_pre)|3|Pattern Similarity|*40674*
r1={AU,UA,GC,CG,GU,UG}
p1=DWD p2=NHN HHDNBN 0...1 p3=YYCYH WDRDNVNDWHNNYYVVNY 2...4  r1~p3[3,0,0] WAHDV 0...2 r1~p2[3,0,0] BNNH 0...2 r1~p1[2,0,0]
##
// Signal : DMRT1_RE|U0053|DMRT1 regulatory element (DMRT1_RE)|3|Pattern Similarity|*33208*
p1=CUGCUGCMGRU[1,0,0]
p1:(YUGNUGYNRRU)
##
// Signal : NRE|U0054|Nanos response element (NRE)|3|Pattern Similarity|*6656*
Y (GUUGU[1,0,0] | GUUGU[0,0,1]) 0...22 RAUUGUAHMHR (0...35 AUUGUA | 0...0)
##
// Signal : tropC_PLS|U0055|troponin C perinuclear localization signal (tropC_PLS)|3|Pattern Similarity|*40674*
VGHRCARRGRR[1,4,0] YYRRGGR[0,3,0] YYYYYY[0,2,0] RRSR[0,2,0] YYYBKGARR[0,3,0] YNYGSNYRTRTY[0,1,0]
##
// Signal : TLS|U0056|K10 Orb transport/localisation signal (TLS)|3|Pattern Similarity|*50557*
r1={AU,UA,CG,GC,GU,UG}
p1=VWWRWUY[0,1,0] p2=SWM p3=UUUUUAR 5...9 r1~p3 (Y|0...0) r1~p2[1,0,0] 1...2 r1~p1
##
// Signal : WLE3|U0057|Wingless localisation element 3 (WLE3)|3|Pattern Similarity|*50557*
r1={AU,UA,CG,GC,GU,UG}
p1=URYUY p2=KCAY 5...6 p3=UUKGS 6...7 r1~p3 4...4 r1~p2[1,0,0] 1...2 r1~p1[2,0,0]
##
// Signal : FLE|U0058|fushi tarazu localization element (FLE)|all|Pattern Similarity|*50557*
r1={AU,UA,CG,GC,GU,UG}
p1=UCUCG RC p2=CAUUUGY HGGUGACAC r1~p2 5...9 r1~p1 30...80 
p3=HACUC p4=CCU ACG p5=GUGCYYAK 3...6 r1~p5[3,2,0] 0...5 r1~p4[1,0,0] 2...4 r1~p3[2,0,0]
##
// Signal : MYC_LE|U0059|MYC mRNA localisation element (MYC_LE)|3|Pattern Similarity|*7742*
r1={AU,UA,CG,GC,GU,UG}
p1=BYTWNR (CA|AA|GA|AG|CAA) p2=NNTTT[1,0,0] GTA p3=TTTAA 6...16 r1~p3[2,0,0] 0...8 r1~p2[2,2,0] 0...1 r1~p1[3,1,0]
