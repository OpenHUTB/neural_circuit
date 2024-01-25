/*  Neuron parameters for makcel.cc simulation script */

/* nval.h */

/* To add parameters, edit "maknval.cc", compile and run "maknval > nval.n",
   then copy nval.n to "nval.h". Remove the param defs from the end of nval.n.
   Remove the nval.n table at the beginning of nval.h.
   Copy nval.h to "nval_var.h", "nval_var.cc", and "nval_var_set.cc"
   and edit these files, then remove this content from nval.h.
   Last, "make clean" and "make retsim".
*/

#define XCONE        0	/* cones */
#define XROD         1	/* rods */
#define HBAT         2	/* hbat */
#define HA           3	/* Type A horizontal cells */
#define HB           4	/* Type B horizontal cells */
#define RBP          5	/* Rod bipolar cells */
#define DBP1         6	/* Depolarizing cone bipolar cell, type 1 */
#define DBP2         7	/* Depolarizing cone bipolar cell, type 2 */
#define DBP3         8	/* Depolarizing cone bipolar cell, type 3 */
#define DBP4         9	/* Depolarizing cone bipolar cell, type 4 */
#define HBP1        10	/* Hyperpolarizing bipolar cell, type 1 */
#define HBP2        11	/* Hyperpolarizing bipolar cell, type 2 */
#define A17         12	/* A17 amacrine cells, feedback to RBP */
#define AII         13	/* AII amacrine cells */
#define SBAC        14	/* Starburst amacrine cells */
#define AM          15	/* Amacrine cell type 1 */
#define AM2         16	/* Amacrine cell type 2 */
#define AM3         17	/* Amacrine cell type 3 */
#define AM4         18	/* Amacrine cell type 4 */
#define AMH         19	/* Amacrine cells, hyperpolarizing */
#define AMH2        20	/* Amacrine cells, hyperpolarizing */
#define AMS         21	/* Amacrine cells, small-field */
#define AMHS        22	/* Amacrine cells, small-field hyperpol */
#define GCA         23	/* Ganglion cells, On-type, alpha */
#define GCB         24	/* Ganglion cells, On-type, beta */
#define DSGC        25	/* Direction-selective ganglion cells */
#define GCAOFF      26	/* Ganglion cells, Off-type, alpha */
#define GCBOFF      27	/* Ganglion cells, Off-type, beta */
#define NCELTYPES   28	/* Number of cell types */

#define MAKE         0	 /* whether to make this cell type */
#define MAKE_DEND    1	 /* whether to make dendrites */
#define MAKE_AXON    2	 /* whether to make axon */
#define MAKE_DIST    3	 /* whether to make axon distal */
#define NMADE        4	 /* number of cells made */
#define MAXNUM       5	 /* maximum number of cells of this type */
#define NCOLOR       6	 /* color of this cell type for display */
#define NCMAP        7	 /* color map for color of this cell type */
#define MAXCOV       8	 /* max coverage factor (for arrays) */
#define MAXSYNI      9	 /* max number of syn input cells */
#define MAXSYNO     10	 /* max number of syn output cells */
#define DENS        11	 /* density of this type (per mm2) */
#define REGU        12	 /* regularity (mean/stdev) of spacing */
#define MORPH       13	 /* morphology (=0 -> file, or artificial) */
#define COMPLAM     14	 /* compartment size (default=complam) */
#define BIOPHYS     15	 /* add biophys properties (chan dens file) */
#define CHNOISE     16	 /* add membrane channel noise properties   */
#define RATIOK      17	 /* set K density values as ratio from Na */
#define VSTART      18	 /* initial resting potential */
#define VREV        19	 /* membrane potential for Rm (VCl) */
#define NRM         20	 /* the cell's Rm, 0 => use default (drm) */
#define NRI         21	 /* the cell's Ri, 0 => use default (dri) */
#define SOMADIA     22	 /* Soma diameter */
#define SOMAZ       23	 /* Z location (x,y loc determ. by array) */
#define SOMAZ2      24	 /* Z location for 2nd morph (x,y loc determ. by array) */
#define DENDARB     25	 /* type of dendritic tree */
#define DENDARBZ    26	 /* dendritic arborization level */
#define DENZDIST    27	 /* dendritic arborization z tolerance */
#define STRATDIA    28	 /* stratif. annulus dia (fract of treedia) */
#define DTIPDIA     29	 /* diameter of dendritic tips */
#define DTIPLEN     30	 /* length of dendritic tips */
#define DTREEDIA    31	 /* diameter of dendritic tree */
#define ARBSCALE    32	 /* scale for dia of real morph dend tree */
#define ARBSCALE2   33	 /* scale for dia of 2nd real morph dend tree */
#define DENDDIA     34	 /* dend dia scale for real morph */
#define DENDDIA2    35	 /* dend dia scale for 2nd real morph */
#define AXARBT      36	 /* type of axonal tree */
#define AXARBZ      37	 /* axonal arborization level */
#define AXTIPDIA    38	 /* diameter of axonal tips */
#define AXARBDIA    39	 /* diameter of axonal arbor */
#define MAXSDIST    40	 /* maximum synaptic distance */
#define TAPERSPC    41	 /* space constant of diameter taper */
#define TAPERABS    42	 /* abs diameter for taper */
#define NDENDR      43	 /* number of first-order dendrites */
#define GROWTHR     44	 /* distance thresh for growth of dendrites */
#define SEGLEN      45	 /* length of dendrite segments */

#define CELPRE1     46	 /* cell type to connect to (neg, no conn) */
#define CONPRE1     47	 /* connection number of presyn cell */
#define CELCONV1    48	 /* number of presyn cells to connect to */
#define GROWPOST1   49	 /* grow when making conn from presyn cell */
#define CELPRE2     50	 /* cell type to connect to (neg, no conn) */
#define CONPRE2     51	 /* connection number of presyn cell */
#define CELCONV2    52	 /* number of presyn cells to connect to */
#define GROWPOST2   53	 /* grow when making conn from presyn cell */
#define CELPRE3     54	 /* cell type to connect to (neg, no conn) */
#define CONPRE3     55	 /* connection number of presyn cell */
#define CELCONV3    56	 /* number of presyn cells to connect to */
#define GROWPOST3   57	 /* grow when making conn from presyn cell */
#define CELPRE4     58	 /* cell type to connect to (neg, no conn) */
#define CONPRE4     59	 /* connection number of presyn cell */
#define CELCONV4    60	 /* number of presyn cells to connect to */
#define GROWPOST4   61	 /* grow when making conn from presyn cell */
#define CELPRE5     62	 /* cell type to connect to (neg, no conn) */
#define CONPRE5     63	 /* connection number of presyn cell */
#define CELCONV5    64	 /* number of presyn cells to connect to */
#define GROWPOST5   65	 /* grow when making conn from presyn cell */
#define CELPRE6     66	 /* cell type to connect to (neg, no conn) */
#define CONPRE6     67	 /* connection number of presyn cell */
#define CELCONV6    68	 /* number of presyn cells to connect to */
#define GROWPOST6   69	 /* grow when making conn from presyn cell */
#define CELPRE7     70	 /* cell type to connect to (neg, no conn) */
#define CONPRE7     71	 /* connection number of presyn cell */
#define CELCONV7    72	 /* number of presyn cells to connect to */
#define GROWPOST7   73	 /* grow when making conn from presyn cell */
#define CELPRE8     74	 /* cell type to connect to (neg, no conn) */
#define CONPRE8     75	 /* connection number of presyn cell */
#define CELCONV8    76	 /* number of presyn cells to connect to */
#define GROWPOST8   77	 /* grow when making conn from presyn cell */
#define CELPRE9     78	 /* cell type to connect to (neg, no conn) */
#define CONPRE9     79	 /* connection number of presyn cell */
#define CELCONV9    80	 /* number of presyn cells to connect to */
#define GROWPOST9   81	 /* grow when making conn from presyn cell */
#define CELPRE10    82	 /* cell type to connect to (neg, no conn) */
#define CONPRE10    83	 /* connection number of presyn cell */
#define CELCONV10   84	 /* number of presyn cells to connect to */
#define GROWPOST10  85	 /* grow when making conn from presyn cell */

#define CELPOST1    86	 /* cell type to connect to (neg, no conn) */
#define CONPOST1    87	 /* connection number for postsyn cell */
#define CELDIV1     88	 /* number of postsyn cells to connect to */
#define GROWPRE1    89	 /* grow when making conn to postsyn cell */
#define SYNREG1     90	 /* synaptic region in presyn dendritic tree */
#define SYNREGP1    91	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC1    92	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI1    93	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO1    94	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI1    95	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO1    96	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG1     97	 /* angle for presynaptic node rel to soma */
#define SYNRNG1     98	 /* range of presyn dendritic angles */
#define USEDYAD1    99	 /* synapse is dyad using preexisting type */
#define DYADTYP1   100	 /* type of dyad synapse to connect with */
#define AUTAPSE1   101	 /* synapse back to presynaptic node */
#define SYNNUM1    102	 /* number of synapses per connection */
#define SENSCA1    103	 /* synaptic release sensitivity calcium */
#define CAEGAIN1   104	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL1   105	 /* synaptic readily releasable pool */
#define SRRPOOLG1  106	 /* synaptic readily releasable pool gain */
#define SMRRPOOL1  107	 /* synaptic max readily releasable pool */
#define SMAXRATE1  108	 /* maximum sustained synaptic release rate */
#define SGAIN1     109	 /* synaptic gain */
#define SVGAIN1    110	 /* synaptic vgain */
#define SDURH1     111	 /* synaptic high pass time const. (ms) */
#define SNFILTH1   112	 /* synaptic high pass nfilt */
#define SHGAIN1    113	 /* synaptic high pass gain */
#define SHOFFS1    114	 /* synaptic high pass offset */
#define SVSIZ1     115	 /* synaptic vesicle size */
#define SCOND1     116	 /* synaptic conductance */
#define SCMUL1     117	 /* synaptic conductance mult for region */
#define SCGRAD1    118	 /* synaptic conductance gradient from soma */
#define SEGRAD1    119	 /* synaptic conductance expon grad fr soma */
#define STHRESH1   120	 /* synaptic threshold */
#define SVNOISE1   121	 /* 1 -> vesicle noise, override, vnoise=0 */
#define SCOV1      122	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR1      123	 /* synaptic event time const. (ms) */
#define SFALL1     124	 /* synaptic event fall time const. (ms) */
#define SNFILT1    125	 /* synaptic vesicle nfilt */
#define STRCONC1   126	 /* synaptic transmitter concentration. */
#define SRESP1     127	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA1      128	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX1   129	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM1     130	 /* synaptic postsyn Ca pump Km. */
#define SVFBG1     131	 /* synaptic voltage feedback gain */
#define SVFBO1     132	 /* synaptic voltage feedback offset */
#define SCNFILT1   133	 /* second mesng. nfilt */
#define SCDUR1     134	 /* second mesng. time const. (ms) */
#define SCGAIN1    135	 /* synaptic second messenger gain */
#define SCOFF1     136	 /* synaptic second messenger offset */
#define SCNOISE1   137	 /* 1 -> channel noise, override cnoise=0 */
#define SNCHAN1    138	 /* number of channels */
#define SUNIT1     139	 /* synaptic channel unitary conductace */
#define SVREV1     140	 /* synaptic reversal potential */

#define CELPOST2   141	 /* cell type to connect to (neg, no conn) */
#define CONPOST2   142	 /* connection number for postsyn cell */
#define CELDIV2    143	 /* number of postsyn cells to connect to */
#define GROWPRE2   144	 /* grow when making conn to postsyn cell */
#define SYNREG2    145	 /* synaptic region in presyn dendritic tree */
#define SYNREGP2   146	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC2   147	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI2   148	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO2   149	 /* outer rad of annulus in presyn dend tree */
#define SYNANPI2   150	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO2   151	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG2    152	 /* angle for presynaptic node rel to soma */
#define SYNRNG2    153	 /* range of presyn dendritic angles */
#define USEDYAD2   154	 /* synapse is dyad using preexisting type */
#define DYADTYP2   155	 /* type of dyad synapse to connect with */
#define AUTAPSE2   156	 /* synapse back to presynaptic node */
#define SYNNUM2    157	 /* number of synapses per connection */
#define SENSCA2    158	 /* synaptic release sensitivity calcium */
#define CAEGAIN2   159	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL2   160	 /* synaptic readily releasable pool */
#define SRRPOOLG2  161	 /* synaptic readily releasable pool gain */
#define SMRRPOOL2  162	 /* synaptic max readily releasable pool */
#define SMAXRATE2  163	 /* maximum sustained synaptic release rate */
#define SGAIN2     164	 /* synaptic gain */
#define SVGAIN2    165	 /* synaptic vgain */
#define SDURH2     166	 /* synaptic high pass time const. (ms) */
#define SNFILTH2   167	 /* synaptic high pass nfilt */
#define SHGAIN2    168	 /* synaptic high pass gain */
#define SHOFFS2    169	 /* synaptic high pass offset */
#define SVSIZ2     170	 /* synaptic vesicle size */
#define SCOND2     171	 /* synaptic conductance */
#define SCMUL2     172	 /* synaptic conductance mult for region */
#define SCGRAD2    173	 /* synaptic conductance grad from soma */
#define SEGRAD2    174	 /* synaptic conductance expon grad fr soma */
#define STHRESH2   175	 /* synaptic threshold */
#define SVNOISE2   176	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV2      177	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR2      178	 /* synaptic event time const. (ms) */
#define SFALL2     179	 /* synaptic event fall time const. (ms) */
#define SNFILT2    180	 /* synaptic vesicle nfilt */
#define STRCONC2   181	 /* synaptic transmitter concentration. */
#define SRESP2     182	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA2      183	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX2   184	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM2     185	 /* synaptic postsyn Ca pump Km. */
#define SVFBG2     186	 /* synaptic voltage feedback gain */
#define SVFBO2     187	 /* synaptic voltage feedback offset */
#define SCNFILT2   188	 /* second mesng. nfilt */
#define SCDUR2     189	 /* second mesng. time const. (ms) */
#define SCGAIN2    190	 /* synaptic second messenger gain */
#define SCOFF2     191	 /* synaptic second messenger offset */
#define SCNOISE2   192	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN2    193	 /* number of channels */
#define SUNIT2     194	 /* synaptic channel unitary conductace */
#define SVREV2     195	 /* synaptic reversal potential */

#define CELPOST3   196	 /* cell type to connect to (neg, no conn) */
#define CONPOST3   197	 /* connection number for postsyn cell */
#define CELDIV3    198	 /* number of postsyn cells to connect to */
#define GROWPRE3   199	 /* grow when making conn to postsyn cell */
#define SYNREG3    200	 /* synaptic region in presyn dendritic tree */
#define SYNREGP3   201	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC3   202	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI3   203	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO3   204	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI3   205	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO3   206	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG3    207	 /* angle for presynaptic node rel to soma */
#define SYNRNG3    208	 /* range of presyn dendritic angles */
#define USEDYAD3   209	 /* synapse is dyad using preexisting type */
#define DYADTYP3   210	 /* type of dyad synapse to connect with */
#define AUTAPSE3   211	 /* synapse back to presynaptic node */
#define SYNNUM3    212	 /* number of synapses per connection */
#define SENSCA3    213	 /* synaptic release sensitivity calcium */
#define CAEGAIN3   214	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL3   215	 /* synaptic readily releasable pool */
#define SRRPOOLG3  216	 /* synaptic readily releasable pool gain */
#define SMRRPOOL3  217	 /* synaptic max readily releasable pool */
#define SMAXRATE3  218	 /* maximum sustained synaptic release rate */
#define SGAIN3     219	 /* synaptic gain */
#define SVGAIN3    220	 /* synaptic vgain */
#define SDURH3     221	 /* synaptic high pass time const. (ms) */
#define SNFILTH3   222	 /* synaptic high pass nfilt */
#define SHGAIN3    223	 /* synaptic high pass gain */
#define SHOFFS3    224	 /* synaptic high pass offset */
#define SVSIZ3     225	 /* synaptic vesicle size */
#define SCOND3     226	 /* synaptic conductance */
#define SCMUL3     227	 /* synaptic conductance mult for region */
#define SCGRAD3    228	 /* synaptic conductance gradient from soma */
#define SEGRAD3    229	 /* synaptic conductance expon grad fr soma */
#define STHRESH3   230	 /* synaptic threshold */
#define SVNOISE3   231	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV3      232	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR3      233	 /* synaptic event time const. (ms) */
#define SFALL3     234	 /* synaptic event fall time const. (ms) */
#define SNFILT3    235	 /* synaptic vesicle nfilt */
#define STRCONC3   236	 /* synaptic transmitter concentration. */
#define SRESP3     237	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA3      238	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX3   239	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM3     240	 /* synaptic postsyn Ca pump Km. */
#define SVFBG3     241	 /* synaptic voltage feedback gain */
#define SVFBO3     242	 /* synaptic voltage feedback offset */
#define SCNFILT3   243	 /* second mesng. nfilt */
#define SCDUR3     244	 /* second mesng. time const. (ms) */
#define SCGAIN3    245	 /* synaptic second messenger gain */
#define SCOFF3     246	 /* synaptic second messenger offset */
#define SCNOISE3   247	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN3    248	 /* number of channels */
#define SUNIT3     249	 /* synaptic channel unitary conductace */
#define SVREV3     250	 /* synaptic reversal potential */

#define CELPOST4   251	 /* cell type to connect to (neg, no conn) */
#define CONPOST4   252	 /* connection number for postsyn cell */
#define CELDIV4    253	 /* number of postsyn cells to connect to */
#define GROWPRE4   254	 /* grow when making conn to postsyn cell */
#define SYNREG4    255	 /* synaptic region in presyn dendritic tree */
#define SYNREGP4   256	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC4   257	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI4   258	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO4   259	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI4   260	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO4   261	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG4    262	 /* angle for presynaptic node rel to soma */
#define SYNRNG4    263	 /* range of presyn dendritic angles */
#define USEDYAD4   264	 /* synapse is dyad using preexisting type */
#define DYADTYP4   265	 /* type of dyad synapse to connect with */
#define AUTAPSE4   266	 /* synapse back to presynaptic node */
#define SYNNUM4    267	 /* number of synapses per connection */
#define SENSCA4    268	 /* synaptic release sensitivity calcium */
#define CAEGAIN4   269	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL4   270	 /* synaptic readily releasable pool */
#define SRRPOOLG4  271	 /* synaptic readily releasable pool gain */
#define SMRRPOOL4  272	 /* synaptic max readily releasable pool */
#define SMAXRATE4  273	 /* maximum sustained synaptic release rate */
#define SGAIN4     274	 /* synaptic gain */
#define SVGAIN4    275	 /* synaptic vgain */
#define SDURH4     276	 /* synaptic high pass time const. (ms) */
#define SNFILTH4   277	 /* synaptic high pass nfilt */
#define SHGAIN4    278	 /* synaptic high pass gain */
#define SHOFFS4    279	 /* synaptic high pass offset */
#define SVSIZ4     280	 /* synaptic vesicle size */
#define SCOND4     281	 /* synaptic conductance */
#define SCMUL4     282	 /* synaptic conductance mult for region */
#define SCGRAD4    283	 /* synaptic conductance gradient from soma */
#define SEGRAD4    284	 /* synaptic conductance expon grad fr soma */
#define STHRESH4   285	 /* synaptic threshold */
#define SVNOISE4   286	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV4      287	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR4      288	 /* synaptic event time const. (ms) */
#define SFALL4     289	 /* synaptic event fall time const. (ms) */
#define SNFILT4    290	 /* synaptic vesicle nfilt */
#define STRCONC4   291	 /* synaptic transmitter concentration. */
#define SRESP4     292	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA4      293	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX4   294	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM4     295	 /* synaptic postsyn Ca pump Km. */
#define SVFBG4     296	 /* synaptic voltage feedback gain */
#define SVFBO4     297	 /* synaptic voltage feedback offset */
#define SCNFILT4   298	 /* second mesng. nfilt */
#define SCDUR4     299	 /* second mesng. time const. (ms) */
#define SCGAIN4    300	 /* synaptic second messenger gain */
#define SCOFF4     301	 /* synaptic second messenger offset */
#define SCNOISE4   302	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN4    303	 /* number of channels */
#define SUNIT4     304	 /* synaptic channel unitary conductace */
#define SVREV4     305	 /* synaptic reversal potential */

#define CELPOST5   306	 /* cell type to connect to (neg, no conn) */
#define CONPOST5   307	 /* connection number for postsyn cell */
#define CELDIV5    308	 /* number of postsyn cells to connect to */
#define GROWPRE5   309	 /* grow when making conn to postsyn cell */
#define SYNREG5    310	 /* synaptic region in presyn dendritic tree */
#define SYNREGP5   311	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC5   312	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI5   313	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO5   314	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI5   315	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO5   316	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG5    317	 /* angle for presynaptic node rel to soma */
#define SYNRNG5    318	 /* range of presyn dendritic angles */
#define USEDYAD5   319	 /* synapse is dyad using preexisting type */
#define DYADTYP5   320	 /* type of dyad synapse to connect with */
#define AUTAPSE5   321	 /* synapse back to presynaptic node */
#define SYNNUM5    322	 /* number of synapses per connection */
#define SENSCA5    323	 /* synaptic release sensitivity calcium */
#define CAEGAIN5   324	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL5   325	 /* synaptic readily releasable pool */
#define SRRPOOLG5  326	 /* synaptic readily releasable pool gain */
#define SMRRPOOL5  327	 /* synaptic max readily releasable pool */
#define SMAXRATE5  328	 /* maximum sustained synaptic release rate */
#define SGAIN5     329	 /* synaptic gain */
#define SVGAIN5    330	 /* synaptic vgain */
#define SDURH5     331	 /* synaptic high pass time const. (ms) */
#define SNFILTH5   332	 /* synaptic high pass nfilt */
#define SHGAIN5    333	 /* synaptic high pass gain */
#define SHOFFS5    334	 /* synaptic high pass offset */
#define SVSIZ5     335	 /* synaptic vesicle size */
#define SCOND5     336	 /* synaptic conductance */
#define SCMUL5     337	 /* synaptic conductance mult for region */
#define SCGRAD5    338	 /* synaptic conductance gradient from soma */
#define SEGRAD5    339	 /* synaptic conductance expon grad fr soma */
#define STHRESH5   340	 /* synaptic threshold */
#define SVNOISE5   341	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV5      342	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR5      343	 /* synaptic event time const. (ms) */
#define SFALL5     344	 /* synaptic event fall time const. (ms) */
#define SNFILT5    345	 /* synaptic vesicle nfilt */
#define STRCONC5   346	 /* synaptic transmitter concentration. */
#define SRESP5     347	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA5      348	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX5   349	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM5     350	 /* synaptic postsyn Ca pump Km. */
#define SVFBG5     351	 /* synaptic voltage feedback gain */
#define SVFBO5     352	 /* synaptic voltage feedback offset */
#define SCNFILT5   353	 /* second mesng. nfilt */
#define SCDUR5     354	 /* second mesng. time const. (ms) */
#define SCGAIN5    355	 /* synaptic second messenger gain */
#define SCOFF5     356	 /* synaptic second messenger offset */
#define SCNOISE5   357	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN5    358	 /* number of channels */
#define SUNIT5     359	 /* synaptic channel unitary conductace */
#define SVREV5     360	 /* synaptic reversal potential */

#define CELPOST6   361	 /* cell type to connect to (neg, no conn) */
#define CONPOST6   362	 /* connection number for postsyn cell */
#define CELDIV6    363	 /* number of postsyn cells to connect to */
#define GROWPRE6   364	 /* grow when making conn to postsyn cell */
#define SYNREG6    365	 /* synaptic region in presyn dendritic tree */
#define SYNREGP6   366	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC6   367	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI6   368	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO6   369	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI6   370	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO6   371	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG6    372	 /* angle for presynaptic node rel to soma */
#define SYNRNG6    373	 /* range of presyn dendritic angles */
#define USEDYAD6   374	 /* synapse is dyad using preexisting type */
#define DYADTYP6   375	 /* type of dyad synapse to connect with */
#define AUTAPSE6   376	 /* synapse back to presynaptic node */
#define SYNNUM6    377	 /* number of synapses per connection */
#define SENSCA6    378	 /* synaptic release sensitivity calcium */
#define CAEGAIN6   379	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL6   380	 /* synaptic readily releasable pool */
#define SRRPOOLG6  381	 /* synaptic readily releasable pool gain */
#define SMRRPOOL6  382	 /* synaptic max readily releasable pool */
#define SMAXRATE6  383	 /* maximum sustained synaptic release rate */
#define SGAIN6     384	 /* synaptic gain */
#define SVGAIN6    385	 /* synaptic vgain */
#define SDURH6     386	 /* synaptic high pass time const. (ms) */
#define SNFILTH6   387	 /* synaptic high pass nfilt */
#define SHGAIN6    388	 /* synaptic high pass gain */
#define SHOFFS6    389	 /* synaptic high pass offset */
#define SVSIZ6     390	 /* synaptic vesicle size */
#define SCOND6     391	 /* synaptic conductance */
#define SCMUL6     392	 /* synaptic conductance mult for region */
#define SCGRAD6    393	 /* synaptic conductance gradient from soma */
#define SEGRAD6    394	 /* synaptic conductance expon grad fr soma */
#define STHRESH6   395	 /* synaptic threshold */
#define SVNOISE6   396	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV6      397	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR6      398	 /* synaptic event time const. (ms) */
#define SFALL6     399	 /* synaptic event fall time const. (ms) */
#define SNFILT6    400	 /* synaptic vesicle nfilt */
#define STRCONC6   401	 /* synaptic transmitter concentration. */
#define SRESP6     402	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA6      403	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX6   404	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM6     405	 /* synaptic postsyn Ca pump Km. */
#define SVFBG6     406	 /* synaptic voltage feedback gain */
#define SVFBO6     407	 /* synaptic voltage feedback offset */
#define SCNFILT6   408	 /* second mesng. nfilt */
#define SCDUR6     409	 /* second mesng. time const. (ms) */
#define SCGAIN6    410	 /* synaptic second messenger gain */
#define SCOFF6     411	 /* synaptic second messenger offset */
#define SCNOISE6   412	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN6    413	 /* number of channels */
#define SUNIT6     414	 /* synaptic channel unitary conductace */
#define SVREV6     415	 /* synaptic reversal potential */

#define CELPOST7   416	 /* cell type to connect to (neg, no conn) */
#define CONPOST7   417	 /* connection number for postsyn cell */
#define CELDIV7    418	 /* number of postsyn cells to connect to */
#define GROWPRE7   419	 /* grow when making conn to postsyn cell */
#define SYNREG7    420	 /* synaptic region in presyn dendritic tree */
#define SYNREGP7   421	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC7   422	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI7   423	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO7   424	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI7   425	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO7   426	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG7    427	 /* angle for presynaptic node rel to soma */
#define SYNRNG7    428	 /* range of presyn dendritic angles */
#define USEDYAD7   429	 /* synapse is dyad using preexisting type */
#define DYADTYP7   430	 /* type of dyad synapse to connect with */
#define AUTAPSE7   431	 /* synapse back to presynaptic node */
#define SYNNUM7    432	 /* number of synapses per connection */
#define SENSCA7    433	 /* synaptic release sensitivity calcium */
#define CAEGAIN7   434	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL7   435	 /* synaptic readily releasable pool */
#define SRRPOOLG7  436	 /* synaptic readily releasable pool gain */
#define SMRRPOOL7  437	 /* synaptic max readily releasable pool */
#define SMAXRATE7  438	 /* maximum sustained synaptic release rate */
#define SGAIN7     439	 /* synaptic gain */
#define SVGAIN7    440	 /* synaptic vgain */
#define SDURH7     441	 /* synaptic high pass time const. (ms) */
#define SNFILTH7   442	 /* synaptic high pass nfilt */
#define SHGAIN7    443	 /* synaptic high pass gain */
#define SHOFFS7    444	 /* synaptic high pass offset */
#define SVSIZ7     445	 /* synaptic vesicle size */
#define SCOND7     446	 /* synaptic conductance */
#define SCMUL7     447	 /* synaptic conductance mult for region */
#define SCGRAD7    448	 /* synaptic conductance gradient from soma */
#define SEGRAD7    449	 /* synaptic conductance expon grad fr soma */
#define STHRESH7   450	 /* synaptic threshold */
#define SVNOISE7   451	 /* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV7      452	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR7      453	 /* synaptic event time const. (ms) */
#define SFALL7     454	 /* synaptic event fall time const. (ms) */
#define SNFILT7    455	 /* synaptic vesicle nfilt */
#define STRCONC7   456	 /* synaptic transmitter concentration. */
#define SRESP7     457	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA7      458	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX7   459	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM7     460	 /* synaptic postsyn Ca pump Km. */
#define SVFBG7     461	 /* synaptic voltage feedback gain */
#define SVFBO7     462	 /* synaptic voltage feedback offset */
#define SCNFILT7   463	 /* second mesng. nfilt */
#define SCDUR7     464	 /* second mesng. time const. (ms) */
#define SCGAIN7    465	 /* synaptic second messenger gain */
#define SCOFF7     466	 /* synaptic second messenger offset */
#define SCNOISE7   467	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN7    468	 /* number of channels */
#define SUNIT7     469	 /* synaptic channel unitary conductace */
#define SVREV7     470	 /* synaptic reversal potential */

#define CELPOST8   471	 /* cell type to connect to (neg, no conn) */
#define CONPOST8   472	 /* connection number for postsyn cell */
#define CELDIV8    473	 /* number of postsyn cells to connect to */
#define GROWPRE8   474	 /* grow when making conn to postsyn cell */
#define SYNREG8    475	 /* synaptic region in presyn dendritic tree */
#define SYNREGP8   476	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC8   477	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI8   478	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO8   479	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI8   480	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO8   481	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG8    482	 /* angle for presynaptic node rel to soma */
#define SYNRNG8    483	 /* range of presyn dendritic angles */
#define USEDYAD8   484	 /* synapse is dyad using preexisting type */
#define DYADTYP8   485	 /* type of dyad synapse to connect with */
#define AUTAPSE8   486	 /* synapse back to presynaptic node */
#define SYNNUM8    487	 /* number of synapses per connection */
#define SENSCA8    488	 /* synaptic release sensitivity calcium */
#define CAEGAIN8   489	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL8   490	 /* synaptic readily releasable pool */
#define SRRPOOLG8  491	 /* synaptic readily releasable pool gain */
#define SMRRPOOL8  492	 /* synaptic max readily releasable pool */
#define SMAXRATE8  493	 /* maximum sustained synaptic release rate */
#define SGAIN8     494	 /* synaptic gain */
#define SVGAIN8    495	 /* synaptic vgain */
#define SDURH8     496	 /* synaptic high pass time const. (ms) */
#define SNFILTH8   497	 /* synaptic high pass nfilt */
#define SHGAIN8    498	 /* synaptic high pass gain */
#define SHOFFS8    499	 /* synaptic high pass offset */
#define SVSIZ8     500	 /* synaptic vesicle size */
#define SCOND8     501	 /* synaptic conductance */
#define SCMUL8     502	 /* synaptic conductance mult for region */
#define SCGRAD8    503	 /* synaptic conductance gradient from soma */
#define SEGRAD8    504	 /* synaptic conductance expon grad fr soma */
#define STHRESH8   505	 /* synaptic threshold */
#define SVNOISE8   506	 /* 1 -> vesicle noise, override, vnoise=0 */
#define SCOV8      507	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR8      508	 /* synaptic event time const. (ms) */
#define SFALL8     509	 /* synaptic event fall time const. (ms) */
#define SNFILT8    510	 /* synaptic vesicle nfilt */
#define STRCONC8   511	 /* synaptic transmitter concentration. */
#define SRESP8     512	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA8      513	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX8   514	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM8     515	 /* synaptic postsyn Ca pump Km. */
#define SVFBG8     516	 /* synaptic voltage feedback gain */
#define SVFBO8     517	 /* synaptic voltage feedback offset */
#define SCNFILT8   518	 /* second mesng. nfilt */
#define SCDUR8     519	 /* second mesng. time const. (ms) */
#define SCGAIN8    520	 /* synaptic second messenger gain */
#define SCOFF8     521	 /* synaptic second messenger offset */
#define SCNOISE8   522	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN8    523	 /* number of channels */
#define SUNIT8     524	 /* synaptic channel unitary conductace */
#define SVREV8     525	 /* synaptic reversal potential */

#define CELPOST9   526	 /* cell type to connect to (neg, no conn) */
#define CONPOST9   527	 /* connection number for postsyn cell */
#define CELDIV9    528	 /* number of postsyn cells to connect to */
#define GROWPRE9   529	 /* grow when making conn to postsyn cell */
#define SYNREG9    530	 /* synaptic region in presyn dendritic tree */
#define SYNREGP9   531	 /* synaptic region in postsyn dendritic tree */
#define SYNSPAC9   532	 /* synaptic spacing in presyn dendritic tree */
#define SYNANNI9   533	 /* inner rad of annulus in presyn dendr tree */
#define SYNANNO9   534	 /* outer rad of annulus in presyn dendr tree */
#define SYNANPI9   535	 /* inner rad of annulus in postsyn dend tree */
#define SYNANPO9   536	 /* outer rad of annulus in postsyn dend tree */
#define SYNANG9    537	 /* angle for presynaptic node rel to soma */
#define SYNRNG9    538	 /* range of presyn dendritic angles */
#define USEDYAD9   539	 /* synapse is dyad using preexisting type */
#define DYADTYP9   540	 /* type of dyad synapse to connect with */
#define AUTAPSE9   541	 /* synapse back to presynaptic node */
#define SYNNUM9    542	 /* number of synapses per connection */
#define SENSCA9    543	 /* synaptic release sensitivity calcium */
#define CAEGAIN9   544	 /* synaptic release calcium exponential sensitivity */
#define SRRPOOL9   545	 /* synaptic readily releasable pool */
#define SRRPOOLG9  546	 /* synaptic readily releasable pool gain */
#define SMRRPOOL9  547	 /* synaptic max readily releasable pool */
#define SMAXRATE9  548	 /* maximum sustained synaptic release rate */
#define SGAIN9     549	 /* synaptic gain */
#define SVGAIN9    550	 /* synaptic vgain */
#define SDURH9     551	 /* synaptic high pass time const. (ms) */
#define SNFILTH9   552	 /* synaptic high pass nfilt */
#define SHGAIN9    553	 /* synaptic high pass gain */
#define SHOFFS9    554	 /* synaptic high pass offset */
#define SVSIZ9     555	 /* synaptic vesicle size */
#define SCOND9     556	 /* synaptic conductance */
#define SCMUL9     557	 /* synaptic conductance mult for region */
#define SCGRAD9    558	 /* synaptic conductance gradient from soma */
#define SEGRAD9    559	 /* synaptic conductance expon grad fr soma */
#define STHRESH9   560	 /* synaptic threshold */
#define SVNOISE9   561	 /* 1 -> vesicle noise, override, vnoise=0 */
#define SCOV9      562	 /* 1=Poisson, <1->more regular, gamma dist */
#define SDUR9      563	 /* synaptic event time const. (ms) */
#define SFALL9     564	 /* synaptic event fall time const. (ms) */
#define SNFILT9    565	 /* synaptic vesicle nfilt */
#define STRCONC9   566	 /* synaptic transmitter concentration. */
#define SRESP9     567	 /* synaptic response (ampa,gaba,gj,etc. */
#define SPCA9      568	 /* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX9   569	 /* synaptic postsyn Ca pump vmax. */
#define SCAKM9     570	 /* synaptic postsyn Ca pump Km. */
#define SVFBG9     571	 /* synaptic voltage feedback gain */
#define SVFBO9     572	 /* synaptic voltage feedback offset */
#define SCNFILT9   573	 /* second mesng. nfilt */
#define SCDUR9     574	 /* second mesng. time const. (ms) */
#define SCGAIN9    575	 /* synaptic second messenger gain */
#define SCOFF9     576	 /* synaptic second messenger offset */
#define SCNOISE9   577	 /* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN9    578	 /* number of channels */
#define SUNIT9     579	 /* synaptic channel unitary conductace */
#define SVREV9     580	 /* synaptic reversal potential */

#define NPARAMS    581	 /* number of neuron parameters */

#define CELPRE       0	/* cell type to connect to (neg, no conn) */
#define CONPRE       1	/* connection number of presyn cell */
#define CELCONV      2	/* number of presyn cells to connect to */
#define GROWPOST     3	/* grow when making conn from presyn cell */
#define NCONNP       4	/* number of connection parameters */

#define CELPOST      0	/* cell type to connect to (neg, no conn) */
#define CONPOST      1	/* connection number for postsyn cell */
#define CELDIV       2	/* number of postsyn cells to connect to */
#define GROWPRE      3	/* grow when making conn to postsyn cell */
#define SYNREG       4	/* synaptic region in presyn dendritic tree */
#define SYNREGP      5	/* synaptic region in postsyn dendritic tree */
#define SYNSPAC      6	/* synaptic spacing in presyn dendritic tree */
#define SYNANNI      7	/* inner dia of annulus in presyn dendr tree */
#define SYNANNO      8	/* outer dia of annulus in presyn dendr tree */
#define SYNANPI      9	/* inner dia of annulus in postsyn dend tree */
#define SYNANPO     10	/* outer dia of annulus in postsyn dend tree */
#define SYNANG      11	/* angle for presynaptic node rel to soma */
#define SYNRNG      12	/* range of presyn dendritic angles */
#define USEDYAD     13	/* synapse is dyad using preexisting type */
#define DYADTYP     14	/* type of dyad synapse to connect with */
#define AUTAPSE     15	/* synapse back to presynaptic node */
#define SYNNUM      16	/* number of synapses per connection */
#define SENSCA      17	/* synaptic release sensitivity calcium */
#define CAEGAIN     18	/* synaptic release calcium exponential sensitivity */
#define SRRPOOL     19	/* synaptic readily releasable pool */
#define SRRPOOLG    20	/* synaptic readily releasable pool gain */
#define SMRRPOOL    21	/* synaptic max readily releasable pool */
#define SMAXRATE    22	/* maximum sustained synaptic release rate */
#define SGAIN       23	/* synaptic gain */
#define SVGAIN      24	/* synaptic vgain */
#define SDURH       25	/* synaptic high pass time const. (ms) */
#define SNFILTH     26	/* synaptic high pass nfilt */
#define SHGAIN      27	/* synaptic high pass gain */
#define SHOFFS      28	/* synaptic high pass offset */
#define SVSIZ       29	/* synaptic vesicle size */
#define SCOND       30	/* synaptic conductance */
#define SCMUL       31	/* synaptic conductance mult for region */
#define SCGRAD      32	/* synaptic conductance gradient from soma */
#define SEGRAD      33	/* synaptic conductance expon grad fr soma */
#define STHRESH     34	/* synaptic threshold */
#define SVNOISE     35	/* 1 -> vesicle noise, override,vnoise=0 */
#define SCOV        36	/* 1=Poisson, <1->more regular, gamma dist */
#define SDUR        37	/* synaptic event time const. (ms) */
#define SFALL       38	/* synaptic event fall time const. (ms) */
#define SNFILT      39	/* synaptic vesicle nfilt */
#define STRCONC     40	/* synaptic transmitter concentration. */
#define SRESP       41	/* synaptic response (ampa,gaba,gj,etc. */
#define SPCA        42	/* synaptic postsyn Ca perm (ampa,nmda,etc. */
#define SCAVMAX     43	/* synaptic postsyn Ca pump vmax. */
#define SCAKM       44	/* synaptic postsyn Ca pump Km. */
#define SVFBG       45	/* synaptic voltage feedback gain */
#define SVFBO       46	/* synaptic voltage feedback offset */
#define SCNFILT     47	/* second mesng. nfilt */
#define SCDUR       48	/* second mesng. time const. (ms) */
#define SCGAIN      49	/* synaptic second messenger gain */
#define SCOFF       50	/* synaptic second messenger offset */
#define SCNOISE     51	/* 1 -> channel noise, override,cnoise=0 */
#define SNCHAN      52	/* number of channels */
#define SUNIT       53	/* synaptic channel unitary conductace */
#define SVREV       54	/* synaptic reversal potential */
#define NSYNP       55	/* number of synaptic parameters */

#define NCONNI      10	/* number of input connection cell types */
#define NCONNO       9	/* number of output connection cell types  */

#define XGLUT        1	/* generic glutamate response */
#define XAMPA        2	/* AMPA (type 1) synaptic response */
#define XAMPA1       3	/* AMPA type 1 synaptic response */
#define XAMPA2       4	/* AMPA type 2 synaptic response */
#define XAMPA3       5	/* AMPA type 3 synaptic response */
#define XAMPA4       6	/* AMPA type 4 synaptic response */
#define XAMPA5       7	/* AMPA type 5 synaptic response */
#define XACH1        8	/* nicotinic ACh synaptic response */
#define XNMDA        9	/* NMDA type 1 synaptic response */
#define XNMDA2      10	/* NMDA type 2 synaptic response */
#define XKAINATE    11	/* Kainate synaptic response */
#define XMGLUR6     12	/* mGluR6 synaptic response */
#define XGABA       13	/* GABA type 1 synaptic response */
#define XGABA1      14	/* GABA type 1 synaptic response */
#define XGABA2      15	/* GABA type 2 synaptic response */
#define XGABA3      16	/* GABA type 3 synaptic response */
#define XGABA4      17	/* GABA type 4 synaptic response */
#define XGLY        18	/* Glycine synaptic response */
#define XGAPJ       19	/* gap junction synaptic response */
#define XDYAD       20	/* dyad synapse (uses other resp type) */
#define XVFB        21	/* voltage feedback (ephaptic or pH) */
#define NRESPTYPES  22	/* number of synaptic types */
