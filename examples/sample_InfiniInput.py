## Chose what kind of fit to perform.
fit_kind="ESP+TDF"
#fit_kind="ESP"
#fit_kind="RESP"
#fit_kind="RG"

#coor="zif8.pdb"
coor="zif8.xyz"


## ESP control section
esp_from_cp2k=True
## first gamma value
g_start=0.6
## increase step
g_step=0.1
## number of gammas
n_gammas=20

## TDF control section (for ESP only fit this will be skipped)
#CP2K TDF data
tdf_from_cp2k=False
tdf_signal_fname = "tdf.tdf"
tdf_coor_fname = "pos.xyz"

## If weight_tdf is uncommented only one loop will be performed for the
## value defined and the nw_loops is overridden
#weight_tdf=0.49

##For Ewald calibration. Default: locut=1.0e-10
#locut=1.0e-30
######################################################################

#########RESTRAINTS SECTION###########################################
#RESP
restraint_weight=0.33
#restraint_weight=[
#("C1", 7.0),
#("C2", 10.0656),
#("C3", 10.0533),
#("H2", 10.1573),
#("H3", 10.0710),
#("H4", 10.0),
#("N" , 0.0),
#("Zn", 100000.0)
#]

restraints=[
("C2",-0.0656),
("C1", 0.2996),
("C3",-0.0533),
("H2", 0.1573),
("H3", 0.0710),
("H4", 0.0696),
("N" ,-0.7403),
("Zn", 1.6786)
]

#restraints=[
#("C1",  -0.5276514 ),
#("C2",  -0.5276514 ),
#("C3",  -0.5276514 ),
#("H2",  -0.34873393),
#("H3",  -0.34873393),
#("H4",  -0.34873393),
#("N" ,  -0.58665038),
#("Zn",  -0.59579913)]

#RG (corresponding to the to Goddard-type weight of REPEAT)
rg_weight=0.1

