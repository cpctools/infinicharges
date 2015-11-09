## Chose what kind of fit to perform.
fit_kind="ESP+TDF"

#coor="zif8.pdb"
coor="zif8.xyz"


## ESP control section
esp_from_cp2k=True
## first gamma value
g_start=1.0
## increase step
g_step=0.2
## number of gammas
n_gammas=3

## TDF control section 
#CP2K TDF data
tdf_from_cp2k=False
tdf_signal_fname = "signal.tdf"
tdf_coor_fname = "pos.xyz"

## If weight_tdf is uncommented only one loop will be performed for the
## value defined and the nw_loops is overridden
#weight_tdf=0.49
