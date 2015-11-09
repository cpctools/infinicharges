## Chose what kind of fit to perform.
fit_kind="RESP"

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

#RESP single weight for each atom type
restraint_weight=0.33

##RESP different weight for each atom type
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

#RESP target values
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

