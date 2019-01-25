list_drillings=['EDC','VK','TALDICE','EDML','NGRIP']

opt_method='leastsq'  #leastsq, leastsq-parallel, none
nb_nodes=6         #Number of nodes for the leastsq-parallel mode

#Defines the colors for the graphs
color_obs='r'       #color for the observations
color_opt='k'       #color for the posterior scenario
color_mod='b'       #color for the prior scenario
color_ci='0.8'      #color for the confidence intervals
color_sigma='m'     #color for the uncertainty
color_di='g'        #color for the dated intervals
show_initial=False  #always put to False for now
color_init='c'      #always put to 'c' for now
scale_ageci=10.     #scaling of the confidence interval in the ice and air age figures
show_figures=False  #whether to show or not the figures at the end of the run
show_airlayerthick=False #whether to show the air layer thickness figure (buggy on anaconda)
