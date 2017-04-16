import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.backends.backend_pdf import PdfPages

readarray=np.loadtxt('output.txt')

IC_depth=readarray[:,0]
IC_age=readarray[:,1]
IC_age_sigma=readarray[:,2]
IC_gage=readarray[:,3]
IC_gage_sigma=readarray[:,4]
IC_accu=readarray[:,5]
IC_accu_sigma=readarray[:,6]
IC_thinning=readarray[:,7]
IC_thinning_sigma=readarray[:,8]
IC_LID=readarray[:,9]
IC_LID_sigma=readarray[:,10]
IC_Ddepth=readarray[:,11]
IC_Ddepth_sigma=readarray[:,12]


readarray=np.loadtxt('AICC2012.txt')

AICC2012_depth=readarray[:,0]
AICC2012_age=readarray[:,1]
AICC2012_age_sigma=readarray[:,2]
AICC2012_gage=readarray[:,3]
AICC2012_gage_sigma=readarray[:,4]
AICC2012_accu=readarray[:,5]
AICC2012_thinning=readarray[:,6]
AICC2012_LID=readarray[:,7]/0.7

mpl.figure('Ice Age')
f=np.interp(IC_depth, AICC2012_depth, AICC2012_age, right=np.nan)
g=np.interp(IC_depth, AICC2012_depth,AICC2012_age_sigma, right=np.nan)
mpl.fill_between(f,-g,+g, color='0.8')
mpl.plot(f,IC_age-f, color='k', label='IceChrono-Datice')
mpl.plot(f,IC_age-f-IC_age_sigma, color='k', linestyle='--', label='IceChrono credible interval')
mpl.plot(f,IC_age-f+IC_age_sigma, color='k', linestyle='--')
mpl.ylabel('Ice age difference (yr)')
mpl.xlabel('AICC2012 age (yr b1950)') 
x1,x2,y1,y2 = mpl.axis()
mpl.legend()
mpl.axis((0, 800000, -6000, 6000))
pp=PdfPages('AICC2012-age.pdf')
pp.savefig(mpl.figure('Ice Age'))
pp.close()
mpl.axis((0, 60000, -3000, 3000))
pp=PdfPages('AICC2012-age-60kyr.pdf')
pp.savefig(mpl.figure('Ice Age'))
pp.close()


mpl.figure('Air Age')
f=np.interp(IC_depth, AICC2012_depth,AICC2012_gage, right=np.nan)
g=np.interp(IC_depth, AICC2012_depth,AICC2012_gage_sigma, right=np.nan)
mpl.fill_between(f,-g,+g, color='0.8')
mpl.plot(f,IC_gage-f, color='k', label='IceChrono-Datice')
mpl.plot(f,IC_gage-f-IC_gage_sigma, color='k', linestyle='--', label='IceChrono credible interval')
mpl.plot(f,IC_gage-f+IC_gage_sigma, color='k', linestyle='--')
mpl.ylabel('Air age difference (yr)')
mpl.xlabel('AICC2012 age (yr b1950)') 
x1,x2,y1,y2 = mpl.axis()
mpl.legend()
mpl.axis((0, 800000, -6000, 6000))
pp=PdfPages('AICC2012-gage.pdf')
pp.savefig(mpl.figure('Air Age'))
pp.close()
mpl.axis((0, 60000, -3000, 3000))
pp=PdfPages('AICC2012-gage-60kyr.pdf')
pp.savefig(mpl.figure('Air Age'))
pp.close()





mpl.show()
