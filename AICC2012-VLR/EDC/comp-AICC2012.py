import numpy as np
import matplotlib.pyplot as mpl
from scipy import interpolate
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

f=interpolate.interp1d(AICC2012_depth,AICC2012_age, bounds_error=False, fill_value=np.nan)
g=interpolate.interp1d(AICC2012_depth,AICC2012_age_sigma, bounds_error=False, fill_value=np.nan)
mpl.figure('Ice Age')
mpl.fill_between(f(IC_depth),-g(IC_depth),+g(IC_depth), color='0.8')
mpl.plot(f(IC_depth),IC_age-f(IC_depth), color='k', label='IceChrono-Datice')
mpl.plot(f(IC_depth),IC_age-f(IC_depth)-IC_age_sigma, color='k', linestyle='--', label='IceChrono credible interval')
mpl.plot(f(IC_depth),IC_age-f(IC_depth)+IC_age_sigma, color='k', linestyle='--')
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


f=interpolate.interp1d(AICC2012_depth,AICC2012_gage, bounds_error=False, fill_value=np.nan)
g=interpolate.interp1d(AICC2012_depth,AICC2012_gage_sigma, bounds_error=False, fill_value=np.nan)
mpl.figure('Air Age')
mpl.fill_between(f(IC_depth),-g(IC_depth),+g(IC_depth), color='0.8')
mpl.plot(f(IC_depth),IC_gage-f(IC_depth), color='k', label='IceChrono-Datice')
mpl.plot(f(IC_depth),IC_gage-f(IC_depth)-IC_gage_sigma, color='k', linestyle='--', label='IceChrono credible interval')
mpl.plot(f(IC_depth),IC_gage-f(IC_depth)+IC_gage_sigma, color='k', linestyle='--')
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


#f=interpolate.interp1d(AICC2012_depth,AICC2012_age, bounds_error=False, fill_value=np.nan)
#mpl.figure('Ice Age')
#mpl.plot(IC_age-f(IC_depth),IC_depth)
#mpl.xlabel('IceChrono - AICC2012 ice age difference (yr)')
#mpl.ylabel('depth (m)') 
#x1,x2,y1,y2 = mpl.axis()
#mpl.axis((x1, x2, y2, y1))
#pp=PdfPages('AICC2012-age.pdf')
#pp.savefig(mpl.figure('Ice Age'))
#pp.close()


#f=interpolate.interp1d(AICC2012_depth,AICC2012_age_sigma, bounds_error=False, fill_value=np.nan)
#mpl.figure('Ice Age credible interval')
#mpl.plot(IC_age_sigma-f(IC_depth),IC_depth)
#mpl.xlabel('IceChrono - AICC2012 ice age credible interval difference (yr)')
#mpl.ylabel('depth (m)') 
#x1,x2,y1,y2 = mpl.axis()
#mpl.axis((x1, x2, y2, y1))
#pp=PdfPages('AICC2012-age_sigma.pdf')
#pp.savefig(mpl.figure('Ice Age credible interval'))
#pp.close()


#g=interpolate.interp1d(AICC2012_depth,AICC2012_gage, bounds_error=False, fill_value=np.nan)
#mpl.figure('Air Age')
#mpl.plot(IC_gage-g(IC_depth),IC_depth)
#mpl.xlabel('IceChrono - AICC2012 air age difference (yr)')
#mpl.ylabel('depth (m)') 
#x1,x2,y1,y2 = mpl.axis()
#mpl.axis((x1, x2, y2, y1))
#pp=PdfPages('AICC2012-gage.pdf')
#pp.savefig(mpl.figure('Air Age'))
#pp.close()


#g=interpolate.interp1d(AICC2012_depth,AICC2012_gage_sigma, bounds_error=False, fill_value=np.nan)
#mpl.figure('Air Age credible interval')
#mpl.plot(IC_gage_sigma-g(IC_depth),IC_depth)
#mpl.xlabel('IceChrono - AICC2012 air age credible interval difference (yr)')
#mpl.ylabel('depth (m)') 
#x1,x2,y1,y2 = mpl.axis()
#mpl.axis((x1, x2, y2, y1))
#pp=PdfPages('AICC2012-gage_sigma.pdf')
#pp.savefig(mpl.figure('Air Age credible interval'))
#pp.close()



mpl.show()
