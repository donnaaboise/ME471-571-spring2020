import numpy
import matplotlib
from matplotlib import pylab, mlab, pyplot
np = numpy
plt = pyplot

from IPython.core.pylabtools import figsize, getfigs

from pylab import *
from numpy import *

from ipywidgets import interactive, fixed
import warnings



def plot_g(fig,g,g_line,it_line,a,kmax,x,x0,alpha):
    
    m = a/(1-alpha)
    g_line.set_ydata(g(x,m))

    tol = 1e-15
    xk = x0
    itlist_x = zeros((2*kmax,1))
    itlist_y = zeros((2*kmax,1))
    for k in range(0,kmax):
        xkp1 = g(xk,m)
        # print("{:3d} {:24.16f} {:12.4e}".format(k,xkp1,abs(xk-xkp1)))
        if abs(xkp1-xk) < tol:
            itlist_x[itlist_x == 0] = numpy.nan
            itlist_y[itlist_y == 0] = numpy.nan
            break

        itlist_x[2*k]   = xk
        itlist_y[2*k]   = xkp1
        itlist_x[2*k+1] = xkp1
        itlist_y[2*k+1] = xkp1
        xk = xkp1
            
    it_line.set_xdata(itlist_x)
    it_line.set_ydata(itlist_y)
    it_line.set_visible(True)

    # plot(itlist_x,itlist_y,'m.-',markersize=10)
    print("{:>22s} {:8d}".format('Number of iterations',k))
    print("{:>22s} {:8.4f}".format("alpha = g'(x)=1-a/m",1-a/m))
    
    fig.canvas.draw()
    return None


def fixed_point(xl,yl,a,b,kmax,x0,plot_f):

    warnings.filterwarnings("ignore",category=FutureWarning)

    fig = figure()
    x = linspace(xl[0],xl[1],3)  
    
    # Plot axes and line y=x
    plot(xl,0*xl,'k')
    plot(0*yl,yl,'k')

    x = linspace(xl[0],xl[1],3);
    plot(x,x,'k',linewidth=2)
    
    def f(x,a,b):
        return b - a*x

    def g(x,m):
        return (1-a/m)*x + b/m

    m0 = a
    g_line, = plot(x,g(x,m0),'g',linewidth=2)
    if plot_f:
        plot(x,f(x,a,b),'b',linewidth=2)
    
    xroot = b/a
    if plot_f:
        plot(xroot,0,'k.',markersize=15)
    plot(xroot,xroot,'k.',markersize=15)
    plot([xroot,xroot],yl,'k--')

    itlist_x = zeros((100,1))
    itlist_y = zeros((100,1))
    it_line, = plot(itlist_x,itlist_y,'m.-',markersize=10)
    it_line.set_visible(False)
        
    xlim(xl)
    ylim(yl)
    gca().set_aspect('equal')

    # alpha = 1-a/m;  -1 < alpha < 1
    alpha_lo = -2
    alpha_hi = 2
    alpha_inc = (alpha_hi - alpha_lo)/100
    
    # needed to create slider
    w = interactive(plot_g, fig=fixed(fig), g = fixed(g), 
                        g_line=fixed(g_line), it_line=fixed(it_line),
                        a = fixed(a),kmax=fixed(kmax),x=fixed(x), x0=fixed(x0),
                    alpha = (alpha_lo,alpha_hi,alpha_inc))
    w.background_color='white'
    return w