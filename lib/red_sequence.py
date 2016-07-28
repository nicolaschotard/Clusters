from ToolBox import Optimizer
import numpy as N
import pylab as P
import seaborn

def red_sequence(mags, fref, f1, f2, clip=2, colorcut=1):
    m, c = mags[fref][(mags[f1]-mags[f2])>colorcut], (mags[f1]-mags[f2])[(mags[f1]-mags[f2])>colorcut]
    op=Optimizer.LinearRegression(m, c, clip=clip)
    fig, ax = P.subplots(ncols=1, figsize=(7, 5))
    ax.plot(N.linspace(min(mags[fref]),max(mags[fref]), 1000),
           N.linspace(min(mags[fref]),max(mags[fref]), 1000)*op.a+op.b, 'g', ls='--')
    ax.plot(N.linspace(min(mags[fref]),max(mags[fref]), 1000),
           N.linspace(min(mags[fref]),max(mags[fref]), 1000)*op.a+op.b+0.2, 'g')
    ax.plot(N.linspace(min(mags[fref]),max(mags[fref]), 1000),
           N.linspace(min(mags[fref]),max(mags[fref]), 1000)*op.a+op.b-0.2, 'g')
    ax.scatter(mags[fref], mags[f1]-mags[f2], s=1, color='k')
    ax.set_xlabel(fref)
    ax.set_ylabel('%s - %s' % (f1, f2))
    print op.a, op.b
    return op

def color_histo(mags):
    filt = (mags['g'] - mags['r']) > 1.2
    for i, f1 in enumerate('gri'):
        for j, f2 in enumerate('gri'):
            if i >= j:
                continue
            fig, ax = P.subplots(ncols=1, figsize=(12, 8))
            ax.hist((mags[f1]-mags[f2])[filt],
                    bins=100, label='%s - %s' % (f1, f2))
            ax.legend(loc='best')
    P.show()

def color_mag_plot(mags):
    filt = (mags['g'] - mags['r']) > 1.2
    for k, fref in enumerate('gri'):
        for i, f1 in enumerate('gri'):
            for j, f2 in enumerate('gri'):
                if i >= j:
                    continue
                fig, ax = P.subplots(ncols=1, figsize=(12, 8))
                ax.scatter(mags[fref][filt], (mags[f1]-mags[f2])[filt],
                           s=1, label='%s - %s' % (f1, f2))
                ax.set_xlabel(fref)
                ax.set_ylabel('%s - %s' % (f1, f2))
                ax.legend(loc='best')
    P.show()
