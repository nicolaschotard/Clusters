import os
import numpy as N
import pylab as P
import seaborn
import subprocess

class LEPHARE:

    def __init__(self, mags, magserr, cname, input=None, filters=None):
        """
        Runs the LEPHARE progam (zphota)
        param list mags:
        param list magserr:
        param string cname:
        """
        self.mags = mags
        self.magserr = magserr
        self.cluster_name = cname
        self.filters = filters
        self.config = '$LEPHAREDIR/config/zphot_megacam.para'
        
        if input is not None:
            self.input = input
            self.output = input.split(".")[0]+"_zphot.out"
        else:
            self.input = cname + "_zphot.in"
            self.output = cname + "_zphot.out"

        self.write_input()

    def write_input(self):
        f = open(self.input, 'w')
        if self.filters is not None:
            f.write("# id " + " ".join(["mag_%s" % filt for filt in self.filters]) + \
                    " " + " ".join(["err_mag_%s" % filt for filt in self.filters]) + "\n")
        for i, mags in enumerate(N.concatenate([self.mags, self.magserr]).T):
            f.write("%i %s\n" % (i, " ".join(["%.3f" % m for m in mags])))
        f.close()
        print "INFO: Input data saved in", self.input
        
    def run(self, config=None):
        """
        Default config file is $LEPHAREDIR/config/zphot_megacam.para
        Can be overwritten with the config argument
        """
        if config is not None:
            if os.path.exists(config):
                self.config = config
            else:
                raise ValueError("%s does not exist" % config)

        # build command line
        cmd = "zphota"
        cmd += " -c " + self.config
        cmd += " -CAT_IN " + self.input
        cmd += " -CAT_OUT " + self.output
        print "INFO: Will run '%s'" % cmd
        self.lephare_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
        print "INFO: LEPHARE output summary (full output in self.lephare_out)" 
        print "\n".join(["   " + zo for zo in self.lephare_out.split("\n")[-6:]])

        self.data_out = LEPHARO(self.output)

    def read_output(self):
        self.header = [l for l in f if l.startswith('#')]
        self.outdata = N.loadtxt(self.output, unpack=True)

class LEPHARO:
    def __init__(self, zphot_output):
        """
        Read the LEPHARe progam Output (zphota)
        """
        self.output = zphot_output
        self.read()

    def read(self):
        f = open(self.output, 'r')
        self.header = [l for l in f if l.startswith('#')]
        f.close()
        self.data_array = N.loadtxt(self.output, unpack=True)
        self.variables = N.loadtxt(os.getenv('LEPHAREDIR')+"/config/zphot_output.para", dtype='string')
        self.data_dict = {v: a for v, a in zip(self.variables, self.data_array)}
        self.nsources = len(self.data_dict['Z_BEST'])

    def hist(self, param, min=None, max=None, nbins=None, xlabel=None, title=None, zclust=None, figname=""):
        pval = self.data_dict[param]
        filt = N.array([1]*len(pval), dtype='bool')
        if min is not None:
            filt &= (pval >= min)
        if max is not None:
            filt &= (pval <= max)
        pval = pval[filt]
        fig = P.figure()
        ax = fig.add_subplot(111, ylabel='#')
        ax.hist(pval, bins=nbins if nbins is not None else 10)
        if xlabel is None:
            xlabel = param
        ax.set_xlabel(xlabel)
        if title is not None:
            ax.set_title(title)
        if zclust is not None:
            ax.axvline(zclust, color='r', label='Cluster redshift (%.4f)' % zclust)
            ax.legend(loc='best')

        fig.savefig(figname+"_"+xlabel+"_zphot_hist.png")

    def plot(self, px, py, minx=None, maxx=None, miny=None, maxy=None,
             xlabel=None, ylabel=None, title=None, figname=""):
        pvalx = self.data_dict[px]
        pvaly = self.data_dict[py]
        filt = N.array([1]*len(pvalx), dtype='bool')
        if minx is not None:
            filt &= (pvalx >= minx)
        if maxx is not None:
            filt &= (pvalx <= maxx)
        if miny is not None:
            filt &= (pvaly >= miny)
        if maxy is not None:
            filt &= (pvaly <= maxy)
        pvalx, pvaly = pvalx[filt], pvaly[filt]
        fig = P.figure()
        ax = fig.add_subplot(111)
        ax.scatter(pvalx, pvaly)
        if xlabel is None:
            xlabel = px
        if ylabel is None:
            ylabel = py        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title is not None:
            ax.set_title(title)

        fig.savefig(figname+"_%s_vs_%s_zphot.png" % (ylabel, xlabel))
        
def dict_to_array(d, filters='ugriz'):
    return N.array([N.array(d[f]) for f in filters])
