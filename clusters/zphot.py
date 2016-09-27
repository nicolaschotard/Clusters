"""Photometric redshift analysis. Includes a wrapper to LEPHARE."""

import os
import sys
import numpy as N
import pylab as P
import seaborn
import subprocess

from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord


class LEPHARE(object):

    """Wrapper to the LEPHARE photometric redshift code."""

    def __init__(self, mags, magserr, cname, input_name=None, filters=None,
                 zpara=None, RA=None, DEC=None, ID=None):
        """
        Run the LEPHARE progam (zphota).

        param list mags: list of magnitude
        param list magserr:
        param string cname:
        param string filters: filter list
        param string zpara: default is $LEPHAREDIR/config/zphot_megacam.para
        param list RA: list of ra
        param list DEC: list of dec
        param list ID: list of ID
        """
        self.mags = mags
        self.magserr = magserr
        self.cluster_name = cname
        self.filters = filters
        self.config = os.environ["LEPHAREDIR"] + "/config/zphot_megacam.para" \
                      if zpara is None else zpara
        self.RA, self.DEC, self.ID = RA, DEC, ID

        if input_name is not None:
            self.input = input_name
            self.output = self.input.split(".")[0] + "_zphot.out"
        else:
            self.input = cname + "_zphot.in"
            self.output = cname + "_zphot.out"
        self.allinput = self.input.replace('.in', '.all')

        # initialize lephare output variables
        self.lephare_out = None
        self.data_out = None

        self.write_input()

    def write_input(self):
        """
        Create and write files needed to run LEPHARE.

        - the input data file for LEPHARE
        - a similare file containing the soruces ID along with their RA DEC.
        """
        f = open(self.input, 'w')
        if self.filters is not None:
            f.write("# id " + " ".join(["mag_%s" % filt for filt in self.filters]) +
                    " " + " ".join(["err_mag_%s" % filt for filt in self.filters]) +
                    "\n")
        for i, mags in enumerate(N.concatenate([self.mags, self.magserr]).T):
            f.write("%i %s\n" % (i, " ".join(["%.3f" % m for m in mags])))
        f.close()
        print "INFO: Input data saved in", self.input

        if self.RA is not None:
            f = open(self.allinput, 'w')
            if self.filters is not None:
                f.write("# id ID RA DEC " +
                        " ".join(["mag_%s" % filt for filt in self.filters]) + " " +
                        " ".join(["err_mag_%s" % filt for filt in self.filters]) +
                        "\n")
            for i, mags in enumerate(N.concatenate([self.mags, self.magserr]).T):
                f.write("%i %i %f %f %s\n" % (i, self.ID[i], self.RA[i], self.DEC[i],
                                              " ".join(["%.3f" % m for m in mags])))
            f.close()
            print "INFO: All data saved in", self.allinput

    def check_config(self, config=None):
        """
        Check that the SED and filters requested for the LePhare run do exist.

        If not: explains where the problem is and aborts.
        """
        if config is not None:
            if os.path.exists(config):
                self.config = config
            else:
                raise ValueError("%s does not exist" % config)
        with open(self.config) as f:
            for line in f:
                if "ZPHOTLIB" in line:
                    libs = line.split()[1].split(",")
        path_to_lib = os.environ["LEPHAREWORK"] + "/lib_mag/"

        for lib in libs:
            lib_tmp = path_to_lib + lib + ".bin"
            if not os.path.exists(lib_tmp):
                print "\nERROR: Requested library %s does not exist " % lib_tmp
                print "INFO: Available SED libraries are:\n ",
                for f in os.listdir(path_to_lib):
                    if f.endswith(".bin"):
                        print f
                print "--> Correct the ZPHOTLIB variable in %s or generate the missing LePhare SED libraries" % self.config
                sys.exit()

        counter = 0
        for lib in libs:
            lib_tmp = path_to_lib + lib + ".doc"
            with open(lib_tmp) as f:
                for line in f:
                    if "FILTER_FILE" in line:
                        filt_tmp = line.split()[1]
                        if counter == 0:
                            filt_ref = filt_tmp
                        if filt_tmp != filt_ref:
                            print "\nERROR: Requested SED libraries %s have been built using different filters. Either change the requested libraries or re-generate them accordingly." % libs
                            sys.exit()
            counter += 1

        path_to_filt = os.environ["LEPHAREWORK"] + "/filt/"
        if not os.path.exists(path_to_filt + filt_ref):
            print "\nERROR: The FILTER_FILE %s used by the SED libraries %s does not exists." % \
                ((path_to_filt + filt_ref), libs)
            sys.exit()

    def run(self, config=None):
        """
        Run LEPHARE.

        Default config file is $LEPHAREDIR/config/zphot_megacam.para.
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

        self.data_out = LEPHARO(self.output, all_input=self.allinput)


class LEPHARO(object):

    """Read LEPHARE output file."""

    def __init__(self, zphot_output, all_input=None):
        """Read the LEPHARe progam Output (zphota output)."""
        self.output = zphot_output
        if all_input is not None:
            self.input = all_input
            self.read_input()
        self.read()

    def read(self):
        """Read the output."""
        f = open(self.output, 'r')
        self.header = [l for l in f if l.startswith('#')]
        f.close()
        self.data_array = N.loadtxt(self.output, unpack=True)
        self.variables = N.loadtxt(os.getenv('LEPHAREDIR') +
                                   "/config/zphot_output.para", dtype='string')
        self.data_dict = {v: a for v, a in zip(self.variables, self.data_array)}
        self.nsources = len(self.data_dict['Z_BEST'])

    def read_input(self):
        """Read the input."""
        data = N.loadtxt(self.input, unpack=True)
        f = open(self.input, 'r')
        l = f.readlines()[0]
        self.input_data = {k: d for k, d in zip(l[2:-1].split(), data)}
        f.close()

    def hist(self, param, minv=None, maxv=None, nbins=None, xlabel=None,
             title=None, zclust=None, figname=""):
        """Plot histograms."""
        pval = self.data_dict[param]
        filt = N.array([1]*len(pval), dtype='bool')
        if minv is not None:
            filt &= (pval >= minv)
        if maxv is not None:
            filt &= (pval <= maxv)
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

        fig.savefig(figname + "_" + xlabel + "_zphot_hist.png")

    def plot(self, px, py, **kwargs):
        """
        Plot x vs. y.

        Possible kwargs are:
        - minx: lower limit of the x axis
        - maxx: upper limit of the x axis
        - miny: lower limit of the y axis
        - maxy: upper limit of the y axis
        - xlabel: label of the x axis
        - ylabel: label of the y axis
        - title: title of the figure
        - figname: name of the figure to save
        """
        pvalx = self.data_dict[px]
        pvaly = self.data_dict[py]
        filt = N.array([1] * len(pvalx), dtype='bool')
        if 'minx' in kwargs and kwargs['minx'] is not None:
            filt &= (pvalx >= kwargs['minx'])
        if 'maxx' in kwargs and kwargs['maxx'] is not None:
            filt &= (pvalx <= kwargs['maxx'])
        if 'miny' in kwargs and kwargs['miny'] is not None:
            filt &= (pvaly >= kwargs['miny'])
        if 'maxy' in kwargs and kwargs['maxy'] is not None:
            filt &= (pvaly <= kwargs['maxy'])
        pvalx, pvaly = pvalx[filt], pvaly[filt]
        fig = P.figure()
        ax = fig.add_subplot(111)
        ax.scatter(pvalx, pvaly)
        if 'xlabel' in kwargs and kwargs['xlabel'] is not None:
            px = kwargs['xlabel']
        if 'ylabel' in kwargs and kwargs['ylabel'] is not None:
            py = kwargs['ylabel']
        ax.set_xlabel(px)
        ax.set_ylabel(py)
        if 'title' in kwargs and kwargs['title'] is not None:
            ax.set_title(kwargs['title'])

        if 'figname' in kwargs and kwargs['figname'] is not None:
            fig.savefig(kwargs['figname'] + "_%s_vs_%s_zphot.png" % (py, px))
        else:
            fig.savefig("%s_vs_%s_zphot.png" % (py, px))

    def plot_map(self, title=None, figname="", zmin=0, zmax=999):
        """Plot the redshift sky-map."""
        if not hasattr(self, 'input_data'):
            print "WARNING: No input data given. Cannot plot the redshift map."
            return

        ra, dec, redshift = self.input_data['RA'], self.input_data['DEC'], self.data_dict['Z_BEST']

        # redshift has to be >= 0
        filt = (redshift >= zmin) & (redshift < zmax)
        ra, dec, redshift = ra[filt], dec[filt], redshift[filt]

        fig = P.figure()
        ax = fig.add_subplot(111, xlabel='RA (deg)', ylabel='DEC (deg)')
        scat = ax.scatter(ra, dec, c=redshift, cmap=(P.cm.jet))
        cb = fig.colorbar(scat)
        cb.set_label('Photometric redshift')
        if title is not None:
            ax.set_title(title)
        ax.set_xlim(xmin=min(ra) - 0.001, xmax=max(ra) + 0.001)
        ax.set_ylim(ymin=min(dec) - 0.001, ymax=max(dec) + 0.001)
        fig.savefig(figname + "_redshift_map.png")


def dict_to_array(d, filters='ugriz'):
    """Transform a dictionnary into a list of arrays."""
    return N.array([N.array(d[f]) for f in filters])


class ZSPEC(object):
    """Compare spectroscopic and photometric redshifts."""

    def __init__(self, sfile, names, unit='deg'):
        """Read the input data file and organize the data.

        :param str sfile: File containing the spectroscopic redshift
        :param list names: A list of names for the columns. You must use 'ra' and 'dec'
         for the coordinates.
        :param str unit: Unit of the given coordinates ('deg', 'rad'). Should be understandable by
         astropy.coordinates.SkyCoord

        The spectroscopic redshift must be named 'zspec'.
        """
        self.sfile = sfile
        self.data = ascii.read(sfile, names=names)
        if 'ra' not in self.data.keys():
            raise IOError("RA coordinate must be called 'ra'")
        elif 'dec' not in self.data.keys():
            raise IOError("DEC coordinate must be called 'dec'")
        elif 'zspec' not in self.data.keys():
            raise IOError("Spectroscopic reshift must be 'zspec'")
        self.skycoords = SkyCoord(self.data['ra'], self.data['dec'], unit=unit)
        self.zphot = self.skycoords_phot = self.match = None

    def load_zphot(self, ra, dec, zphot, unit='deg'):
        """Load the photometric informations and match them to the spectro ones.

        :param list ra: List of RA coordinates
        :param list dec: List of DEC coordinates
        :param list zphot: List of photometric redshift
        :param list unit: List of RA coordinates

        All lists must have the same length.
        """
        assert len(ra) == len(dec) == len(zphot)
        self.zphot = Table([ra, dec, zphot], names=['ra', 'dec', 'zphot'])
        self.skycoords_phot = SkyCoord(ra, dec, unit=unit)
        self.match = self.skycoords.match_to_catalog_sky(self.skycoords_phot)

    def plot(self, cut=0.05):
        """Plot a sky-map of the matches.

        TODO:

         - convert degree to mas
         - apply quality cuts on photometric redshift measurments
         - bug when loading 0018.zcat
        """
        if self.match is None:
            raise IOError("ERROR: You must load the photometric data first (load_zphot).")

        zspec = self.data['zspec']
        zphot = self.zphot['zphot'][self.match[0]]
        sdist = N.array(self.match[1])

        filt = sdist < cut
        zspec, zphot, sdist = [x[filt] for x in [zspec, zphot, sdist]]

        fig = P.figure(figsize=(15, 8))

        # z-phot as a function of z-spec
        ax = fig.add_subplot(121, xlabel='Z-spec', ylabel='Z-phot')
        scat = ax.scatter(zspec, zphot, c=sdist, cmap=(P.cm.copper))
        cb = fig.colorbar(scat)
        cb.set_label('On-sky distance (deg)')
        ax.set_title("%i galaxies" % len(self.match[0]))

        # z-phot - zspec as a function of sources on-sky distances
        ax = fig.add_subplot(122, xlabel='On-sky distance', ylabel='(Z-phot - Z-spec)')
        scat = ax.scatter(sdist, zphot - zspec, color='k')
        ax.set_title("%i galaxies" % len(self.match[0]))

        P.show()
