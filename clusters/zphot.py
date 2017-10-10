"""Photometric redshift analysis. Includes a wrapper to LEPHARE and BPZ.

- LEPHARE: http://www.cfht.hawaii.edu/~arnouts/LEPHARE/lephare.html
- BPZ: http://www.stsci.edu/~dcoe/BPZ
"""

from __future__ import print_function
import os
import sys
import subprocess
import numpy as N
import pylab as P
from scipy.optimize import curve_fit
from astropy.io import ascii
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from . import utils as cutils


class LEPHARE(object):

    """Wrapper to the LEPHARE photometric redshift code.

    http://www.cfht.hawaii.edu/~arnouts/LEPHARE/lephare.html
    """

    def __init__(self, magnitudes, errors, zpara=None, spectro_file=None, **kwargs):
        """
        Run the LEPHARE progam (zphota).

        :param list magnitudes: Magnitudes. A list of list
        :param list errors: Error on magnitudes. Same shape as the magnitude list.
        :param string zpara: default is $LEPHAREDIR/config/zphot_megacam.para
        :param dictionnary files: A dictioannry conating the names of the

        Kwargs include the following list of possible parameters

        :param string basename: Base name for all created files
        :param string cname: Name of the studied cluster. If basename if not given,
         it will be used as base name for all created file
        :param string filters: filter list
        :param list ra: list of ra
        :param list dec: list of dec
        :param list id: list of ID
        """
        self.data = {'mag': magnitudes, 'err': errors}
        self.kwargs = kwargs
        self.config = os.environ["LEPHAREDIR"] + \
            "/config/zphot_megacam.para" if zpara is None else zpara

        self.spectro_file = spectro_file
        # Name of created files?
        prefix = ""
        if 'basename' in kwargs:
            prefix = kwargs['basename'] + "_"
        elif 'cname' in kwargs:
            prefix = kwargs['cname'] + "_"
        self.files = {}
        self.files['input'] = prefix + "zphot.in"
        self.files['output'] = prefix + "zphot.out"
        self.files['pdz_output'] = prefix + "zphot"
        self.files['all_input'] = self.files['input'].replace('.in', '.all')

        # Initialize lephare output variables
        self.lephare_out = None
        self.data_out = None

        self.write_input()

    def write_input(self):
        """
        Create and write files needed to run LEPHARE.

        - the input data file for LEPHARE
        - a similar file containing the sources ID along with their RA DEC.
        """
        f = open(self.files['input'], 'w')
        if self.spectro_file is None:
            # No spectroscopic redshift file provided in config.yaml
            # --> only need the SHORT format for LePhare input file
            if 'filters' in self.kwargs:
                f.write("# id " + " ".join(["mag_%s" % filt for filt in self.kwargs['filters']]) +
                        " " + " ".join(["err_mag_%s" % filt for filt in self.kwargs['filters']]) +
                        "\n")
                for i, mags in enumerate(N.concatenate([self.data['mag'], self.data['err']]).T):
                    f.write("%i %s\n" %
                            (i, " ".join(["%.3f" % m for m in mags])))

                f.close()
                print("INFO: Input data saved in", self.files['input'])
        else:
            # Spectroscopic redshift file provided in config.yaml
            # --> Need to write LePhare input file in the LONG format,
            # i.e with 'context' and 'spectroz' of matching galaxies
            zspec = ZSPEC(self.spectro_file, names=[
                          'object', 'ra', 'dec', 'zspec'])
            zspec.skycoords = SkyCoord(
                zspec.data['ra'], zspec.data['dec'], unit='deg')
            skycoords_cat = SkyCoord(
                self.kwargs['ra'], self.kwargs['dec'], unit='deg')
            idx, d2d, d3d = skycoords_cat.match_to_catalog_sky(zspec.skycoords)
            zp = zspec.data['zspec'][idx]
            # identify galaxies with bad match, i.e. dist > 300 mas
            bad = N.where(d2d.mas > 300)
            zp[bad] = -99
            print("INFO: Using " + str(len(idx) - N.size(bad)) +
                  " galaxies for spectroz training")
            if 'filters' in self.kwargs:
                f.write("# id " + " ".join(["mag_%s" % filt for filt in self.kwargs['filters']]) +
                        " " + " ".join(["err_mag_%s" % filt for filt in self.kwargs['filters']]) +
                        " context" + " zspec" + "\n")
                # tells LePhare to run using the u, g, r, i and z bands.
                context = 31
                for i, mags in enumerate(N.concatenate([self.data['mag'], self.data['err']]).T):
                    f.write("%i %s %s\n" % (i, " ".join(["%.3f" % m for m in mags]),
                                            " ".join(("%i" % context, "%.3f" % zp[i]))))
                f.close()
        if 'ra' in self.kwargs:
            f = open(self.files['all_input'], 'w')
            if 'filters' in self.kwargs is not None:
                f.write("# id ID RA DEC " +
                        " ".join(["mag_%s" % filt for filt in self.kwargs['filters']]) + " " +
                        " ".join(["err_mag_%s" % filt for filt in self.kwargs['filters']]) +
                        "\n")
            for i, mags in enumerate(N.concatenate([self.data['mag'], self.data['err']]).T):
                f.write("%i %i %f %f %s\n" % (i, self.kwargs['id'][i],
                                              self.kwargs['ra'][i], self.kwargs['dec'][i],
                                              " ".join(["%.3f" % m for m in mags])))
            f.close()
            print("INFO: All data saved in", self.files['all_input'])

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
                print("\nERROR: Requested library %s does not exist " % lib_tmp)
                print("INFO: Available SED libraries are:\n ")
                for f in os.listdir(path_to_lib):
                    if f.endswith(".bin"):
                        print(f)
                        print("--> Correct the ZPHOTLIB variable in %s" % self.config +
                              "or generate the missing LePhare SED libraries")
                sys.exit()

        for counter, lib in enumerate(libs):
            lib_tmp = path_to_lib + lib + ".doc"
            with open(lib_tmp) as f:
                for line in f:
                    if "FILTER_FILE" in line:
                        filt_tmp = line.split()[1]
                        if counter == 0:
                            filt_ref = filt_tmp
                        if filt_tmp != filt_ref:
                            print("\nERROR: Requested SED libraries %s" % libs +
                                  " have been built using different filters. Either change the " +
                                  "requested libraries or re-generate them accordingly.")
                            sys.exit()

        path_to_filt = os.environ["LEPHAREWORK"] + "/filt/"
        if not os.path.exists(path_to_filt + filt_ref):
            print("\nERROR: The FILTER_FILE %s used by the SED libraries %s does not exists." %
                  (path_to_filt + filt_ref, libs))
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
        cmd += " -CAT_IN " + self.files['input']
        cmd += " -CAT_OUT " + self.files['output']
        cmd += " -PDZ_OUT " + self.files['pdz_output']
        print("INFO: Will run '%s'" % cmd)

        self.lephare_out = subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, shell=True)
        print("INFO: LEPHARE output summary (full output in self.lephare_out)")
        print(
            "\n".join(["   " + zo for zo in self.lephare_out.split("\n")[-6:]]))

        self.data_out = ZPHOTO(self.files['output'], self.files['pdz_output'], zcode_name='lephare',
                               all_input=self.files['all_input'], **self.kwargs)


class BPZ(object):

    """
    Wrapper to the BPZ photometric redshift code.

    http://www.stsci.edu/~dcoe/BPZ
    """

    def __init__(self, magnitudes, errors, zpara=None, spectro_file=None, **kwargs):
        """
        Run the BPZ progam (zphota).

        :param list magnitudes: Magnitudes. A list of list
        :param list errors: Error on magnitudes. Same shape as the magnitude list.

        Kwargs include the following list of possible parameters

        :param string basename: Base name for all created files
        :param string cname: Name of the studied cluster. If basename if not given,
                             it will be used as base name for all created file
        :param string filters: filter list
        :param list ra: list of ra
        :param list dec: list of dec
        :param list id: list of ID
        """
        self.data = {'mag': magnitudes, 'err': errors}
        self.kwargs = kwargs
        self.config = zpara
        self.spectro_file = spectro_file

        # Name of created files?
        prefix = ""
        if 'basename' in kwargs:
            prefix = kwargs['basename'] + "_"
        elif 'cname' in kwargs:
            prefix = kwargs['cname'] + "_"
        self.files = {}
        self.files['input'] = prefix + "bpz.in"
        self.files['flux_comparison'] = prefix + "bpz.flux_comparison"
        self.files['output'] = prefix + "bpz.bpz"
        self.files['pdz_output'] = prefix + "bpz.probs"
        self.files['columns'] = prefix + "bpz.columns"
        self.files['all_input'] = self.files['input'].replace('.in', '.all')

        # Initialize BPZ output variables
        self.bpz_out = None
        self.data_out = None

        self.write_input()
        self.build_columns_file()

    def write_input(self):
        """
        Create and write files needed to run BPZ.

        - the input data file for BPZ
        - a similar file containing the sources ID along with their RA DEC.
        """
        f = open(self.files['input'], 'w')
        if 'filters' in self.kwargs:
            f.write("# id " + " ".join(["mag_%s" % filt for filt in self.kwargs['filters']]) +
                    " " + " ".join(["err_mag_%s" % filt for filt in self.kwargs['filters']]) +
                    "\n")
            for i, mags in enumerate(N.concatenate([self.data['mag'], self.data['err']]).T):
                f.write("%i %s\n" % (i, " ".join(["%.3f" % m for m in mags])))

            f.close()
            print("INFO: Input data saved in", self.files['input'])
        if 'ra' in self.kwargs:
            f = open(self.files['all_input'], 'w')
            if 'filters' in self.kwargs is not None:
                f.write("# id ID RA DEC " +
                        " ".join(["mag_%s" % filt for filt in self.kwargs['filters']]) + " " +
                        " ".join(["err_mag_%s" % filt for filt in self.kwargs['filters']]) +
                        "\n")
            for i, mags in enumerate(N.concatenate([self.data['mag'], self.data['err']]).T):
                f.write("%i %i %f %f %s\n" % (i, self.kwargs['id'][i],
                                              self.kwargs['ra'][i], self.kwargs['dec'][i],
                                              " ".join(["%.3f" % m for m in mags])))
            f.close()
            print("INFO: All data saved in", self.files['all_input'])

    def build_columns_file(self, prefix='CFHT_megacam_', sufix='p',
                           filters=None, ref='i', z_s=False):
        """Build and write the 'columns' file.

        Hardcoded for test purpose.
        """
        if filters is None:  # set the default values as the CFHT ones
            filters = ['u', 'g', 'r', 'i', 'z']
        f = open(self.files['columns'], 'w')
        f.write("# Filter  columns  AB/Vega  zp_error  zp_offset\n")
        for i, filt in enumerate(filters):
            f.write("%s%s%s     %i, %i   AB        0.01      0.00\n" %
                    (prefix, filt, sufix, i + 2, i + len(filters) + 2))
        if z_s:
            f.write("Z_S                  %i\n" % (len(filters) * 2 + 2))
        f.write("M_0                 %i\n"
                    % [j + 2 for j, filt in enumerate(filters) if filt == ref][0])
        f.write("ID                    1\n")
        f.close()

    def run(self):
        """
        Run BPZ.

        Configuration file must exist in the current directory.

        .. todo:: Build the configuration file on the fly (the .columns)
        """
        if not os.getenv('BPZPATH'):
            mess = "BPZPATH if not defined. Please install BPZ "
            mess += "(see http://www.stsci.edu/~dcoe/BPZ/install.html)."
            raise IOError(mess)

        if not os.path.exists(self.files['columns']):
            raise IOError("%s does not exist" % self.files['columns'])

        # build command line options from param file
        if self.config is not None:
            opt_arr = N.genfromtxt(self.config, dtype=None)
            option = ['-' + opt_arr[i, 0] + ' ' + opt_arr[i, 1]
                      for i in N.arange(len(opt_arr))]
            options = ' '.join(e for e in option)
        else:
            options = ''

        # build command line
        cmd = "python $BPZPATH/bpz.py %s " % self.files['input'] + options
        print("INFO: Will run '%s'" % cmd)

        self.bpz_out = subprocess.check_output(
            cmd, stderr=subprocess.STDOUT, shell=True)
        print("INFO: BPZ output summary (full output in self.bpz_out)")
        print("\n".join(["   " + zo for zo in self.bpz_out.split("\n")[:20]]))

        self.data_out = ZPHOTO(self.files['output'], self.files['pdz_output'],
                               zcode_name='bpz', all_input=self.files['all_input'],
                               **self.kwargs)


class ZPHOTO(object):

    """Read photoz code (LePhare, BPZ) output file and creates/saves astropy tables."""

    def __init__(self, zphot_output, zphot_pdz_output, zcode_name=None, all_input=None, **kwargs):
        """Read the photoz progam (LePhare, BPZ, ...) output."""
        self.files = {}
        self.files['output'] = zphot_output
        self.files['pdz_output'] = zphot_pdz_output
        self.kwargs = kwargs
        if all_input is not None:
            self.files['input'] = all_input
            self.read_input()
        self.code = zcode_name
        self.read()

    def read(self):
        """Read the output."""
        f = open(self.files['output'], 'r')
        self.data_array = N.loadtxt(self.files['output'], unpack=True)

        if self.code == 'lephare':
            self.header = [l for l in f if l.startswith('#')]
            f.close()
            self.variables = N.loadtxt(os.getenv('LEPHAREDIR') +
                                       "/config/zphot_output.para", dtype='string')
            self.data_dict = {v: a for v, a in zip(
                self.variables, self.data_array)}
            self.nsources = len(self.data_dict['Z_BEST'])
            self.pdz_zbins = N.loadtxt(
                self.files['pdz_output'] + '.zph', unpack=True)
            self.pdz_val = N.loadtxt(
                self.files['pdz_output'] + '.pdz', unpack=True)

        elif self.code == 'bpz':
            self.header = [l for l in f if l.startswith('##')]
            f.seek(0)
            self.variables = [l[4:].replace(' ', '').split(
                '\n')[0] for l in f if l.startswith('# ')]
            f.seek(0)
            # BPZ does not provide a zbins file.
            # Needs to create it from zmin, zmax and dz specified in output file
            zmin = float([line.split('=')[1]
                          for line in f if 'ZMIN' in line][0])
            f.seek(0)
            zmax = float([line.split('=')[1]
                          for line in f if 'ZMAX' in line][0])
            f.seek(0)
            dz = float([line.split('=')[1] for line in f if 'DZ' in line][0])
            f.close()

            self.data_dict = {v: a for v, a in zip(
                self.variables, self.data_array)}
            self.nsources = len(self.data_dict['Z_B'])
            self.pdz_zbins = N.arange(zmin, zmax + dz, dz)
            self.pdz_val = N.loadtxt(self.files['pdz_output'], unpack=True,
                                     usecols=N.arange(1, len(self.pdz_zbins) + 1))

    def save_zphot(self, file_out, path_output, overwrite=False):
        """Save the output of photoz code (z_best, chi^2, pdz) into astropy table."""
        # Duplicates the zbins vector for each object.
        # It is redundant information but astropy tables need each field to have the
        # same size. Or maybe I'm missing something.
        zbins = N.tile(self.pdz_zbins, (len(self.kwargs['id']), 1))

        # Converts LePhare or BPZ likelihood to actual probability density
        for i in N.arange(len(self.pdz_val.T)):
            norm = N.trapz(self.pdz_val[:, i], self.pdz_zbins)
            new_pdz_val = self.pdz_val[:, i] / norm
            self.pdz_val[:, i] = new_pdz_val

        # Creates astropy table to be saved in path_output of file_out
        new_tab = hstack([Table([self.kwargs['id']]), Table([self.kwargs['ra']]),
                          Table([self.kwargs['dec']]), Table(self.data_dict),
                          Table([zbins], names=['zbins']),
                          Table([self.pdz_val.T], names=['pdz'])],
                         join_type='inner')

        # Rename BPZ Z_B to Z_BEST to match LePhare
        if 'Z_B' in new_tab.keys():
            new_tab.rename_column('Z_B', 'Z_BEST')

        # overwrite keyword of data.write(file,path) does not only overwrites
        # the data in path, but the whole file, i.e. (we lose all other
        # paths in the process) --> overwrite_or_append (see above)

        cutils.overwrite_or_append(file_out, path_output, new_tab, overwrite=overwrite)

        print("INFO: ", self.code, "data saved in", file_out, "as", path_output)

    def read_input(self):
        """Read the input."""
        data = N.loadtxt(self.files['input'], unpack=True)
        f = open(self.files['input'], 'r')
        l = f.readlines()[0]
        self.input_data = {k: d for k, d in zip(l[2:-1].split(), data)}
        f.close()

    def hist(self, param, **kwargs):
        """Plot histograms.

        Possible kwargs

        :params float minv: Lower value of the histogram
        :params float maxv: Upper value of the histogram
        :params int nbins: Number of bins. Default is 10.
        :params string xlabel: An xlbal for the figure
        :params string title: A title for the figure
        :params float zclust: Redshift of the studies cluster
        """
        # Get the data and apply some cut if asked
        pval = self.data_dict[param]
        filt = N.array([1] * len(pval), dtype='bool')
        if 'minv' in kwargs:
            filt &= (pval >= kwargs['minv'])
        if 'maxv' in kwargs:
            filt &= (pval <= kwargs['maxv'])
        pval = pval[filt]

        # Plot the histogram
        fig = P.figure()
        ax = fig.add_subplot(111, ylabel='#')
        ax.hist(pval, bins=kwargs['nbins'] if 'nbins' in kwargs else 10)
        xlabel = kwargs['xlabel'] if 'xlabel' in kwargs else param
        ax.set_xlabel(xlabel)
        if 'title' in kwargs:
            ax.set_title(kwargs['title'])
        if 'zclust' in kwargs:
            ax.axvline(kwargs['zclust'], color='r',
                       label='Cluster redshift (%.4f)' % kwargs['zclust'])
            ax.legend(loc='best')

        # Save the figure
        fig.savefig(self.files['output'].replace(
            '.out', '') + "_" + xlabel + "_zphot_hist.png")

    def plot(self, px, py, **kwargs):
        """
        Plot x vs. y.

        Possible kwargs are:
        :params float minx: lower limit of the x axis
        :params float maxx: upper limit of the x axis
        :params float miny: lower limit of the y axis
        :params float maxy: upper limit of the y axis
        :params string xlabel: label of the x axis
        :params string ylabel: label of the y axis
        :params string title: title of the figure
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
            fig.savefig(self.files['output'].replace(
                '.out', '') + "_%s_vs_%s_zphot.png" % (py, px))
        else:
            fig.savefig("%s_vs_%s_zphot.png" % (py, px))

    def plot_map(self, title=None, zmin=0, zmax=999):
        """Plot the redshift sky-map."""
        if not hasattr(self, 'input_data'):
            print("WARNING: No input data given. Cannot plot the redshift map.")
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
        fig.savefig(self.files['output'].replace(
            '.out', '') + "_redshift_map.png")


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

        # Check whether duplicate galaxies exist in the spectroz sample
        # and remove them. Ideally would average out the various occurances.
        # (To be done later)
        ra = ['{:.20}'.format(x) for x in self.data['ra']]
        dec = ['{:.20}'.format(x) for x in self.data['dec']]
        radec = N.core.defchararray.add(ra, dec)
        unique_radec, good = N.unique(radec, return_index=True)
        if len(unique_radec) < len(radec):
            print("INFO: There are " + str(len(radec) - len(unique_radec)) +
                  " duplicate galaxies in spectroz sample. They are removed.")
        bad = N.delete(N.arange(len(self.data)), good)
        self.data.remove_rows(bad)

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
        idx, d2d, d3d = self.skycoords.match_to_catalog_sky(
            self.skycoords_phot)
        self.match = Table([idx, d2d.marcsec, d3d], names=['idx', 'd2d', 'd3d'],
                           meta={'description': ['Indices into first catalog.',
                                                 'On-sky separation between the '
                                                 'closest match for each element',
                                                 '3D distance between the closest'
                                                 'match for each element'],
                                 'unit': ['', 'marcsec', '']})

    def plot(self, cut=300, path_to_png=None):
        """Plot a sky-map of the matches."""
        if self.match is None:
            raise IOError(
                "ERROR: You must load the photometric data first (load_zphot).")

        zspec = self.data['zspec']
        zphot = self.zphot['zphot'][self.match['idx']]
        sdist = N.array(self.match['d2d'])
        filt = sdist < cut
        zspec, zphot, sdist = [x[filt] for x in [zspec, zphot, sdist]]

        fig = P.figure(figsize=(15, 8))

        # z-phot as a function of z-spec
        ax = fig.add_subplot(121, xlabel='Z-spec', ylabel='Z-phot')
        ax.set_xlim([0, 1.5])
        ax.set_ylim([0, 3.5])
        scat = ax.scatter(zspec, zphot, c=sdist, cmap=(P.cm.copper))
        ax.plot(N.arange(int(max(zspec + 1))),
                N.arange(int(max(zspec + 1))), c='red')
        cb = fig.colorbar(scat)
        cb.set_label('On-sky distance (marcsec)')
        ax.set_title("%i galaxies" % len(self.match[filt]))

        # z-phot - zspec as a function of sources on-sky distances
        ax = fig.add_subplot(122, xlabel='On-sky distance',
                             ylabel='(Z-phot - Z-spec)')
        scat = ax.scatter(sdist, zphot - zspec, color='k')
        ax.set_title("%i galaxies" % len(self.match[filt]))
        ax.set_xscale('log')
        ax.set_xlim([1., 1.e6])
        ax.set_ylim([-1., 3.5])

        if path_to_png is not None:
            fig.savefig(path_to_png)

        fig = P.figure()

        # radec scatter plot of all catalogue and zspec galaxies
        ax = fig.add_subplot(121, xlabel='ra', ylabel='dec')
        ax.scatter(self.skycoords_phot.ra, self.skycoords_phot.dec,
                   color='k', label='Photo-z', s=15)
        ax.scatter(self.skycoords.ra, self.skycoords.dec,
                   color='r', label='Spectro-z', s=12)

        # radec scatter plot of matched catalogue and zspec galaxies, within the cut criterion
        ax = fig.add_subplot(122, xlabel='ra', ylabel='dec')
        ax.scatter(self.skycoords_phot.ra[self.match['idx'][filt]],
                   self.skycoords_phot.dec[self.match['idx'][filt]],
                   color='k', label='Photo-z', s=5)
        ax.scatter(self.skycoords.ra[filt], self.skycoords.dec[filt],
                   color='r', label='Spectro-z', s=5)

        if path_to_png is not None:
            fig.savefig(path_to_png.replace('.png', '_map.png'))

        P.show()

    def scatter(self, zclust, cluster=None, cut=0.1, stability=False):
        """Redshift scatter in the cluster.

        Plot the spectroscopic redshift distribution and apply a gaussian fit.
        """
        # Apply selection
        zspec = self.data['zspec'][N.abs(self.data['zspec'] - zclust) < cut]

        # The gaussian fit
        hist, bin_edges = N.histogram(zspec, bins=len(zspec) / 4)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Use initial guess for the fitting coefficients (A, mu and sigma above)
        coeff = curve_fit(gauss, bin_centers, hist, p0=[
                          N.max(hist[0]), zclust, 0.3])[0]

        # Plot
        fig = P.figure()
        ax = fig.add_subplot(111, xlabel='Z-spec', ylabel='#',
                             title='' if cluster is None else cluster +
                             ', %i object included (cut=%.2f)' % (len(zspec), cut))
        ax.hist(zspec, bins=len(zspec) / 4, histtype='stepfilled')
        ax.axvline(zclust, color='r', label='Z-cluster (%.3f)' % zclust, lw=2)
        ax.plot(bin_centers, gauss(bin_centers, *coeff),
                label='Gaussian fit (mean, std)=(%.3f, %.3f)' %
                (coeff[1], N.abs(coeff[2])), lw=2, color='k')
        ax.legend(loc='best')

        if stability:
            cuts = [0.5, 0.1, 0.05, 0.03]
            z, ze = zip(*[self.scatter(zclust, cluster=cluster, cut=c)
                          for c in cuts])
            fig = P.figure()
            ax = fig.add_subplot(111, xlabel='Selection cut around cluste rredshift',
                                 ylabel='Estimated redshift',
                                 title='' if cluster is None else cluster +
                                 ', %i object included (cut=%.2f)' % (len(zspec), cut))
            ax.axhline(N.mean(z), label='Average value', color='r', lw=2)
            ax.errorbar(cuts, z, yerr=ze, label='Individual fits', color='k',
                        capsize=20, elinewidth=3)
            ax.legend(loc='best')
            print("INFO: Stability over the followed range of selection cuts:")
            print("       ", cuts)
            print("INFO: Input redshift: %.4f" % zclust)
            print("INFO: Average redshift: %.4f =/- %.4f" %
                  (N.mean(z), N.sqrt(N.std(z)**2 + N.mean(ze)**2)))
            P.show()
            return N.mean(z), N.sqrt(N.std(z)**2 + N.mean(ze)**2)
        P.show()
        return coeff[1], N.abs(coeff[2])


def gauss(x, *p):
    """Model function to be used to fit a gaussian distribution."""
    A, mu, sigma = p
    return A * N.exp(- (x - mu) ** 2 / (2. * sigma ** 2))
