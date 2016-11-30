#!/usr/bin/env python
####################
# Example driver for maxlike masses in Clusters project
####################

import maxlike_controller
import maxlike_bentstep_voigt
import maxlike_masses
import astropytable_filehandler as atf


#######################

makeTestingController = lambda: maxlike_controller.Controller(modelbuilder=maxlike_masses.LensingModel(),
                                                              filehandler=atf.AstropyTableFilehandler(),
                                                              runmethod=maxlike_masses.ScanModelToFile())

makeController = lambda: maxlike_controller.Controller(modelbuilder=maxlike_bentstep_voigt.BentVoigtShapedistro(),
                                                       filehandler=atf.AstropyTableFilehandler(),
                                                       runmethod=maxlike_masses.SampleModelToFile())

controller = makeController()


if __name__ == '__main__':

    controller.run_all()

